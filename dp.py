from sys import exit
from argparse import ArgumentParser

from utility import (
    get_var,
    get_ref,
    get_seq_from_fastq,
    is_same_arm,
    load_cytobands,
    rev_comp,
    run_shell_cmd,
)

"""class DP is used for executing DP. It first generates ref. seq. based on
    the target SV candidates, and then executes external DP executables
"""

class DpOptions:
    def __init__(self, sv_str, query_name, query_strand, fastq_prefix, gap_buf_size=None, query_seq=None):
        self.sv_str = sv_str
        self.query_seq = None
        self.query_name = query_name
        self.query_strand = query_strand
        self.fastq_prefix = fastq_prefix
        self.gap_buf_size = gap_buf_size
        self.query_seq = query_seq

        self.preview_size = None


class DP:
    def __init__(self, options):
        self.load_config()

        self.options = options
        if options.gap_buf_size:
            self.gap_buf_size = options.gap_buf_size

        self.preview_size = None

        if options.sv_str:
            arr = options.sv_str.split('_')
            sv_type = arr[4]

            if sv_type in ['DUP']:
                self.bp = {
                    'chrom2': arr[0],
                    'end': int(arr[1]),
                    'chrom': arr[2],
                    'start': int(arr[3]),
                    'sv_type': arr[4]
                }
            elif sv_type in ['DEL', 'INV', 'TRA', 'CHIMERIC', 'REF']:
                self.bp = {
                    'chrom': arr[0],
                    'start': int(arr[1]),
                    'chrom2': arr[2],
                    'end': int(arr[3]),
                    'sv_type': arr[4]
                }
            elif sv_type in ['DEL-FROM-INS-FROM', 'DEL-TO-INS-FROM', 'DEL-FROM-INS-TO', 'DEL-TO-INS-TO']:
                self.bp = {
                    'chrom': arr[0],
                    'start': int(arr[1]),
                    'chrom2': arr[2],
                    'end': int(arr[3]),
                    'sv_type': arr[4]
                }
            else:
                print('Error: unsupported sv_type=%s. Program aborted' % (sv_type))
                exit(0)

        if options.query_seq:
            if self.preview_size:
                options.query_seq = (
                    options.query_seq[500:500+self.preview_size] +
                    options.query_seq[-(self.preview_size+500):-int(self.preview_size)]
                )
            self.buf_size = len(options.query_seq) + self.gap_buf_size

        load_cytobands()

    def load_config(self):
        self.dp_exe = get_var('dp', 'dp_exe')
        self.gap_buf_size = int(get_var('dp', 'gap_buf_size'))

    def swap_value(self, list, a, b):
        list[a], list[b] = list[b], list[a]

    def dp_pos_to_genomic_pos(self, dp_pos):
        zone, chrom, genomic_pos = -1, -1, -1

        for _idx, entry in enumerate(self.dp_to_genomic_map_table):
            if entry['start'] <= dp_pos < entry['end']:
                dir = 1 if entry['genomic_start'] < entry['genomic_end'] else -1
                genomic_pos = entry['genomic_start'] + (dp_pos - entry['start']) * dir
                zone = entry['zone']
                chrom = entry['chrom']
                break

        return [zone, chrom, genomic_pos]

    def get_gap_type(self, gap_start_zone, gap_end_zone):
        if self.bp['sv_type'] in ['DEL','DUP','TRA','DEL-FROM-INS-FROM','DEL-TO-INS-FROM','DEL-FROM-INS-TO','DEL-TO-INS-TO','CHIMERIC','REF']:
            return self.bp['sv_type']

        elif self.bp['sv_type'] == 'INV':
            if abs(gap_start_zone-gap_end_zone) == 1:
                return 'INV'
            else:
                return 'DEL'

        return ''

    def run_dp(self):
        self.gen_ref_seq()

        cmd = '%s' % (self.dp_exe)
        stdin_input = '%s %s' % (self.options.query_seq, self.ref_seq)
        output = run_shell_cmd(cmd, stdin_input)

        arr = [int(s) for s in output.split()]

        if len(arr) < 12:
            print('Error in run_dp. output', output)
            return {}

        score, query_start, query_end, ref_start, ref_end, gap_start, gap_end, query_pos, \
            genomic_ref_start, genomic_ref_end, genomic_gap_start, genomic_gap_end = arr[0:12]

        _ref_start_zone, _ref_start_chrom, genomic_ref_start = self.dp_pos_to_genomic_pos(ref_start)
        _ref_end_zone, _ref_end_chrom, genomic_ref_end = self.dp_pos_to_genomic_pos(ref_end)
        gap_start_zone, gap_start_chrom, genomic_gap_start = self.dp_pos_to_genomic_pos(gap_start)
        gap_end_zone, gap_end_chrom, genomic_gap_end = self.dp_pos_to_genomic_pos(gap_end)
        gap_type = self.get_gap_type(gap_start_zone, gap_end_zone)

        if self.bp['sv_type'] != gap_type:
            genomic_gap_start = -1
            genomic_gap_end = -1

        result = {
            'score': score,
            'query_start': query_start,
            'query_end': query_end,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'gap_start': gap_start,
            'gap_end': gap_end,
            'query_pos': query_pos,
            'genomic_ref_start': genomic_ref_start,
            'genomic_ref_end': genomic_ref_end,
            'genomic_gap_start': genomic_gap_start,
            'genomic_gap_end': genomic_gap_end,
            'query_len': len(self.options.query_seq),
            'ref_len': len(self.ref_seq),
            'chrom': gap_start_chrom,
            'chrom2': gap_end_chrom,
            'qname': self.options.query_name,
            'gap_type': gap_type,
            'sv_type': self.bp['sv_type'],
        }

        if self.bp['sv_type'] == 'DUP':
            self.swap_value(result, 'query_start', 'query_end')
            self.swap_value(result, 'ref_start', 'ref_end')
            self.swap_value(result, 'gap_start', 'gap_end')
            self.swap_value(result, 'genomic_ref_start', 'genomic_ref_end')
            self.swap_value(result, 'genomic_gap_start', 'genomic_gap_end')

        self.result = result

        return result

    def gen_ref_seq(self):
        options = self.options

        gap_buf_size = self.gap_buf_size
        buf_size = self.buf_size
        bp = self.bp

        bp_zone = []
        if bp['sv_type'] in ['DEL', 'DUP']:
            # need to check if this fix can be applied to other sv_type
            if bp['sv_type'] in ['DUP']:
                #gap_buf_size += len(options.query_seq)
                gap_buf_size += 5000

            bp_zone = [{}, {}]
            bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'] = bp['chrom'], bp['start']-buf_size, bp['start']+gap_buf_size
            bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'] = bp['chrom2'], bp['end']-gap_buf_size, bp['end']+buf_size

            #print('bp', bp)
            #print('bp_zone[0]', bp_zone[0])
            #print('bp_zone[1]', bp_zone[1])

            ref_seq = get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'])
            ref_seq += get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'])
        elif bp['sv_type'] in ['INV']:
            bp_zone = [{}, {}, {}]
            bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'] = bp['chrom'], bp['start']-gap_buf_size, bp['start']+buf_size
            bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'] = bp['chrom2'], bp['end']-buf_size, bp['end']+buf_size
            bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end'] = bp['chrom'], bp['start']-buf_size, bp['start']+gap_buf_size

            ref_seq = rev_comp(get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end']))
            ref_seq += get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'])
            ref_seq += rev_comp(get_ref(bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end']))

            bp_zone[0]['start'], bp_zone[0]['end'] = bp_zone[0]['end'], bp_zone[0]['start']
            bp_zone[2]['start'], bp_zone[2]['end'] = bp_zone[2]['end'], bp_zone[2]['start']
        elif bp['sv_type'] in ['TRA']:
            bp_zone = [{}, {}, {}]
            bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'] = bp['chrom'], bp['start']-buf_size, bp['start']+gap_buf_size
            bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'] = bp['chrom2'], bp['end']-buf_size, bp['end']+buf_size
            bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end'] = bp['chrom'], bp['start']-gap_buf_size, bp['start']+buf_size

            if is_same_arm(bp['chrom'], bp['start'], bp['chrom2'], bp['end']):
                print('same arm')

                ref_seq = get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'])
                ref_seq += get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'])
                ref_seq += get_ref(bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end'])
            else:
                print('different arm')

                ref_seq = get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'])
                ref_seq += rev_comp(get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end']))
                ref_seq += get_ref(bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end'])

                bp_zone[1]['start'], bp_zone[1]['end'] = bp_zone[1]['end'], bp_zone[1]['start']
        elif bp['sv_type'] in ['CHIMERIC']:
            bp_zone = [{}, {}, {}]
            bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'] = bp['chrom'], bp['start']-buf_size, bp['start']+gap_buf_size
            bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'] = bp['chrom2'], bp['end']-buf_size, bp['end']+buf_size
            bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end'] = bp['chrom'], bp['start']-gap_buf_size, bp['start']+buf_size

            if options.query_strand in ['++','--']:
                ref_seq = get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'])
                ref_seq += get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'])
                ref_seq += get_ref(bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end'])
            else:
                ref_seq = get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'])
                ref_seq += rev_comp(get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end']))
                ref_seq += get_ref(bp_zone[2]['chrom'], bp_zone[2]['start'], bp_zone[2]['end'])

                bp_zone[1]['start'], bp_zone[1]['end'] = bp_zone[1]['end'], bp_zone[1]['start']
        elif bp['sv_type'] in ['DEL-FROM-INS-FROM','DEL-TO-INS-FROM','DEL-FROM-INS-TO','DEL-TO-INS-TO']:
            bp_zone = [{}, {}]

            bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'] = bp['chrom'], bp['start']-buf_size, bp['start']+gap_buf_size
            bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'] = bp['chrom2'], bp['end']-gap_buf_size, bp['end']+buf_size

            """
            bp_zone[0]['start'], bp_zone[0]['end'] = bp['start']-100000, bp['start']+100000
            bp_zone[1]['start'], bp_zone[1]['end'] = bp['end']-100000, bp['end']+100000
            """

            if 'DEL-FROM' in bp['sv_type']:
                ref_seq = get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'])
            else:
                ref_seq = rev_comp(get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end']))
                bp_zone[0]['start'], bp_zone[0]['end'] = bp_zone[0]['end'], bp_zone[0]['start']

            if 'INS-FROM' in bp['sv_type']:
                ref_seq += get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end'])
            else:
                ref_seq += rev_comp(get_ref(bp_zone[1]['chrom'], bp_zone[1]['start'], bp_zone[1]['end']))
                bp_zone[1]['start'], bp_zone[1]['end'] = bp_zone[1]['end'], bp_zone[1]['start']

        elif bp['sv_type'] in ['REF']:
            bp_zone = [{}]
            bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'] = bp['chrom'], bp['start']-buf_size, bp['end']+buf_size

            ref_seq = get_ref(bp_zone[0]['chrom'], bp_zone[0]['start'], bp_zone[0]['end'])

        self.ref_seq = ref_seq

        map_table = []
        total_length = 0
        for i in range(len(bp_zone)):
            curr_length = abs(bp_zone[i]['end'] - bp_zone[i]['start'])
            map_table.append({
                'zone': i,
                'chrom': bp_zone[i]['chrom'],
                'start': total_length,
                'end': total_length + curr_length,
                'genomic_start': bp_zone[i]['start'],
                'genomic_end': bp_zone[i]['end'],
            })
            total_length += curr_length

        self.dp_to_genomic_map_table = map_table

    def run(self):
        if self.options.query_name and self.options.fastq_prefix and not self.options.query_seq:
            self.options.query_seq = get_seq_from_fastq(self.options.query_name, self.options.fastq_prefix)

            if self.preview_size:
                self.options.query_seq = (
                    self.options.query_seq[500:500+self.preview_size] +
                    self.options.query_seq[-(self.preview_size+500):-int(self.preview_size)]
                )

            if self.options.query_strand[0] == '-':
                self.options.query_seq = rev_comp(self.options.query_seq)

            self.buf_size = len(self.options.query_seq) + self.gap_buf_size

        result = self.run_dp()

        return result


if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-sv_str', '--sv_str', help='SV str [chr_start_chr2_end_svtype]', required=True)

    parser.add_argument('-query_seq', '--query_seq', help='query seq', required=False)

    parser.add_argument('-query_name', '--query_name', help='query name', required=False)
    parser.add_argument('-fastq_prefix', '--fastq_prefix', help='fastq prefix', required=False)
    parser.add_argument('-query_strand', '--query_strand', help='query strand', required=False)

    parser.add_argument('-preview_size', '--preview_size', help='preview size', required=False, type=int)

    options = parser.parse_args()
    dp = DP(options)
    result = dp.run()

    print('result', result)
