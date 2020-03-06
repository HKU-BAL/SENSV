import logging
import csv
from argparse import ArgumentParser
from multiprocessing import Pool

from dp import DpOptions, DP
from utility import (
    get_var,
    get_short_name,
    get_seq_from_fastq,
    get_normal_avg_depth,
    rev_comp,
    run_shell_cmd,
    init_logger,
)

"""
class ChainCaller finds target SV regions based on chain info.
"""


class ChainCallerOptions:
    def __init__(self, depth_file, chain_file, fastq_file, bam_file, output_prefix, min_sv_size_by_depth, disable_dp_filter, gender, nprocs):
        self.depth_file = depth_file
        self.chain_file = chain_file
        self.fastq_file = fastq_file
        self.bam_file = bam_file
        self.output_prefix = output_prefix
        self.min_sv_size_by_depth = min_sv_size_by_depth
        self.disable_dp_filter = disable_dp_filter
        self.gender = gender
        self.nprocs = nprocs


class ChainCaller:

    def __init__(self, options):
        self.load_config()
        self.options = options

        self.depth_region = []
        self.target_chain_list = []

        self.dp_file = f'{options.output_prefix}.dp'
        self.bed2_file = f'{options.output_prefix}.dp.bed2'

    def load_config(self):
        self.config = {
            'buf_size': int(get_var('chain_caller', 'buf_size')),
            'max_chain_count': int(get_var('chain_caller', 'max_chain_count')),
            'gap_buf_size': int(get_var('dp', 'gap_buf_size')),
        }

    def load_depth_region(self, do_filter=0):
        options = self.options

        with open(options.depth_file) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                arr = line.strip().split()
                if arr[0] == '23':
                    arr[0] = 'X'
                elif arr[0] == '24':
                    arr[0] = 'Y'

                """
                region_type = 'DEL'
                if len(arr) >= 4:
                    region_type = arr[3]
                """
                region_type = arr[4]

                chrom, start, end = arr[0], int(float(arr[1])), int(float(arr[2]))
                sv_str = f'{chrom}_{start}_{end}'

                if do_filter:
                    avg_depth = self.get_avg_depth(sv_str)
                    normal_avg_depth = self.get_avg_depth(sv_str, 1)

                    is_skip = 0
                    if normal_avg_depth < 0.5:
                        is_skip = 1
                    elif options.gender == 'f' or chrom not in ['X', 'Y']:
                        if options.gender == 'f' and chrom == 'Y':
                            is_skip = 1

                        # if end - start >= options.min_sv_size_by_depth and avg_depth < 1:
                        #    is_skip = 1
                    """
                    else:

                        if end - start >= options.min_sv_size_by_depth and avg_depth > 1:
                            is_skip = 1
                    """

                    if is_skip:
                        print(f'gender = {options.gender}, avg depth = {avg_depth}, depth region skipped: {sv_str}')
                        continue

                self.depth_region.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'region_type': region_type,
                })

        print('loaded depth_region count', len(self.depth_region))

    def get_bp_depth_region_info(self, info):
        config = self.config
        options = self.options

        bp_depth_region_info = {'depth_region': {}}

        for region in self.depth_region:
            chrom, start, end, region_type = region['chrom'], region['start'], region['end'], region['region_type']

            start1, start2 = start - config['buf_size'], start + config['buf_size']
            end1, end2 = end - config['buf_size'], end + config['buf_size']

            in_depth_region = False

            if 'strand' in info and 'strand2' in info:
                if (
                    region_type in ['DEL', 'DUP'] and
                    info['strand'] != info['strand2']
                ) or (
                    region_type in ['INV'] and
                    info['strand'] == info['strand2']
                ):
                    continue

            if (
                info['bp_chrom'] == chrom and
                info['bp_chrom2'] == chrom and
                start1 <= info['bp_start'] <= start2 and
                end1 <= info['bp_end'] <= end2
            ):
                if region_type in ['DEL', 'DUP'] and info['bp_end'] - info['bp_start'] < options.min_sv_size_by_depth:
                    continue

                in_depth_region = True

            if in_depth_region:
                bp_depth_region_info['depth_region'] = {
                    'depth_chrom': chrom,
                    'depth_start': start,
                    'depth_end': end,
                    'region_type': region_type
                }

        return bp_depth_region_info

    def find_target_chain_list(self):
        config = self.config
        options = self.options

        with open(options.chain_file) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                chain = []
                prev_qname = None

                for l in [line]:
                    arr = l.strip().split()
                    qname, strand, n_chain, seq_len = arr[0], arr[1], int(arr[2]), int(arr[3])

                    # the two lines must be referring to the same reads
                    if prev_qname and qname != prev_qname:
                        print(f'Error: reads in chain files not appear in '+' and'-' pair! ({prev_qname},{qname})')
                        continue

                    prev_qname = qname

                    count = 0
                    for i in range(n_chain):
                        if count >= config['max_chain_count']:
                            break

                        idx = 4 + i * 8
                        _chain_id, seed_count, query_start, query_end, ref_chrom, ref_start, ref_end, score = \
                            int(arr[idx]), int(arr[idx+1]), int(arr[idx+2]), int(arr[idx+3]), arr[idx+4], \
                            int(arr[idx+5]), int(arr[idx+6]), int(arr[idx+7])

                        if not options.disable_dp_filter:
                            if 1.0 * score / (query_end - query_start) < 0.3:
                                # if 1.0 * score / (query_end - query_start) < 0.25:
                                continue

                        # if seed_count < 200:
                        #    continue

                        chain.append({
                            'qname': get_short_name(qname),
                            'strand': strand,
                            'query_start': query_start,
                            'query_end': query_end,
                            'ref_chrom': ref_chrom,
                            'ref_start': ref_start,
                            'ref_end': ref_end,
                            'seq_len': seq_len,
                            'score': score,
                            'seed_count': seed_count
                        })

                        count += 1

                for i in range(len(chain)-1):
                    for j in range(i+1, len(chain)):
                        if chain[i]['ref_end'] < chain[j]['ref_start']:
                            chain1, chain2 = chain[i], chain[j]
                        else:
                            chain1, chain2 = chain[j], chain[i]

                        info = {
                            'qname': qname,
                            'strand': chain1['strand'],
                            'strand2': chain2['strand'],
                            'bp_chrom': chain1['ref_chrom'],
                            'bp_start': chain1['ref_end'],
                            'bp_chrom2': chain2['ref_chrom'],
                            'bp_end': chain2['ref_start'],
                            'query_start': chain1['query_start'],
                            'query_end': chain2['query_end'],
                            'seq_len': chain2['seq_len'],
                        }

                        bp_depth_region_info = self.get_bp_depth_region_info(info)
                        have_target_gap = 0
                        need_dp = 0

                        if bp_depth_region_info['depth_region']:
                            depth_region = bp_depth_region_info['depth_region']
                            depth_start = depth_region['depth_start']
                            depth_end = depth_region['depth_end']
                            region_type = depth_region['region_type']

                            have_target_gap = 1

                            # debug
                            #options.disable_dp_filter = True

                            if region_type == 'DEL':
                                if options.disable_dp_filter or (
                                    # chain1['query_start'] < chain2['query_start'] and
                                    # chain1['query_end'] <= chain2['query_end'] and
                                    chain1['query_end'] < chain2['query_end'] and

                                    chain1['query_start'] < 150 and
                                    chain2['seq_len'] - chain2['query_end'] < 300 and

                                        chain1['query_end'] + 500 >= chain2['query_start']):

                                    #print('chain1', chain1)
                                    #print('chain2', chain2)

                                    need_dp = 1
                            elif region_type == 'DUP':
                                if options.disable_dp_filter or (
                                        chain2['query_start'] <= chain1['query_start'] and
                                        chain2['query_end'] < chain1['query_end']):
                                    need_dp = 1
                            elif region_type == 'INV':
                                # Todo: add checking
                                # if options.disable_dp_filter or (chain1['query_end'] + 500 >= chain2['query_start'] \
                                #     and chain1['query_end'] < chain2['query_end']):
                                if True:
                                    need_dp = 1
                            else:
                                need_dp = 1

                            self.target_chain_list.append({
                                'have_target_gap': have_target_gap,
                                'need_dp': need_dp,
                                'qname': qname,
                                'strand': strand,
                                'bp_chrom': info['bp_chrom'],
                                'bp_start': info['bp_start'],
                                'bp_chrom2': info['bp_chrom2'],
                                'bp_end': info['bp_end'],
                                'sv_type': region_type,
                                'depth_start': depth_start,
                                'depth_end': depth_end,
                                # 'query_start':info['query_start'],
                                # 'query_end':info['query_end'],
                                # 'score':info['score']}),
                                'query_start': info['query_start'],
                                'query_end': info['query_end'],
                                'seq_len': info['seq_len'],
                            })

        """
        self.target_chain_list.sort(key=lambda x: x['seq_len'], reverse=True)

        # debug
        for chain_list in self.target_chain_list:
            print('%s: %d %s %d %s %d %s' % (chain_list['qname'], chain_list['seq_len'], \
                chain_list['bp_chrom'], chain_list['bp_start'], chain_list['bp_chrom2'], \
                chain_list['bp_end'], chain_list['sv_type']))
        """

    def __call__(self, chain_list):
        options = self.options
        result = {}

        if chain_list['need_dp'] != 1:
            return result

        query_seq = get_seq_from_fastq(chain_list['qname'], options.fastq_file, '')

        if chain_list['sv_type'] in ['DEL', 'DUP']:
            if chain_list['strand'] == '1':
                query_seq = rev_comp(query_seq)
        elif chain_list['sv_type'] in ['INV']:
            if chain_list['strand'] == '0':
                query_seq = rev_comp(query_seq)

        sv_str = f"{chain_list['bp_chrom']}_{chain_list['bp_start']}_{chain_list['bp_chrom2']}_{chain_list['bp_end']}_{chain_list['sv_type']}"

        return DP(
            DpOptions(
                sv_str,
                chain_list['qname'],
                None,
                None,
                None,
                query_seq
            )
        ).run()

    def do_dp(self):
        config = self.config
        options = self.options

        nprocs = self.options.nprocs
        pool = Pool(processes=nprocs)
        results = pool.map(self, self.target_chain_list)
        pool.close()
        pool.join()

        # remove empty entries
        #results = filter(None, results)
        results = filter(lambda x: x and x['chrom'] != -1, results)

        toIntStr = lambda text: text if text.isdigit() else '23' if text == 'X' else '24'
        sorted_results = sorted(
            results,
            key=lambda k: (
                toIntStr(k['chrom']).rjust(2) +
                str(k['genomic_gap_start']).rjust(9) +
                str(k['genomic_gap_end']).rjust(9)
            )
        )

        called_sv = {}
        with open(self.dp_file, 'w') as f_dp, open(self.bed2_file, 'w') as f_bed2:
            dp_writer = csv.writer(f_dp, delimiter='\t')
            bed2_writer = csv.writer(f_bed2, delimiter='\t')

            dp_header=[
                'chrom',
                'genomic_gap_start',
                'chrom2',
                'genomic_gap_end',
                'qname',
                'query_len',
                'ref_len',
                'score',
                'sv_type'
            ]
            bed_header = [
                'chrom',
                'genomic_gap_start',
                'chrom2',
                'genomic_gap_end',
                'sv_type',
                'qname',
                'query_len',
                'ref_len',
                'score',
                'depth_start',
                'depth_end'
            ]

            dp_writer.writerow(dp_header)

            for result in sorted_results:
                if not result:
                    continue

                output = [str(result[i]) for i in dp_header]
                dp_writer.writerow(output)

                # if result['score'] < result['query_len'] * 1.1:
                if result['score'] < min(result['query_len'], config['gap_buf_size']*2) * 1.1:
                    continue

                info = {
                    'bp_chrom': result['chrom'],
                    'bp_start': result['genomic_gap_start'],
                    'bp_chrom2': result['chrom2'],
                    'bp_end': result['genomic_gap_end']
                }
                bp_depth_region_info = self.get_bp_depth_region_info(info)

                if bp_depth_region_info['depth_region']:
                    depth_region = bp_depth_region_info['depth_region']
                    region_type = depth_region['region_type']

                    result['sv_type'] = region_type
                    result['depth_start'] = depth_region['depth_start']
                    result['depth_end'] = depth_region['depth_end']

                    if region_type in ['DEL', 'DUP', 'INV']:
                        if result['genomic_gap_end'] - result['genomic_gap_start'] < options.min_sv_size_by_depth:
                            continue

                    bed_output = [str(result[i]) for i in bed_header]
                    bed_key = '_'.join(bed_output)

                    if bed_key in called_sv:
                        continue

                    called_sv[bed_key] = 1

                    bed2_writer.writerow(bed_output)

    # debug
    def print_chain_list(self):
        header = ['qname', 'bp_chrom', 'bp_start', 'bp_chrom2', 'bp_end']

        for chain in self.target_chain_list:
            if not chain['need_dp']:
                continue

            output = [str(chain[h]) for h in header]
            print('\t'.join(output))

    def get_avg_depth(self, sv_str, use_normal_male_bam=0):
        options = self.options

        if use_normal_male_bam:
            #bam_file = config['normal_male_bam']
            return get_normal_avg_depth(sv_str)['normal_avg_depth']
        else:
            bam_file = options.bam_file

        arr = sv_str.split('_')
        chrom, start, end = arr[0], int(arr[1]), int(arr[2])

        # cmd = "samtools depth -a %s -r %s:%d-%d | awk '{ sum += $3 } END { if (NR > 0) print sum / NR }'" % \
        #      (options.bam_file, chrom, start, end)

        max_depth = 8
        cmd = "samtools depth -a -a %s -r %s:%d-%d | awk '{ max=%d; depth = ($3<max?$3:max); sum += depth } END { if (NR > 0) print sum / NR }'" % \
              (bam_file, chrom, start, end, max_depth)

        #print('cmd', cmd)

        output = run_shell_cmd(cmd).split('\n')

        return float(output[0])

    def run(self):
        self.load_depth_region()
        self.find_target_chain_list()

        have_target_gap = 0
        need_dp = 0
        for chain_list in self.target_chain_list:
            if chain_list['have_target_gap'] == 1:
                have_target_gap += 1

            if chain_list['need_dp'] == 1:
                need_dp += 1
                """
                logging.info('qname: %s, strand: %s, bp_chrom: %s, bp_start: %d, bp_end: %d' %
                             (chain_list['qname'], chain_list['strand'], chain_list['bp_chrom'],
                              chain_list['bp_start'], chain_list['bp_end']))
                """

        logging.info('have_target_gap: %d, need_dp: %d' % (have_target_gap, need_dp))

        # debug
        # self.print_chain_list()

        self.do_dp()


if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-depth_file', '--depth_file', help='depth file', required=True)
    parser.add_argument('-chain_file', '--chain_file', help='chain file', required=True)
    parser.add_argument('-fastq_file', '--fastq_file', help='fastq file', required=True)
    parser.add_argument('-bam_file', '--bam_file', help='bam file', required=True)
    parser.add_argument('-output_prefix', '--output_prefix', help='output prefix', required=True)
    parser.add_argument('-min_sv_size_by_depth', '--min_sv_size_by_depth', help='min sv size by depth', required=True, type=int)
    parser.add_argument('-gender', '--gender', help='gender', required=True)

    parser.add_argument('-disable_dp_filter', '--disable_dp_filter', help='disable DP filter', required=False, type=int)

    init_logger()

    options = parser.parse_args()
    chain_caller = ChainCaller(options)
    chain_caller.run()
