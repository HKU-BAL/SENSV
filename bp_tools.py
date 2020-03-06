from argparse import ArgumentParser
from multiprocessing import Pool

from dp import DpOptions, DP
from utility import (
    get_var,
    get_seq_from_fastq,
    get_sorted_sv_str_list,
    rev_comp,
    init_logger
)

merge_dist = 100


class BpToolsOptions:
    def __init__(self, input_bed2, output_bed2, fastq_prefix, min_sv_size, max_sv_size, nprocs):
        self.input_bed2 = input_bed2
        self.output_bed2 = output_bed2
        self.fastq_prefix = fastq_prefix
        self.min_sv_size = min_sv_size
        self.max_sv_size = max_sv_size
        self.nprocs = nprocs


class BpTools:
    def __init__(self, options):
        self.options = options

        self.gap_buf_size = int(get_var('bp_tools', 'gap_buf_size'))
        self.buf_size = int(get_var('bp_tools', 'buf_size'))
        self.min_score_read_length_ratio = float(get_var('bp_tools', 'min_score_read_length_ratio'))

    def __call__(self, param_str):
        options = self.options

        arr = param_str.split('|')
        sv_str, query_name, query_strand, bp_query_pos = arr[0], arr[1], arr[2], int(arr[3])

        query_seq = get_seq_from_fastq(query_name, options.fastq_prefix)

        if query_strand == '-':
            query_seq = rev_comp(query_seq)
            bp_query_pos = len(query_seq) - bp_query_pos

        query_start = 0 if bp_query_pos < self.buf_size else bp_query_pos - self.buf_size
        query_end = bp_query_pos + self.buf_size
        query_seq = query_seq[query_start:query_end]

        dp_options = DpOptions(sv_str, query_name, query_strand, options.fastq_prefix, self.gap_buf_size, query_seq)
        dp = DP(dp_options)
        result = dp.run()
        result['bp_query_pos'] = bp_query_pos
        result['local_bp_query_pos'] = query_end-query_start-self.buf_size
        result['bp_buf_size'] = self.buf_size

        return result

    def merge_bp(self):
        options = self.options
        sv_str_dict = {}

        with open(options.input_bed2) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                arr = line.strip().split()
                chrom, start, chrom2, end, type = arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4]
                score_ratio = 0
                if len(arr) == 10:
                    score_ratio = arr[8]
                sv_str_dict[f'{chrom}_{start}_{chrom2}_{end}_{type}'] = {
                    'value': arr[5:],
                    'score_ratio': score_ratio
                }

        sorted_sv_str_list = get_sorted_sv_str_list(sv_str_dict)

        merged_sv_dict = {}
        merged_sv_str = None
        for sv_str in sorted_sv_str_list:
            arr = sv_str.split('_')
            value = sv_str_dict[sv_str]['value']
            score_ratio = sv_str_dict[sv_str]['score_ratio']
            curr_sv = {
                'chrom': arr[0],
                'start': int(arr[1]),
                'chrom2': arr[2],
                'end': int(arr[3]),
                'type': arr[4],
                'value': value,
                'score_ratio': score_ratio
            }

            if curr_sv['type'] in ['DEL', 'DUP']:
                if options.min_sv_size and curr_sv['end'] - curr_sv['start'] < options.min_sv_size:
                    continue
                if options.max_sv_size and curr_sv['end'] - curr_sv['start'] > options.max_sv_size:
                    continue

            if merged_sv_str:
                arr = merged_sv_str.split('_')
                merged_sv = {
                    'chrom': arr[0],
                    'start': int(arr[1]),
                    'chrom2': arr[2],
                    'end': int(arr[3]),
                    'type': arr[4]
                }
                if (
                    curr_sv['type'] == merged_sv['type'] and
                    curr_sv['chrom'] == merged_sv['chrom'] and
                    curr_sv['chrom2'] == merged_sv['chrom2'] and
                    abs(curr_sv['start'] - merged_sv['start']) < merge_dist and
                    abs(curr_sv['end'] - merged_sv['end']) < merge_dist
                ):
                    if curr_sv['score_ratio'] > merged_sv_dict[merged_sv_str]['score_ratio']:
                        del merged_sv_dict[merged_sv_str]
                        merged_sv_str = sv_str
                        merged_sv_dict[merged_sv_str] = {
                            'value': curr_sv['value'],
                            'score_ratio': curr_sv['score_ratio'],
                        }
                else:
                    merged_sv_str = sv_str
                    merged_sv_dict[merged_sv_str] = {
                        'value': curr_sv['value'],
                        'score_ratio': curr_sv['score_ratio'],
                    }
            else:
                merged_sv_str = sv_str
                merged_sv_dict[merged_sv_str] = {
                    'value': value,
                    'score_ratio': score_ratio,
                }

        sorted_list = get_sorted_sv_str_list(merged_sv_dict)
        with open(options.output_bed2, 'w') as f:
            for key in sorted_list:
                value = merged_sv_dict[key]['value']
                output = key.split('_')
                output.extend(value)
                print('\t'.join(output), file=f)

    def load_param_str_list(self):
        options = self.options

        self.param_str_list = []
        with open(options.input_bed2) as f:
            for line in f:
                if line[0] == '#':
                    continue

                arr = line.strip().split()
                if len(arr) == 9:
                    chrom, start, chrom2, end, sv_type, _supp_type, query_name, query_strand, query_pos = \
                        arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4], arr[5], arr[6], arr[7], arr[8]

                    sv_str = f'{chrom}_{start}_{chrom2}_{end}_{sv_type}'
                    self.param_str_list.append(f'{sv_str}|{query_name}|{query_strand}|{query_pos}')
                elif len(arr) == 7:
                    chrom, start, chrom2, end, sv_type, query_name, query_strand = \
                        arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4], arr[5], arr[6]

                    sv_str = f'{chrom}_{start}_{chrom2}_{end}_{sv_type}'
                    self.param_str_list.append(f'{sv_str}|{query_name}|{query_strand}')

    def refine_bp(self):
        options = self.options

        self.load_param_str_list()

        nprocs = min(len(self.param_str_list), self.options.nprocs)
        pool = Pool(processes=nprocs)
        results = pool.map(self, self.param_str_list)
        pool.close()
        pool.join()

        header = [
            'chrom',
            'genomic_gap_start',
            'chrom2',
            'genomic_gap_end',
            'sv_type',
            'qname',
            'query_len',
            'score',
            'score_read_length_ratio',
            'query_pos',
        ]
        with open(options.output_bed2, 'w') as f:
            for result in results:
                # adjustment of query_pos
                result['query_pos'] = result['bp_query_pos'] + (result['query_pos'] - result['local_bp_query_pos'])

                #print('result', result)
                """
                if result['sv_type'] in ['DEL', 'DUP']:
                    if abs(int(result['genomic_gap_start'])-int(result['genomic_gap_end'])) < 10000:
                        continue
                """

                if result['sv_type'] in ['DUP']:
                    if int(result['genomic_gap_start']) > int(result['genomic_gap_end']):
                        result['genomic_gap_start'], result['genomic_gap_end'] = result['genomic_gap_end'], result['genomic_gap_start']

                # if result['sv_type'] in ['DEL']:
                #    if abs(int(result['genomic_gap_start']) - int(result['genomic_gap_end'])) < 10000:
                #        continue

                """
                if result['sv_type'] not in ['DEL', 'DUP']:
                    continue
                """

                result['score_read_length_ratio'] = '%.2f' % (1. * result['score'] / result['query_len'])

                # if result['genomic_gap_start'] > 0 and result['genomic_gap_end'] > 0:
                row = [str(result[x]) for x in header]
                is_filtered = False
                # if float(result['score_read_length_ratio']) < self.min_score_read_length_ratio:
                if False:
                    is_filtered = True
                elif result['sv_type'] == 'DEL':
                    if (
                        options.min_sv_size and
                        result['genomic_gap_end'] - result['genomic_gap_start'] < options.min_sv_size
                    ) or (
                        options.max_sv_size and
                        result['genomic_gap_end'] - result['genomic_gap_start'] > options.max_sv_size
                    ):
                        is_filtered = True

                if not is_filtered:
                    print('\t'.join(row), file=f)
                else:
                    print('#' + '\t'.join(row), file=f)


if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-input_bed2', '--input_bed2', help='input bed2', required=True)
    parser.add_argument('-output_bed2', '--output_bed2', help='output bed2', required=True)
    parser.add_argument('-fastq_prefix', '--fastq_prefix', help='fastq prefix', required=True)
    parser.add_argument('-min_sv_size', '--min_sv_size', help='min sv size', required=True, type=int)
    parser.add_argument('-max_sv_size', '--max_sv_size', help='max sv size', required=True, type=int)

    parser.add_argument('-action', '--action', help='action', required=False)

    init_logger()

    options = parser.parse_args()
    bp_tools = BpTools(options)
    if options.action and options.action == 'merge_bp':
        bp_tools.merge_bp()
    else:
        bp_tools.refine_bp()
