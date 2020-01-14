from __future__ import print_function
from argparse import ArgumentParser

def run(options):
    with open(options.in_bed2, 'r') as f_in:
        with open(options.out_bed, 'w') as f_out:
            for line in f_in:
                if not line.strip() or line[0] == '#':
                    continue

                arr = line.strip().split()
                chrom, start, chrom2, end, sv_type = arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4]

                """
                if sv_type not in ["DEL-FROM-INS-FROM", "DEL-TO-INS-FROM", "DEL-FROM-INS-TO", "DEL-TO-INS-TO"]:
                    row = [chrom, str(start - options.search_size), str(end + options.search_size)]
                    print('\t'.join(row), file=f_out)
                else:
                    row = [chrom, str(start - options.search_size), str(start + options.search_size)]
                    print('\t'.join(row), file=f_out)

                    row = [chrom2, str(end - options.search_size), str(end + options.search_size)]
                    print('\t'.join(row), file=f_out)
                """

                """
                if sv_type in ['DEL']:
                    row = [chrom, str(start - options.search_size), str(end + options.search_size)]
                    print('\t'.join(row), file=f_out)
                else:
                    row = [chrom, str(start - options.search_size), str(start + options.search_size)]
                    print('\t'.join(row), file=f_out)

                    row = [chrom2, str(end - options.search_size), str(end + options.search_size)]
                    print('\t'.join(row), file=f_out)

                """

                if sv_type in ['DEL'] and end - start < options.search_size * 2:
                    temp = 1 if start <= options.search_size else start - options.search_size
                    row = [chrom, str(temp), str(end + options.search_size)]
                    print('\t'.join(row), file=f_out)
                else:
                    temp = 1 if start <= options.search_size else start - options.search_size
                    row = [chrom, str(temp), str(start + options.search_size)]
                    print('\t'.join(row), file=f_out)

                    temp = 1 if end <= options.search_size else end - options.search_size
                    row = [chrom2, str(temp), str(end + options.search_size)]
                    print('\t'.join(row), file=f_out)

if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-in_bed2', '--in_bed2', help='in bed2', required=True)
    parser.add_argument('-out_bed', '--out_bed', help='out bed', required=True)
    parser.add_argument('-search_size', '--search_size', help='search size', required=True, type=int)

    options = parser.parse_args()
    run(options)
