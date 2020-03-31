import sys
from utility import run_shell_cmd
from argparse import ArgumentParser


def x(row, z):
    return float(row[ord(z) - ord('a')])

def get_chromsize(chrom, chrom_dict):
    if chrom in chrom_dict:
        return (chrom_dict[chrom], chrom_dict)
    else:
        cmd = (f'faidx -i chromsizes {args.ref_file_path} {chrom}')
        line = run_shell_cmd(cmd)
        line = line.strip().split()
        chrom_dict[line[0]] = line[1]

        return (chrom_dict[chrom], chrom_dict)


if __name__ == "__main__":
    parser = ArgumentParser(description='A script to see see what split reads look like >.0')
    parser.add_argument('--header_file_path', help='header file path', required=True, type=str, default='')
    parser.add_argument('--ref_file_path', help='fasta file path', required=True, type=str, default='')
    args = parser.parse_args()

    chromsize_dict = {}

    INF = 999999

    header_file_path = args.header_file_path
    print(f'header_file_path: #{header_file_path}#', file=sys.stderr)

    f = open(header_file_path).read().splitlines()
    for ff in f:
        print(ff)

    for row in sys.stdin:
        row = row.split()
        if len(row[0].split('_')) < 5:
            continue
        #if row[0][-4:] != "_DEL" and row[0][-4:] != "_DUP" and row[0][-3:] != "-TO" and row[0][-5] != "-FROM":
            # print(f'row continue: #{row}#', file=sys.stderr)
        #    continue

        # TYPE = row[0][-3:]
        # if row.find("DEL") == -1: continue
        chromosome, start_position, end_chromosome, end_position, TYPE = row[0].split('_')
        if TYPE == "TERM-DEL":
            if start_position == "0":
                SVLEN = int(end_position)-1
                print(
                    f"{chromosome}\t1\t.\tN\t<DEL>\t0\tPASS\tSVLEN={SVLEN};SVTYPE=DEL;END={end_position}\tGT\t./."
                )
            else:
                chromosome_length, chromsize_dict = get_chromsize(chromosome, chromsize_dict)
                SVLEN = int(chromosome_length)-int(start_position)
                print(
                    f"{chromosome}\t{start_position}\t.\tN\t<DEL>\t0\tPASS\tSVLEN={SVLEN};SVTYPE=DEL;END={chromosome_length}\tGT\t./."
                )
            continue
        start_position = int(start_position)
        end_position = int(end_position)

        QUAL = 0
        """sv!sv_length!read_name!length!origAS!AS!AS_inc(%)!AS_left!AS_right!AS-to-readlen!bp_start/bp_end!left_clip_length!right_clip_length!left_query_span!right_query_span!left_ref_span!right_ref_span!left_match_span!right_match_span!orig_MAPQ!MAPQ!MAPQ_left!MAPQ_right!INDEL_left_length!INDEL_left_count!INDEL_left_span!INDEL_right_length!INDEL_right_count!INDEL_right_span"""

        # print row
        mapq = x(row, 't')
        if mapq < 60:
            QUAL = -INF

        t = abs(x(row, 'p') - x(row, 'n')) + abs(x(row, 'q') - x(row, 'o'))
        if t >= 2500:
            QUAL = -INF
        else:
            QUAL += (2500-t) / 100

        t1 = 0 if x(row, 'p') == 0 else x(row, 'r') / x(row, 'p')
        t2 = 0 if x(row, 'q') == 0 else x(row, 's') / x(row, 'q')
        #t = 2000.0/((1e10 if t1 ==0 else 1.0/t1)+(1e10 if t2==0 else 1.0/t2))
        QUAL += (t1 + t2 - 1.50) * 100
        if QUAL < 0:
            QUAL = 0

        #QUAL += t
        #QUAL = QUAL / 10000
        QUAL = int(QUAL)
        SVLEN = end_position - start_position

        if TYPE == "DEL" or TYPE == "DUP":
            print(
                f"{chromosome}\t{start_position}\t.\tN\t<{TYPE}>\t{QUAL}\tPASS\tSVLEN=-{SVLEN};SVTYPE={TYPE};END={end_position}\tGT\t./."
            )
        else:
            TYPE = TYPE.split('-')
            if TYPE[1] == "TO" and TYPE[3] == "TO":
                ALT = "N]" + end_chromosome + ":" + str(end_position) + "]"
            elif TYPE[1] == "TO" and TYPE[3] == "FROM":
                ALT = "N[" + end_chromosome + ":" + str(end_position) + "["
            elif TYPE[1] == "FROM" and TYPE[3] == "TO":
                ALT = "]" + end_chromosome + ":" + str(end_position) + "]N"
            elif TYPE[1] == "FROM" and TYPE[3] == "FROM":
                ALT = "[" + end_chromosome + ":" + str(end_position) + "[N"

            print(
                f"{chromosome}\t{start_position}\t.\tN\t{ALT}\t{QUAL}\tPASS\tSVLEN=0;SVTYPE=BND\tGT\t./."
            )
