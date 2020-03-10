import os
import sys
from pathlib import Path

from utility import data_directory


def x(row, z):
    return float(row[ord(z) - ord('a')])


if __name__ == "__main__":
    INF = 999999

    header_file_path = data_directory() / 'output' / 'header'
    print(f'header_file_path: #{header_file_path}#', file=sys.stderr)

    f = open(header_file_path).read().splitlines()
    for ff in f:
        print(ff)

    for row in sys.stdin:
        row = row.split()
        if row[0][-4:] != "_DEL" and row[0][-4:] != "_DUP":
            # print(f'row continue: #{row}#', file=sys.stderr)
            continue

        TYPE = row[0][-3:]
        # if row.find("DEL") == -1: continue
        chromosome, start_position, _, end_position, _ = row[0].split('_')
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

        print(
            f"{chromosome}\t{start_position}\t.\tN\t<{TYPE}>\t{QUAL}\tPASS\tSVLEN=-{SVLEN};SVTYPE={TYPE};END={end_position}\tGT\t./."
        )
