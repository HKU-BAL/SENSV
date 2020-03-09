from __future__ import print_function

from src.utility import *

class FindSingleEndRead:
    cType = get_var('denovo_common','cType')

    def __init__(self, options):
        self.bam = options.bam
        self.fastq = options.fastq
        self.pos_chrom = options.pos_chrom
        self.pos_start = int(options.pos_start)
        self.pos_end = int(options.pos_end)
        self.min_clip_size = int(options.min_clip_size)

        self.out_fastq = options.out_fastq
        self.out_fasta = options.out_fasta
        self.out_txt = options.out_txt

        if options.clip_side:
            self.clip_side = options.clip_side
        else:
            self.clip_side = 'end'

    def gen(self):
        header = ['name','strand','clip_size','breakpoint','sa_str','clip_seq']
        #header = ['name','strand','clip_size','breakpoint','clip_seq']

        fastq_read = {}
        bam = pysam.AlignmentFile(self.bam, 'rb')

        if self.out_fasta:
           fout_fasta = open(self.out_fasta, 'w')
        if self.out_fastq:
           fout_fastq = open(self.out_fastq, 'w')
        if self.out_txt:
           fout_txt = open(self.out_txt, 'w')

        for read in bam.fetch(self.pos_chrom, self.pos_start, self.pos_end):
            if read.is_secondary:
                continue

            pos = read.reference_end if self.clip_side == 'end' else read.reference_start

            if self.pos_start <= pos <= self.pos_end:
                tuples = read.cigartuples

                clip = get_clip_seq(read, self.fastq[:-9])
                if clip[self.clip_side]['length'] > self.min_clip_size:
                    sa_str = read.get_tag("SA").strip(';') if read.has_tag('SA') else ''
                    strand = '+' if not read.is_reverse else '-'
                    #clip_seq = clip[self.clip_side]['seq'] if self.clip_side == 'end' else clip[self.clip_side]['seq'][::-1]
                    clip_seq = clip[self.clip_side]['seq']
                    breakpoint = read.reference_end if self.clip_side == 'end' else read.reference_start

                    #print(read.reference_start, read.reference_end, breakpoint, clip_seq[:200])

                    """
                    if sa_str:
                        print('name', read.query_name, 'sa_str', sa_str)
                        sa_arr = sa_str.split(';')

                        is_ok = 1
                        for sa in sa_arr:
                            arr = sa.split(',')
                            sa_chrom = arr[0]
                            sa_start = int(arr[1])

                            if sa_chrom == self.pos_chrom:
                                if abs(sa_start - breakpoint) < 1000000:
                                   is_ok = 0
                                   break

                        if not is_ok:
                            continue
                    """

                    if self.out_fasta:
                        print('>%s' % read.query_name, file=fout_fasta)
                        print('%s' % clip_seq, file=fout_fasta)

                    if self.out_txt:
                        output = {}
                        output['name'] = read.query_name
                        output['strand'] = strand
                        output['clip_size'] = clip[self.clip_side]['length']
                        output['breakpoint'] = breakpoint
                        output['sa_str'] = sa_str
                        output['clip_seq'] = clip_seq

                        print('\t'.join([str(output[i]) for i in header]), file=fout_txt)

                    if self.out_fastq:
                        if read.query_name in fastq_read:
                            continue
                        fastq_read[read.query_name] = 1

                        rows = get_seq_from_fastq(read.query_name, self.fastq, '', True)
                        print('\n'.join(rows), file=fout_fastq)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='run')
    parser.add_argument('-bam', '--bam', help='bam', required=True)
    parser.add_argument('-fastq', '--fastq', help='fastq', required=True)
    parser.add_argument('-pos_chrom', '--pos_chrom', help='pos chrom', required=True)
    parser.add_argument('-pos_start', '--pos_start', help='pos start', required=True)
    parser.add_argument('-pos_end', '--pos_end', help='pos end', required=True)
    parser.add_argument('-min_clip_size', '--min_clip_size', help='min clip size', required=True)

    parser.add_argument('-clip_side', '--clip_side', help='clip side', required=False)

    parser.add_argument('-out_fasta', '--out_fasta', help='fasta', required=False)
    parser.add_argument('-out_fastq', '--out_fastq', help='fastq', required=False)
    parser.add_argument('-out_txt', '--out_txt', help='txt', required=False)

    options = parser.parse_args()
    find_single_end_read = FindSingleEndRead(options)
    find_single_end_read.gen()
