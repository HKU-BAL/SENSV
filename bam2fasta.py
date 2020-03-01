import pysam
from argparse import ArgumentParser
from multiprocessing import Pool

from utility import *


class Bam2Fasta:
    def __init__(self, options):
        self.input_bam = options.input_bam
        self.input_fastq = options.input_fastq
        self.output_fasta = options.output_fasta

    def __call__(self, read_name):
        read_seq = get_seq_from_fastq(read_name, self.input_fastq, '')

        return {'read_name': read_name, 'read_seq': read_seq}

    def read_bam(self):
        read_info = {}

        bam = pysam.AlignmentFile(self.input_bam, 'rb')
        for read in bam.fetch():
            # ** need to include secondary reads
            # if read.is_secondary:
            #    continue

            if read.query_name not in read_info:
                if 'H' not in read.cigarstring:
                    read_info[read.query_name] = read.query_sequence
                else:
                    read_info[read.query_name] = ''
            elif not read_info[read.query_name] and 'H' not in read.cigarstring:
                read_info[read.query_name] = read.query_sequence

        self.read_info = read_info

    def read_fastq(self):
        missing_seq = {}
        for read_name in self.read_info:
            if not self.read_info[read_name]:
                missing_seq[read_name] = ''

        pool_size = 100
        pool = Pool(pool_size)
        results = pool.map(self, missing_seq.keys())
        pool.close()
        pool.join()

        for result in results:
            self.read_info[result['read_name']] = result['read_seq']

    def gen_fasta(self):
        with open(self.output_fasta, 'w') as fasta:
            for read_name in self.read_info:
                read_seq = self.read_info[read_name]

                print('@%s' % (read_name), file=fasta)
                print('%s' % (read_seq), file=fasta)

    def run(self):
        self.read_bam()
        self.read_fastq()
        self.gen_fasta()


if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-input_bam', '--input_bam', help='in bam', required=True)
    parser.add_argument('-input_fastq', '--input_fastq', help='in fastq', required=True)
    parser.add_argument('-output_fasta', '--output_fasta', help='out fasta', required=True)

    options = parser.parse_args()
    bam2fasta = Bam2Fasta(options)
    bam2fasta.run()
