import os
from utility import *
import csv


class TermSvCallerOptions:
    def __init__(self, name, in_bam, depth_file, fai_file, output_prefix, term_threshold, min_clip_size, fastq, term_seq_file_basepath, ref, out_bed2_file, buf_size):
        self.name = name
        self.in_bam = in_bam
        self.depth_file = depth_file
        self.fai_file = fai_file
        self.output_prefix = output_prefix
        self.term_threshold = term_threshold
        self.min_clip_size = min_clip_size
        self.fastq = fastq
        self.term_seq_file_basepath = term_seq_file_basepath
        self.ref = ref
        self.out_bed2_file = out_bed2_file
        self.buf_size = buf_size


class TermSvCaller:
    def __init__(self, options):
        self.options = options

        output_path = os.path.dirname(options.output_prefix)
        try:
            os.makedirs(output_path)
        except:
            pass

        self.term_region = []
        self.bp_reads = {'start': {}, 'end': {}}
        self.called_reads = {}

    def load_depth_region(self):
        options = self.options

        self.depth_region = []
        with open(options.depth_file) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                arr = line.strip().split()
                arr[0] = 'X' if arr[0] == '23' else 'Y' if arr[0] == '24' else arr[0]

                chrom, start, end, score, region_type = arr[0], int(arr[1]), int(arr[2]), int(float(arr[3])), arr[4]

                if region_type not in ['DEL']:
                    continue

                self.depth_region.append({'chrom': chrom, 'start': start, 'end': end, 'region_type': region_type})

        #print('loaded depth_region', self.depth_region)

    def load_chrom_info(self):
        options = self.options

        self.chrom_info = {}
        with open(options.fai_file) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                arr = line.strip().split()

                chrom, size = arr[0], int(arr[1])

                self.chrom_info[chrom] = size

    def find_term_region(self):
        options = self.options

        self.term_region = []

        for region in self.depth_region:
            chrom = region['chrom']
            start = region['start']
            end = region['end']

            chrom_size = self.chrom_info[chrom]

            if start < options.term_threshold:
                self.term_region.append({'chrom': chrom, 'pos': end, 'side': 'start'})
            if end > chrom_size - options.term_threshold:
                self.term_region.append({'chrom': chrom, 'pos': start, 'side': 'end'})

    def get_region_file_name(self, region, type):
        options = self.options

        if type != 'ref':
            prefix = '%s_%s_%s_%s' % (options.output_prefix, region['chrom'], region['pos'], region['side'])
            return '%s.%s' % (prefix, type)
        else:
            return 'assembled_ref/HG001_assembled.fa'

    def gen_fasta(self, term_region):
        options = self.options

        in_bam = pysam.AlignmentFile(options.in_bam, 'rb')
        fasta_file = self.get_region_file_name(term_region, 'fasta')
        fasta = open(fasta_file, 'w')

        term_chrom = term_region['chrom']
        term_start = 1 if term_region['pos'] <= options.buf_size else term_region['pos'] - options.buf_size
        term_end = term_region['pos'] + options.buf_size
        term_side = term_region['side']

        reads = {}
        for read in in_bam.fetch(term_chrom, term_start, term_end):
            if read.is_secondary:
                continue

            pos = read.reference_end if term_side == 'end' else read.reference_start

            if term_start <= pos <= term_end:
                tuples = read.cigartuples

                clip = get_clip_seq(read, options.fastq[:-9])
                if clip[term_side]['length'] > options.min_clip_size:
                    sa_str = read.get_tag("SA").strip(';') if read.has_tag('SA') else ''
                    strand = '+' if not read.is_reverse else '-'
                    clip_seq = clip[term_side]['seq']
                    #breakpoint = read.reference_end if term_side == 'end' else read.reference_start
                    read_name = read.query_name

                    if read_name not in reads:
                        reads[read_name] = 1
                    else:
                        reads[read_name] += 1

                    new_read_name = '%s_%s' % (read_name, reads[read_name])
                    self.bp_reads[term_side][new_read_name] = {'bp': pos}

                    print('>%s' % (new_read_name), file=fasta)
                    print('%s' % (clip_seq), file=fasta)

    def align(self, region):
        options = self.options

        samtools = get_var('common', 'samtools')
        # term_seq_file = '%s/%s_%s.fa' % (options.term_seq_file_basepath, region['chrom'], region['side'])

        fasta = self.get_region_file_name(region, 'fasta')
        ref = self.get_region_file_name(region, 'ref')
        out_bam = self.get_region_file_name(region, 'bam')

        cmd = "%s -Y -t 48 -a %s %s | %s sort -@ 24 -o %s - && %s index -@ 24 %s &> %s.log" % \
            (get_var('common', 'minimap2'), ref, fasta, samtools, out_bam, samtools, out_bam, out_bam)

        run_shell_cmd(cmd)

    def get_read_support_info(self, read):
        options = self.options

        support_info = {'name': read.query_name, 'is_support': False, 'is_assembled_seq': False, 'seq_chrom': None, 'seq_pos': 0}

        # read_name = read.query_name
        ref_name = read.reference_name

        arr = ref_name.split('_')
        if arr[0] == 'assembled':
            chrom = arr[1]
            pos = int(arr[2])
            support_info['is_support'] = True
            support_info['is_assembled_seq'] = True
            support_info['seq_chrom'] = chrom
            support_info['seq_pos'] = pos
            support_info['AS'] = int(read.get_tag('AS'))
        else:
            chrom = ref_name
            chrom_size = self.chrom_info[chrom]

            pos = read.reference_start
            if pos < options.term_threshold or pos > chrom_size - options.term_threshold:
                support_info['is_support'] = True
                support_info['is_assembled_seq'] = False
                support_info['seq_chrom'] = chrom
                support_info['seq_pos'] = pos
                support_info['AS'] = int(read.get_tag('AS'))

        return support_info

    def check_result(self, region):
        term_chrom, term_pos, term_side = region['chrom'], region['pos'], region['side']

        bam_file = self.get_region_file_name(region, 'bam')

        bam = pysam.AlignmentFile(bam_file, 'rb')

        for read in bam.fetch():
            if read.is_secondary:
                continue

            read_support_info = self.get_read_support_info(read)
            if read_support_info['is_support']:
                read_name = read.query_name

                if read_name in self.bp_reads[term_side]:
                    if read_name not in self.called_reads or read_support_info['AS'] > self.called_reads[read_name]['AS']:
                        read_length = read.infer_read_length()

                        if term_side == 'end':
                            self.called_reads[read_name] = {
                                'chrom': term_chrom,
                                'bp1': self.bp_reads[term_side][read_name]['bp'],
                                'chrom2': term_chrom,
                                'bp2': 0,
                                'AS': read_support_info['AS'],
                                'read_name': get_short_name(read_name),
                                'read_length': read_length,
                            }
                        else:
                            self.called_reads[read_name] = {
                                'chrom': term_chrom,
                                'bp1': 0,
                                'chrom2': term_chrom,
                                'bp2': self.bp_reads[term_side][read_name]['bp'],
                                'AS': read_support_info['AS'],
                                'read_name': get_short_name(read_name),
                                'read_length': read_length,
                            }

    def write_results(self):
        options = self.options
        called_reads = self.called_reads

        with open(options.out_bed2_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            for read_name in called_reads:
                called_read = called_reads[read_name]
                row = [
                    called_read['chrom'],
                    called_read['bp1'],
                    called_read['chrom2'],
                    called_read['bp2'],
                    'TERM-DEL',
                    called_read['read_name'],
                    called_read['read_length'],
                    called_read['AS'],
                ]
                writer.writerow(row)

    def run(self):
        self.load_depth_region()
        self.load_chrom_info()
        self.find_term_region()

        for region in self.term_region:
            self.gen_fasta(region)
            self.align(region)
            self.check_result(region)

        self.write_results()

        return self.called_reads


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='run')
    parser.add_argument('-name', '--name', help='name', required=True)
    parser.add_argument('-bam', '--bam', help='bam', required=True)
    parser.add_argument('-depth_file', '--depth_file', help='depth file', required=True)
    parser.add_argument('-fai_file', '--fai_file', help='ref. index file', required=True)
    parser.add_argument('-output_prefix', '--output_prefix', help='output prefix', required=True)
    parser.add_argument('-term_threshold', '--term_threshold', help='term threshold', required=True, type=int)
    parser.add_argument('-min_clip_size', '--min_clip_size', help='min clip size', required=True, type=int)
    parser.add_argument('-fastq', '--fastq', help='fastq', required=True)
    parser.add_argument('-term_seq_file_basepath', '--term_seq_file_basepath', help='telomere seq file basepath', required=True)
    parser.add_argument('-ref', '--ref', help='ref', required=True)
    parser.add_argument('-buf_size', '--buf_size', help='buf size', required=True, type=int)

    options = parser.parse_args()
    term_sv_caller_options = TermSvCallerOptions(
        options.name,
        options.bam,
        options.depth_file,
        options.fai_file,
        options.output_prefix,
        options.term_threshold,
        options.min_clip_size,
        options.fastq,
        options.term_seq_file_basepath,
        options.ref,
        'term_sv_caller.bed2',
        options.buf_size,
    )
    term_sv_caller = TermSvCaller(term_sv_caller_options)
    term_sv_caller.run()
