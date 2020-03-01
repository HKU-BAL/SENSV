import os
import shutil
import csv
from argparse import ArgumentParser
from os.path import dirname

from utility import *

from split_read_caller import *
from chain_caller import *
from altref import *
from depth_info import *
from term_sv_caller import *
from bp_tools import *

"""
class SENSV is the integration class for all the steps of SV calling
"""


class SENSV:
    FASTQ_INDEX_SCRIPT = './index_fastq.sh'
    MERGE_SV_SCRIPT = './merge_sv.sh'

    def __init__(self, options):
        self.load_config()

        config = self.config

        self.sample_name = options.sample_name
        self.fastq_file_orig = options.fastq_file
        self.output_prefix = options.output_prefix
        #self.sensitive_mode = options.sensitive_mode

        self.min_sv_size = options.min_sv_size if options.min_sv_size else int(config['default_min_sv_size'])
        self.min_sv_size_by_depth = int(config['default_min_sv_size_by_depth'])
        self.max_sv_size = options.max_sv_size if options.max_sv_size else int(config['default_max_sv_size'])

        self.fastq_file = f'{options.output_prefix}.fastq.gz'
        self.fastq_prefix = options.output_prefix
        self.depth_file = '%s.depth' % (options.output_prefix)
        self.depth_file_filtered = '%s.depth_filtered' % (options.output_prefix)
        self.minimap2_bam_file = '%s_minimap2.bam' % (options.output_prefix)
        self.minimap2_filter_bam_file = '%s_minimap2_filter.bam' % (options.output_prefix)
        self.minimap2_raw_bam_file = '%s_minimap2_raw.bam' % (options.output_prefix)
        #self.ngmlr_bam_file = '%s_ngmlr.bam' % (options.output_prefix)
        self.filter_read_file = '%s_filter_read.fastq' % (options.output_prefix)
        self.filter_read_id_file = '%s_filter_read_id_file' % (options.output_prefix)
        self.filter_read_bam_file = '%s_filter_read.bam' % (options.output_prefix)
        self.chain_file = '%s.chain' % (options.output_prefix)

        # working files
        self.bed2_file = '%s.bed2' % (options.output_prefix)
        self.split_bed2_file = '%s.split.bed2' % (options.output_prefix)
        self.dp_bed2_file = '%s.dp.bed2' % (options.output_prefix)
        self.validate_bed2_file = '%s.validate.bed2' % (options.output_prefix)
        self.merged_bed2_file = '%s.merge.bed2' % (self.output_prefix)

        self.term_bed2_file = '%s.term.bed2' % (self.output_prefix)

        # options for calling modified_minimap2
        self.minimap2_opt1 = '-Y -t 48 -z 200 --MD -a %s %s.fastq.gz' % (config['ref'], self.output_prefix)
        self.minimap2_opt2 = '-5 %s' % (self.chain_file)

        # make working dir
        outputPath = os.path.dirname(self.output_prefix)
        try:
            os.makedirs(outputPath)
        except:
            pass

        self.disable_dp_filter = options.disable_dp_filter
        self.disable_gen_altref_bam = options.disable_gen_altref_bam
        self.target_sv_type = options.target_sv_type.split(',') if options.target_sv_type else 'DUP,DEL'
        self._gender = None

        # less than 10k bed2
        self.lt_10k_bed2 = '%s_lt_10k.bed2' % (self.output_prefix)
        self.lt_10k_merge_bed2 = '%s_lt_10k.merge.bed2' % (self.output_prefix)

        # for chrom info
        self.fai_file = '%s.fai' % (config['ref'])

        self.final_result = '%s_final.result' % (self.output_prefix)
        self.final_bed2 = '%s_final.bed2' % (self.output_prefix)

        self.disable_depth_analysis = False
        if self.max_sv_size and self.max_sv_size <= self.min_sv_size_by_depth:
            self.disable_depth_analysis = True

        if self.disable_depth_analysis:
            print('!!! disabled depth analysis')

            # init files
            files = [self.depth_file_filtered, '%s.dp.bed2' % (self.output_prefix)]
            for file in files:
                cmd = 'rm -f %s; touch %s' % (file, file)
                run_shell_cmd(cmd)

    def load_config(self):
        self.config = {
            'ref_ver': get_var('common', 'ref_ver'),
            'samtools': get_var('common', 'samtools'),
            'ref': get_var('common', 'ref_%s' % get_var('common', 'ref_ver')),
            'default_min_sv_size': int(get_var('default_value', 'min_sv_size')),
            'default_min_sv_size_by_depth': int(get_var('default_value', 'min_sv_size_by_depth')),
            'default_max_sv_size': int(get_var('default_value', 'max_sv_size')),
            'sv_type': get_var('default_value', 'sv_type'),
            'my_minimap2': get_var('dp', 'my_minimap2'),
            'chrom_list': get_var('denovo_common', 'chrom_list'),
            'depth_ref': get_var('depth', 'depth_ref'),
        }

    @log(message="step 0 (convert fastq.gz to bgzip and index it)")
    def index_fastq(self):
        run_shell_cmd(f'{self.FASTQ_INDEX_SCRIPT} {self.fastq_file_orig} {self.output_prefix}')

    @log(message="step 1 (generate bam and chain files)")
    def gen_bam_and_chain_file(self):
        """
        generate bam file and chain file
        - input: fastq
        - output: bam, chain file
        """
        my_minimap2 = self.config['my_minimap2']
        samtools = self.config['samtools']

        cmd = "%s %s %s | %s sort -@ 48 -o %s - && %s index -@ 48 %s" % \
              (my_minimap2, self.minimap2_opt1, self.minimap2_opt2, samtools,
               self.minimap2_bam_file, samtools, self.minimap2_bam_file)
        # logging.info(f'cmd: #{cmd}#')

        run_shell_cmd(cmd)

    @property
    def gender(self):
        if self._gender is not None:
            return self._gender

        def reads_count_for_chromosome(chromosome):
            samtools = self.config['samtools']
            bam = self.minimap2_bam_file
            return int(run_shell_cmd(f'{samtools} view -F 260 {bam} {chromosome} | wc -l'))

        X_count = reads_count_for_chromosome("X")
        Y_count = reads_count_for_chromosome("Y")

        self._gender = 'f' if X_count / (Y_count + 0.001) > 100 else 'm'
        return self._gender

    @log(message="step 2 (generate depth file)")
    def gen_depth_file(self):
        """
        generate depth file
        - input: bam
        - output: depth file
        """
        basedir = os.path.dirname(os.path.realpath('__file__'))
        depth_path = basedir + "/depth"

        current_dir = os.getcwd()

        ref_file = self.config['depth_ref']

        if self.output_prefix[0] == "/":
            cmd = 'cd %s && python depth_normalize_gender_thread.py %s %s %s %s' % (
                depth_path, self.minimap2_bam_file, ref_file, self.gender, self.output_prefix)
        else:
            cmd = 'cd %s && python depth_normalize_gender_thread.py %s/%s %s %s %s/%s' % (
                depth_path, current_dir, self.minimap2_bam_file, ref_file, self.gender, current_dir, self.output_prefix)
        # logging.info(f"cmd: #{cmd}#")

        run_shell_cmd(cmd)

    @log(message="step 3 (find bp from cigarstring and split reads)")
    def find_bp_by_cigar_str_and_split_read(self):
        term_threshold = 300000
        chromosome_list = self.config['chrom_list']

        SplitReadCaller(
            SplitReadCallerOptions(
                self.minimap2_bam_file,
                chromosome_list,
                self.min_sv_size,
                self.max_sv_size,
                self.target_sv_type,
                self.bed2_file,
                self.depth_file_filtered,
                term_threshold,
                self.lt_10k_bed2,
                self.fai_file,
            )
        ).run()

    @log(message="step 4 (filter by depth and sv size)")
    def filter_by_depth_and_sv_size(self):
        filtered_sv = []
        min_sv_size_for_depth = int(get_var('depth', 'min_sv_size_for_depth'))
        max_deviate_size = int(get_var('depth', 'max_deviate_size'))

        depth_region = []
        with open(self.depth_file_filtered) as f:
            for line in f:
                if not line.strip():
                    continue

                arr = line.strip().split()
                chrom, start, end = arr[0], int(float(arr[1])), int(float(arr[2]))

                depth_region.append({
                    'chrom': chrom,
                    'start1': start - max_deviate_size,
                    'start2': start + max_deviate_size,
                    'end1': end-max_deviate_size,
                    'end2': end+max_deviate_size
                })

        with open(self.bed2_file) as f:
            for line in f:
                if not line.strip():
                    continue

                arr = line.strip().split()
                chrom, start, chrom2, end, sv_type, supp_type, read_name, read_strand, query_pos = \
                    arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4], arr[5], arr[6], arr[7], int(arr[8])

                svSize = end-start
                if supp_type != "adhoc" and svSize >= min_sv_size_for_depth:
                    match = False
                    for depth in depth_region:
                        if (
                            chrom == depth['chrom'] and
                            depth['start1'] <= start <= depth['start2'] and
                            depth['end1'] <= end <= depth['end2']
                        ):
                            match = True
                            break

                    if match:
                        filtered_sv.append({
                            'chrom': chrom,
                            'start': start,
                            'chrom2': chrom2,
                            'end': end,
                            'sv_type': sv_type,
                            'supp_type': supp_type,
                            'read_name': read_name,
                            'read_strand': read_strand,
                            'query_pos': query_pos,
                        })
                else:
                    filtered_sv.append({
                        'chrom': chrom,
                        'start': start,
                        'chrom2': chrom2,
                        'end': end,
                        'sv_type': sv_type,
                        'supp_type': supp_type,
                        'read_name': read_name,
                        'read_strand': read_strand,
                        'query_pos': query_pos,
                    })

        os.rename(self.bed2_file, '%s_raw' % (self.bed2_file))

        with open(self.bed2_file, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['chrom', 'start', 'chrom2', 'end', 'sv_type',
                                                   'supp_type', 'read_name', 'read_strand', 'query_pos'], delimiter='\t')
            for sv in filtered_sv:
                writer.writerow(sv)

        print('len(filtered_sv)', len(filtered_sv))

    @log(message="step 5 (refine bp by dp)")
    def refine_bp_by_dp(self):
        tmp_bed2 = '%s_tmp' % (self.split_bed2_file)

        options = BpToolsOptions(self.bed2_file, tmp_bed2, self.fastq_prefix, self.min_sv_size, self.max_sv_size)
        bp_tools = BpTools(options)
        bp_tools.refine_bp()

        options = BpToolsOptions(tmp_bed2, self.split_bed2_file, self.fastq_prefix, self.min_sv_size, self.max_sv_size)
        bp_tools = BpTools(options)
        bp_tools.merge_bp()

    @log(message="step 6 (validate bp by local realign)")
    def validate_bp_by_local_realign(self):
        action = 'local_validate'
        options = AltrefOptions(self.sample_name, self.split_bed2_file, self.fastq_file,
                                self.minimap2_bam_file, self.split_bed2_file, action)
        altref = Altref(options)
        altref.run()

        cmd = 'mv %s_filtered.bed2 %s' % (self.split_bed2_file, self.validate_bed2_file)
        run_shell_cmd(cmd)

    @log(message="step 7 (find bp from chain and dp)")
    def find_bp_by_chain_and_dp(self):
        ChainCaller(
            ChainCallerOptions(
                self.depth_file_filtered,
                self.chain_file,
                self.fastq_file,
                self.minimap2_bam_file,
                self.output_prefix,
                self.min_sv_size_by_depth,
                self.disable_dp_filter,
                self.gender,
            )
        ).run()

    @log(message="step 8 (cluster and merge cigarRead in bed2_file)")
    def merge_bed2(self):
        merged_bed2_file_tmp = '%s_tmp' % (self.merged_bed2_file)

        shutil.copy2(self.validate_bed2_file, merged_bed2_file_tmp)
        with open(merged_bed2_file_tmp, 'a') as outfile:
            # read cigarRead and adhoc in bed2_file
            with open(self.bed2_file) as infile:
                for line in infile:
                    arr = line.strip().split()
                    supp_type = arr[5]

                    if supp_type not in ['cigar', 'adhoc']:
                        continue

                    outfile.write(line)

            with open(self.dp_bed2_file) as infile:
                outfile.write(infile.read())

            if os.path.getsize(self.lt_10k_bed2):
                cmd = '%s %s %s' % (self.MERGE_SV_SCRIPT, self.lt_10k_bed2, self.lt_10k_merge_bed2)

                run_shell_cmd(cmd)

                with open(self.lt_10k_merge_bed2) as infile:
                    outfile.write(infile.read())

            # if os.path.isfile(self.term_bed2_file) and os.path.getsize(self.term_bed2_file):
            #    with open(self.term_bed2_file) as infile:
            #        outfile.write(infile.read())

        options = BpToolsOptions(merged_bed2_file_tmp, self.merged_bed2_file,
                                 self.fastq_prefix, self.min_sv_size, self.max_sv_size)
        bp_tools = BpTools(options)
        bp_tools.merge_bp()

    @log(message="step 9 (validate and generate stats by altref)")
    def validate_bp_by_altref(self):
        action = 'all' if not self.disable_gen_altref_bam else None
        options = AltrefOptions(self.sample_name, '%s_altref' %
                                (self.output_prefix), self.fastq_file,
                                self.minimap2_bam_file, self.merged_bed2_file, action)
        altref = Altref(options)
        altref.run()

    @log(message="step ? (detect term SV)")
    def call_term_sv(self):
        config = self.config
        output_prefix = '%s/term_sv_caller_temp/%s' % (os.path.dirname(self.output_prefix), self.sample_name)
        term_threshold = 300000
        buf_size = 200000
        min_clip_size = 100
        term_seq_file_basepath = ''

        options = TermSvCallerOptions(
            self.sample_name,
            self.minimap2_bam_file,
            self.depth_file_filtered,
            self.fai_file,
            output_prefix,
            term_threshold,
            min_clip_size,
            self.fastq_file,
            term_seq_file_basepath,
            config['ref'],
            self.term_bed2_file,
            buf_size
        )
        term_sv_caller = TermSvCaller(options)
        self.term_sv = term_sv_caller.run()

    def filter_depth_file(self):
        filter_depth_file(self.gender, self.depth_file, self.depth_file_filtered)

    @log(message="step ? (output final result)")
    def gen_final_result(self):
        altref_result = '%s_altref_filtered.result' % (self.output_prefix)
        shutil.copy2(altref_result, self.final_result)

        # merge term sv
        with open(self.final_result, 'a') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            with open(self.term_bed2_file, 'r') as infile:
                for line in infile:
                    arr = line.strip().split()
                    chrom, start, chrom2, end, read_name, read_length = arr[0], int(
                        arr[1]), arr[2], int(arr[3]), arr[5], int(arr[6])
                    """
                    if start == 0:
                        start = 'pter'
                    if end == 0:
                        end = 'qter'
                    """
                    sv_str = '%s_%s_%s_%s_TERM-DEL' % (chrom, str(start), chrom2, str(end))
                    row = [sv_str, '-', read_name, str(read_length), '?', '?']
                    writer.writerow(row)

        final_bed2_tmp = '%s_tmp' % (self.final_bed2)
        cmd = "cat %s | tail -n +2 | cut -f1 | grep -v '^-' | sort -k1V | uniq | tr '_' '\t' > %s" % (
            self.final_result, final_bed2_tmp)
        #print('cmd', cmd)
        run_shell_cmd(cmd)

        options = BpToolsOptions(final_bed2_tmp, self.final_bed2, None, self.min_sv_size, self.max_sv_size)
        bp_tools = BpTools(options)
        bp_tools.merge_bp()

    def run(self):
        self.index_fastq()
        self.gen_bam_and_chain_file()
        if not self.disable_depth_analysis:
            self.gen_depth_file()
            self.filter_depth_file()
        self.find_bp_by_cigar_str_and_split_read()
        self.filter_by_depth_and_sv_size()
        self.refine_bp_by_dp()  # self.refine_bp_by_realign()
        self.validate_bp_by_local_realign()  # self.validate_bp_by_realign()
        if not self.disable_depth_analysis:
            self.find_bp_by_chain_and_dp()
        self.call_term_sv()
        self.merge_bed2()
        self.validate_bp_by_altref()
        self.gen_final_result()


if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-sample_name', '--sample_name', help='sample name', required=True)
    parser.add_argument('-fastq', '--fastq_file', help='fastq file', required=True)
    parser.add_argument('-output_prefix', '--output_prefix', help='output prefix', required=True)

    parser.add_argument('-min_sv_size', '--min_sv_size', help='min Sv Size', required=False, type=int)
    parser.add_argument('-max_sv_size', '--max_sv_size', help='max Sv Size', required=False, type=int)
    parser.add_argument('-disable_dp_filter', '--disable_dp_filter', help='disable DP filter', required=False, type=int)
    parser.add_argument('-disable_gen_altref_bam', '--disable_gen_altref_bam',
                        help='disable gen altref bam', required=False, type=int)
    #parser.add_argument('-sensitive_mode', help='eneable sensitive mode', nargs='?', default=False, type=bool)
    parser.add_argument('-target_sv_type', '--target_sv_type', help='target sv type', required=False)

    init_logger()

    options = parser.parse_args()
    senSV = SENSV(options)
    senSV.run()
