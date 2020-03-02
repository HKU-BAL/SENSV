import pysam
import csv
import math
from sys import exit
from os import makedirs, path
from multiprocessing import Pool
from collections import defaultdict

from utility import (
    get_var,
    run_shell_cmd,
    get_ref,
    rev_comp,
    is_same_arm,
    load_common_config,
    load_cytobands,
    get_short_name,
    get_sorted_sv_str_list,
    get_seq_from_fastq,
    init_logger,
)

"""
class Altref generates contigs of alt. ref for each of the SV candidates.
It also generates supporting information (if any) of the candidates
"""


class AltrefOptions:
    def __init__(self, name, output_prefix, fastq, orig_bam, input_bed2, action):
        self.name = name
        self.output_prefix = output_prefix
        self.fastq = fastq
        self.orig_bam = orig_bam
        self.input_bed2 = input_bed2
        self.action = action


class Altref:
    def __init__(self, options):
        self.load_config()

        try:
            makedirs(self.outputPath)
        except:
            pass

        self.output_prefix = options.output_prefix
        self.fastq = options.fastq
        self.orig_bam = options.orig_bam
        self.input_bed2 = options.input_bed2
        self.action = options.action

        output_prefix = options.output_prefix
        self.name = options.name
        self.altref = '%s.altref' % (output_prefix)
        self.localref = '%s.localref' % (output_prefix)
        self.out_ref = '%s.fa' % (output_prefix)
        self.out_bam = '%s_minimap2.bam' % (output_prefix)
        self.out_result = '%s.result' % (output_prefix)
        self.out_filtered_result = '%s_filtered.result' % (output_prefix)
        self.out_filtered_bed2 = '%s_filtered.bed2' % (output_prefix)
        self.filtered_read_fasta = '%s_filtered.fasta' % (output_prefix)

        load_cytobands()

    def load_config(self):
        config = {
            'pool_size': int(get_var('altref', 'pool_size')),

            'buf_size': int(get_var('altref', 'buf_size')),
            'ref_buf_size': int(get_var('altref', 'ref_buf_size')),
            'min_mapq': int(get_var('altref', 'min_mapq')),
            'max_indel_span': int(get_var('altref', 'max_indel_span')),
            'search_size': int(get_var('altref', 'search_size')),
            'chrom_max_size': int(get_var('altref', 'chrom_max_size')),

            'filter_max_clip': int(get_var('altref', 'filter_max_clip')),
            'filter_max_clip_ratio': float(get_var('altref', 'filter_max_clip_ratio')),
            'filter_min_avg_as': float(get_var('altref', 'filter_min_avg_as')),
            'filter_min_ref_span': float(get_var('altref', 'filter_min_ref_span')),
        }
        config.update(load_common_config())

        self.config = config

    # load input sv candidates from bed2 file
    def load_bed2(self):
        sv_str_dict = {}
        with open(self.input_bed2) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                arr = line.strip().split()
                supp_type = None
                chrom, start, chrom2, end, type = arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4]
                if len(arr) >= 6:
                    if arr[5] in ['cigar']:
                        supp_type = arr[5]
                sv_str_dict[f'{chrom}_{start}_{chrom2}_{end}_{type}'] = {
                    'supp_type': supp_type
                }

        self.sv_str_dict = sv_str_dict
        self.sv_str_list = list(sv_str_dict.keys())

    # generate contig of alt. ref based on SV candidates
    def gen_altref_file(self, sv_str=None):
        if sv_str:
            sv_str_list = [sv_str]
            fa_file = '%s_%s' % (self.altref, sv_str)
        else:
            sv_str_list = self.sv_str_list
            fa_file = self.altref

        ref_info_list = [self.get_altref_info(sv_str) for sv_str in sv_str_list]

        # gen fa
        with open(fa_file, 'w') as out_file:
            for ref_info in ref_info_list:
                for info in ref_info:
                    ref_seq_name = info['ref_seq_name']
                    ref_seq = info['ref_seq']
                    desc = info['desc']

                    print('>%s %s' % (ref_seq_name, desc), file=out_file)
                    seq_arr = [ref_seq[i: i + 70] for i in range(0, len(ref_seq), 70)]
                    for seq in seq_arr:
                        print('%s' % seq, file=out_file)

    def gen_localref_file(self, target_sv_str=None):
        if target_sv_str:
            sv_str_list = [target_sv_str]
            fa_file = '%s_%s' % (self.localref, target_sv_str)
        else:
            sv_str_list = self.sv_str_list
            fa_file = self.localref

        ref_info_list = []
        for sv_str in sv_str_list:
            ref_info = self.get_localref_info(sv_str)
            ref_info_list.append(ref_info)

        # gen localref fa
        out_file = open(fa_file, 'w')
        for ref_info in ref_info_list:
            for info in ref_info:
                ref_seq_name = info['ref_seq_name']
                ref_seq = info['ref_seq']
                desc = info['desc']

                print('>%s %s' % (ref_seq_name, desc), file=out_file)
                seq_arr = [ref_seq[i: i + 70] for i in range(0, len(ref_seq), 70)]
                for seq in seq_arr:
                    print('%s' % seq, file=out_file)

        out_file.close()

        # gen out ref
        if target_sv_str:
            altref = '%s_%s' % (self.altref, target_sv_str)
            out_ref = '%s_%s.fa' % (self.altref, target_sv_str)
        else:
            altref = self.altref
            out_ref = self.out_ref

        cmd = "cat %s %s > %s" % (self.localref, altref, out_ref)
        #print('cmd', cmd)
        run_shell_cmd(cmd)

    # genrate sequence of alt. ref.
    # format of sv_str: '1_start_1_end_type'
    def get_altref_info(self, sv_str):
        config = self.config

        bp_chrom, bp_start, bp_chrom2, bp_end, bp_type = sv_str.split('_')
        bp_start, bp_end = int(bp_start), int(bp_end)

        ref_info = []
        seq, desc = '', ''

        if bp_type in ['DEL-FROM-INS-FROM', 'DEL-FROM-INS-TO', 'DEL-TO-INS-FROM', 'DEL-TO-INS-TO']:
            if bp_type == 'DEL-FROM-INS-FROM':
                seq = get_ref(bp_chrom, bp_start-config['buf_size'], bp_start)
                seq += get_ref(bp_chrom2, bp_end+1, bp_end+1+config['buf_size'])
            elif bp_type == 'DEL-FROM-INS-TO':
                seq = get_ref(bp_chrom, bp_start-config['buf_size'], bp_start)
                seq += rev_comp(get_ref(bp_chrom2, bp_end-config['buf_size'], bp_end))
            elif bp_type == 'DEL-TO-INS-FROM':
                seq = rev_comp(get_ref(bp_chrom, bp_start, bp_start+config['buf_size']))
                seq += get_ref(bp_chrom2, bp_end+1, bp_end+1+config['buf_size'])
            elif bp_type == 'DEL-TO-INS-TO':
                seq = rev_comp(get_ref(bp_chrom, bp_start, bp_start+config['buf_size']))
                seq += rev_comp(get_ref(bp_chrom2, bp_end-config['buf_size'], bp_end))
        elif bp_type in ['DEL', 'INV', 'DUP']:
            if bp_type == 'DEL':
                seq = get_ref(bp_chrom, bp_start-1-config['buf_size'], bp_start-1)
                seq += get_ref(bp_chrom, bp_end+1, bp_end+1+config['buf_size'])
            elif bp_type == 'INV':
                # seq1
                seq = get_ref(bp_chrom, bp_start-1-config['buf_size'], bp_start-1)
                seq += rev_comp(get_ref(bp_chrom, bp_end-config['buf_size'], bp_end))

                # seq2
                seq += rev_comp(get_ref(bp_chrom, bp_start, bp_start + config['buf_size']))
                seq += get_ref(bp_chrom, bp_end, bp_end + config['buf_size'])
            elif bp_type == 'DUP':
                seq = get_ref(bp_chrom, bp_end-config['buf_size'], bp_end)
                seq += get_ref(bp_chrom, bp_start, bp_start+config['buf_size'])
            # elif bp_type == 'INS':
            #    seq = 'N' * (bp_end - bp_start)
        elif bp_type in ['TRA']:
            if is_same_arm(bp_chrom, bp_start, bp_chrom2, bp_end):
                # seq 1
                seq = get_ref(bp_chrom, bp_start-config['buf_size'], bp_start-1)
                seq += get_ref(bp_chrom2, bp_end+1, bp_end+config['buf_size'])

                # seq 2
                seq += get_ref(bp_chrom2, bp_end-config['buf_size'], bp_end-1)
                seq += get_ref(bp_chrom, bp_start+1, bp_start+config['buf_size'])
            else:
                # seq 1
                seq = get_ref(bp_chrom, bp_start-config['buf_size'], bp_start-1)
                seq += rev_comp(get_ref(bp_chrom2, bp_end-config['buf_size'], bp_end-1))

                # seq 2
                seq += rev_comp(get_ref(bp_chrom, bp_start, bp_start+config['buf_size']))
                seq += get_ref(bp_chrom2, bp_end+1, bp_end+config['buf_size'])
        else:
            print('Error. Invalid type: %s' % (bp_type))
            exit(0)

        ref_info.append({'ref_seq_name': 'altref_%s' % sv_str, 'ref_seq': seq, 'desc': desc})

        return ref_info

    def get_localref_info(self, sv_str):
        config = self.config

        bp_chrom, bp_start, _bp_chrom2, bp_end, _bp_type = sv_str.split('_')
        bp_start, bp_end = int(bp_start), int(bp_end)

        ref_info = []
        seq, desc = '', ''

        seq = get_ref(bp_chrom, bp_start-config['ref_buf_size'], bp_start+config['ref_buf_size'])
        seq += get_ref(bp_chrom, bp_end-config['ref_buf_size'], bp_end+config['ref_buf_size'])

        ref_info.append({'ref_seq_name': 'ref_%s' % sv_str, 'ref_seq': seq, 'desc': desc})

        return ref_info

    # get info of original alignment
    def get_orig(self, read_name, sv_str):
        config = self.config
        search_size = config['search_size']

        bp_chrom, bp_start, bp_chrom2, bp_end, _bp_type = sv_str.split('_')
        bp_start, bp_end = int(bp_start), int(bp_end)

        orig = {'AS': 0, 'MAPQ': 0}
        bam = pysam.AlignmentFile(self.orig_bam, 'rb')

        for chr, pos in [[bp_chrom, bp_start], [bp_chrom2, bp_end]]:
            for read in bam.fetch(chr, pos - search_size if pos > search_size else 1, pos + search_size):
                if read.is_secondary:
                    continue
                if get_short_name(read.query_name) != read_name:
                    continue

                AS = int(read.get_tag('AS'))
                if AS > orig['AS']:
                    orig = {'AS': AS, 'MAPQ': read.mapping_quality}

        return orig

    # get indel info
    def get_indel_info(self, read, span):
        max_indel_span = self.config['max_indel_span']
        cType = self.config['cType']
        cTypeConsumeQuery = self.config['cTypeConsumeQuery']
        cTypeConsumeRef = self.config['cTypeConsumeRef']

        indel_info = {
            'all': 0,
            'left': {
                'indel_count': 0,
                'indel_length': 0,
                'span': 0,
            },
            'right': {
                'indel_count': 0,
                'indel_length': 0,
                'span': 0,
            },
        }

        ref_pos = read.reference_start
        query_pos = 0
        cigar_tuples = read.cigartuples

        indel_info['left']['span'] = min(max_indel_span, span['left'])
        indel_info['right']['span'] = min(max_indel_span, span['right'])

        start = span['left'] - indel_info['left']['span']
        end = span['left'] + indel_info['right']['span']

        for type, length in cigar_tuples:
            if start <= query_pos <= end:
                side = 'left' if query_pos <= span['left'] else 'right'
                if cType[type] in 'ID':
                    indel_info[side]['indel_count'] += 1
                    indel_info[side]['indel_length'] += length

            if cType[type] in cTypeConsumeQuery:
                query_pos += length
            if cType[type] in cTypeConsumeRef:
                ref_pos += length

        indel_info['all'] = indel_info['left']['indel_length'] + indel_info['right']['indel_length']

        return indel_info

    # get alignment score info
    def get_AS(self, read, span):
        cType = self.config['cType']
        cTypeConsumeQuery = self.config['cTypeConsumeQuery']
        cTypeConsumeRef = self.config['cTypeConsumeRef']

        aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
        aligned_pairs_dict = dict((x, {'pos': y, 'seq': z}) for x, y, z in aligned_pairs)

        AS = {'all': 0, 'left': 0, 'right': 0}

        ref_pos = read.reference_start
        query_pos = 0
        cigar_tuples = read.cigartuples
        for type, length in cigar_tuples:
            side = 'left' if query_pos <= span['left'] else 'right'
            if cType[type] in 'M=':
                for i in range(length):
                    if aligned_pairs_dict[query_pos + i]['seq'] in 'ACGT':
                        AS[side] += 2
                    if aligned_pairs_dict[query_pos + i]['seq'] in 'acgt':
                        AS[side] -= 4
            elif cType[type] in 'ID':
                AS[side] -= min(4 + length * 2, 24 + length * 1)

            if cType[type] in cTypeConsumeQuery:
                query_pos += length
            if cType[type] in cTypeConsumeRef:
                ref_pos += length

        AS['all'] = AS['left'] + AS['right']

        return AS

    # generate ref file with contigs of alt. ref
    def gen_ref_file(self, sv_str=None):
        ref = self.config['ref']

        if sv_str:
            altref = '%s_%s' % (self.altref, sv_str)
            out_ref = '%s_%s.fa' % (self.altref, sv_str)
        else:
            altref = self.altref
            out_ref = self.out_ref

        cmd = "cat %s %s > %s" % (ref, altref, out_ref)
        run_shell_cmd(cmd)

    # alignment
    def align(self, sv_str=None):
        minimap2 = self.config['minimap2']
        samtools = self.config['samtools']

        if sv_str:
            out_ref = '%s_%s.fa' % (self.altref, sv_str)
            fastq = '%s_%s.fasta' % (self.altref, sv_str)
            out_bam = '%s_%s.bam' % (self.altref, sv_str)
        else:
            out_ref = self.out_ref
            fastq = self.filtered_read_fasta
            out_bam = self.out_bam

        ref_size = math.ceil(1. * path.getsize(out_ref) / 1024 / 1024 / 1024 + .5)
        if ref_size > 4:
            cmd = "%s -Y -I %dG -t 48 --MD -a %s %s | %s sort -@ 48 -o %s - && %s index -@ 48 %s" % \
                (minimap2, ref_size, out_ref, fastq, samtools, out_bam, samtools, out_bam)
        else:
            cmd = "%s -Y -t 48 --MD -a %s %s | %s sort -@ 48 -o %s - && %s index -@ 48 %s" % \
                (minimap2, out_ref, fastq, samtools, out_bam, samtools, out_bam)
        #print('cmd', cmd)
        run_shell_cmd(cmd)

    def get_mapq(self, sv_str, read_list):
        mapq = defaultdict(dict)

        out_bam = '%s_%s.bam' % (self.altref, sv_str)
        bam = pysam.AlignmentFile(out_bam, 'rb')
        for read in bam.fetch('altref_%s' % sv_str):
            if read.is_secondary:
                continue

            arr = read.query_name.split('_')
            read_name, part = arr[0], arr[1]

            mapq[read_name][part] = read.mapping_quality

        return mapq

    # get span of alignment
    def get_span(self, read, bp_start):
        cType = self.config['cType']
        cTypeConsumeQuery = self.config['cTypeConsumeQuery']
        cTypeConsumeRef = self.config['cTypeConsumeRef']

        aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
        aligned_pairs_dict = dict((x, {'pos': y, 'seq': z}) for x, y, z in aligned_pairs)

        ref_pos = read.reference_start
        cigar_tuples = read.cigartuples
        span = {
            'left': 0,
            'right': 0,
            'ref_left': 0,
            'ref_right': 0,
            'match_left': 0,
            'match_right': 0,
        }

        query_pos = 0
        for type, length in cigar_tuples:
            side = 'left' if ref_pos <= bp_start else 'right'
            if cType[type] in 'M=':
                #span['match_' + side ] += length
                for i in range(length):
                    if aligned_pairs_dict[query_pos + i]['seq'] in 'ACGT':
                        span['match_' + side] += 1

            if cType[type] in cTypeConsumeQuery:
                if cType[type] not in 'HS':
                    span[side] += length

                query_pos += length

            if cType[type] in cTypeConsumeRef:
                ref_pos += length

        span['ref_left'] = bp_start - read.reference_start
        span['ref_right'] = read.reference_end - bp_start

        return span

    def gen_sv_result(self, sv_str):
        buf_size = self.config['buf_size']
        min_mapq = self.config['min_mapq']
        cType = self.config['cType']

        _bp_chrom, bp_start, _bp_chrom2, bp_end, bp_type = sv_str.split('_')
        bp_start, bp_end = int(bp_start), int(bp_end)

        sv_length = bp_end - bp_start

        clip = {'start': 0, 'end': 0}

        sv_bps = []
        if bp_type in ['DEL', 'INS', 'DEL-FROM-INS-FROM', 'DEL-FROM-INS-TO', 'DEL-TO-INS-FROM', 'DEL-TO-INS-TO']:
            #sv_bps = [bp_start]
            sv_bps = [buf_size]
        elif bp_type in ['DUP']:
            #sv_bps = [bp_end-bp_start+buf_size]
            sv_bps = [buf_size]
        elif bp_type in ['INV', 'TRA']:
            sv_bps = [buf_size, buf_size * 3]
        else:
            print('Error. Unknown bp_type %s' % (bp_type))
            exit(0)

        supp_read = {}
        bam = pysam.AlignmentFile(self.out_bam, 'rb')

        for sv_bp in sv_bps:
            # for read in bam.fetch('altref', sv_bp-1, sv_bp):
            for read in bam.fetch('altref_%s' % sv_str, sv_bp-1, sv_bp):
                # if read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapq:
                if read.is_secondary or read.mapping_quality < min_mapq:
                    continue

                read_name = get_short_name(read.query_name)
                if read_name in supp_read:
                    continue

                # check for exceptionally high depth region => ignore them
                if len(supp_read) > 8:
                    supp_read = {}
                    break

                ref_start = read.reference_start + 1
                ref_end = read.reference_end

                orig = self.get_orig(read_name, sv_str)
                AS = {'all': read.get_tag('AS'), 'left': 0, 'right': 0}

                # if AS['all'] < orig['AS']:
                #    continue

                supp_read[read_name] = {}
                supp_read[read_name]['strand'] = 1 if read.is_reverse else 0
                supp_read[read_name]['sv_str'] = sv_str
                supp_read[read_name]['sv_length'] = sv_length

                if sv_bp == buf_size:
                    supp_read[read_name]['bp_pos'] = 'bp_start'
                else:
                    supp_read[read_name]['bp_pos'] = 'bp_end'
                supp_read[read_name]['ref_start'] = ref_start
                supp_read[read_name]['ref_end'] = ref_end
                supp_read[read_name]['clip'] = clip
                #supp_read[read_name]['MAPQ'] = read.mapping_quality
                supp_read[read_name]['AS'] = AS
                supp_read[read_name]['origAS'] = orig['AS']
                supp_read[read_name]['orig_MAPQ'] = orig['MAPQ']
                supp_read[read_name]['length'] = read.infer_read_length()
                span = self.get_span(read, sv_bp)
                supp_read[read_name]['span'] = span

                cigar_tuples = read.cigartuples
                supp_read[read_name]['clip'] = {'start': 0, 'end': 0}
                type, length = cigar_tuples[0]
                if cType[type] in 'HS':
                    supp_read[read_name]['clip']['start'] = length
                type, length = cigar_tuples[-1]
                if cType[type] in 'HS':
                    supp_read[read_name]['clip']['end'] = length

                supp_read[read_name]['AS'] = self.get_AS(read, span)
                as_inc_str = ''
                if supp_read[read_name]['origAS']:
                    as_inc_str = '{0:.0%}'.format(1.*(supp_read[read_name]['AS']['all']-supp_read[read_name]['origAS'])/supp_read[read_name]['origAS'])
                else:
                    as_inc_str = 'N/A'
                supp_read[read_name]['AS_inc'] = as_inc_str
                supp_read[read_name]['AS_avg'] = 1.*supp_read[read_name]['AS']['all']/supp_read[read_name]['length']
                supp_read[read_name]['indel'] = self.get_indel_info(read, span)

        return {
            'sv_str': sv_str,
            'supp_read': supp_read,
        }

    def __call__(self, sv_str):
        return self.gen_sv_result(sv_str)

    # generate results
    def gen_result(self):
        pool_size = self.config['pool_size']

        pool = Pool(pool_size)
        results = pool.map(self, self.sv_str_list)
        pool.close()
        pool.join()

        all_supp_read = {}
        for result in results:
            all_supp_read[result['sv_str']] = result['supp_read']

        sorted_sv_str_list = get_sorted_sv_str_list(all_supp_read)
        with open(self.out_result, 'w') as f_result:
            writer = csv.writer(f_result, delimiter='\t')

            header = [
                'sv',
                'sv_length',
                'read_name',
                'read_length',
                'origAS',
                'AS',
                'AS_inc(%)',
                'AS_left',
                'AS_right',
                'AS-to-readlen',
                'bp_start/bp_end',
                'left_clip_length',
                'right_clip_length',
                'left_query_span',
                'right_query_span',
                'left_ref_span',
                'right_ref_span',
                'left_match_span',
                'right_match_span',
                'orig_MAPQ',
                'MAPQ',
                'MAPQ_left',
                'MAPQ_right',
                'INDEL_left_length',
                'INDEL_left_count',
                'INDEL_left_span',
                'INDEL_right_length',
                'INDEL_right_count',
                'INDEL_right_span',
            ]

            writer.writerow(header)

            # for sv_str in all_supp_read:
            for sv_str in sorted_sv_str_list:

                #print('sv_str', sv_str)

                supp_read = all_supp_read[sv_str]
                mapq = self.gen_mapq(sv_str, supp_read)

                for read_name in supp_read:
                    read = supp_read[read_name]

                    for part in ['all', 'left', 'right']:
                        if part in mapq[read_name]:
                            supp_read[read_name]['MAPQ_%s' % (part)] = mapq[read_name][part]
                        else:
                            supp_read[read_name]['MAPQ_%s' % (part)] = '-'

                    output = [
                        sv_str,
                        str(read['sv_length']),
                        read_name,
                        str(read['length']),
                        str(read['origAS']),
                        str(read['AS']['all']),
                        str(read['AS_inc']),
                        str(read['AS']['left']),
                        str(read['AS']['right']),
                        '{:.2f}'.format(read['AS_avg']),
                        read['bp_pos'],
                        str(read['clip']['start']),
                        str(read['clip']['end']),
                        str(read['span']['left']),
                        str(read['span']['right']),
                        str(read['span']['ref_left']),
                        str(read['span']['ref_right']),
                        str(read['span']['match_left']),
                        str(read['span']['match_right']),
                        str(read['orig_MAPQ']),
                        str(read['MAPQ_all']),
                        str(read['MAPQ_left']),
                        str(read['MAPQ_right']),
                        str(read['indel']['left']['indel_length']),
                        str(read['indel']['left']['indel_count']),
                        str(read['indel']['left']['span']),
                        str(read['indel']['right']['indel_length']),
                        str(read['indel']['right']['indel_count']),
                        str(read['indel']['right']['span']),
                    ]

                    writer.writerow(output)

                    sv_str = '-'

        self.all_supp_read = all_supp_read

    def gen_fasta(self, sv_str, read_info_list):
        fasta_file = '%s_%s.fasta' % (self.altref, sv_str)

        with open(fasta_file, 'w') as f:
            for read_name in read_info_list:
                read_info = read_info_list[read_name]
                span = read_info['span']

                seqs = {}
                seqs['all'] = get_seq_from_fastq(read_name, self.fastq, '')
                if read_info['strand'] == 1:
                    seqs['all'] = rev_comp(seqs['all'])
                seqs['left'] = seqs['all'][0:span['left']]
                seqs['right'] = seqs['all'][span['left']:]

                for part in seqs:
                    seq = seqs[part]

                    print('>%s_%s' % (read_name, part), file=f)
                    print('%s' % seq, file=f)

    def gen_mapq(self, sv_str, supp_read):
        mapq = {}

        for read_name in supp_read:
            mapq[read_name] = {'all': '-', 'left': '-', 'right': '-'}

        # problem: take long time, as 1 alignment for 1 candidate
        # commented now. to be fixed
        """
        # gen realign fa
        self.gen_altref_file(sv_str)
        self.gen_ref_file(sv_str)

        # gen fasta
        self.gen_fasta(sv_str, supp_read)

        # realign
        self.align(sv_str)

        # get mapq
        mapq = self.get_mapq(sv_str, supp_read.keys())
        """

        return mapq

    # generate resulting, filtered bed2 file
    def gen_filtered_bed2(self):
        with open(self.out_filtered_result, 'w') as f_filtered_result, open(self.out_filtered_bed2, 'w') as f_filtered_bed2:
            writer = csv.writer(f_filtered_result, delimiter='\t')
            writer_bed2 = csv.writer(f_filtered_bed2, delimiter='\t')

            header = [
                'sv',
                'sv_length',
                'read_name',
                'length',
                'origAS',
                'AS',
                'AS_inc(%)',
                'AS_left',
                'AS_right',
                'AS-to-readlen',
                'bp_start/bp_end',
                'left_clip_length',
                'right_clip_length',
                'left_query_span',
                'right_query_span',
                'left_ref_span',
                'right_ref_span',
                'left_match_span',
                'right_match_span',
                'orig_MAPQ',
                'MAPQ',
                'MAPQ_left',
                'MAPQ_right',
                'INDEL_left_length',
                'INDEL_left_count',
                'INDEL_left_span',
                'INDEL_right_length',
                'INDEL_right_count',
                'INDEL_right_span',
            ]

            writer.writerow(header)

            sorted_sv_str_list = get_sorted_sv_str_list(self.all_supp_read)
            for sv_str in sorted_sv_str_list:
                supp_read = self.all_supp_read[sv_str]

                # filtered = 0
                # n_supp = len(supp_read)

                for read_name in supp_read:
                    read = supp_read[read_name]

                    # filter
                    # if read['clip']['start'] > config['filter_max_clip'] or read['clip']['end'] > config['filter_max_clip'] or \
                    """
                    if 1. * read['clip']['start'] / read['length'] > config['filter_max_clip_ratio'] or \
                        1. * read['clip']['end'] / read['length'] > config['filter_max_clip_ratio'] or \
                        read['AS_avg'] < config['filter_min_avg_as']:
                    """
                    """
                    if 1. * read['clip']['start'] / read['length'] > config['filter_max_clip_ratio'] or \
                        1. * read['clip']['end'] / read['length'] > config['filter_max_clip_ratio']:
                        continue
                    """

                    """
                    if n_supp <= 1 and self.sv_str_dict[sv_str]['supp_type'] not in ['cigar'] and \
                        (1. * read['clip']['start'] / read['length'] > config['filter_max_clip_ratio'] or \
                        1. * read['clip']['end'] / read['length'] > config['filter_max_clip_ratio'] or \
                        read['AS_avg'] < config['filter_min_avg_as']):
                        continue
                    """

                    """
                    if read['span']['ref_left'] < config['filter_min_ref_span'] or read['span']['ref_right'] < config['filter_min_ref_span']:
                        continue
                    """

                    """
                    if read['clip']['start'] > config['filter_max_clip'] or read['clip']['end'] > config['filter_max_clip']:
                        continue
                    """
                    output = [
                        sv_str,
                        str(read['sv_length']),
                        read_name,
                        str(read['length']),
                        str(read['origAS']),
                        str(read['AS']['all']),
                        str(read['AS_inc']),
                        str(read['AS']['left']),
                        str(read['AS']['right']),
                        '{:.2f}'.format(read['AS_avg']),
                        read['bp_pos'],
                        str(read['clip']['start']),
                        str(read['clip']['end']),
                        str(read['span']['left']),
                        str(read['span']['right']),
                        str(read['span']['ref_left']),
                        str(read['span']['ref_right']),
                        str(read['span']['match_left']),
                        str(read['span']['match_right']),
                        str(read['orig_MAPQ']),
                        str(read['MAPQ_all']),
                        str(read['MAPQ_left']),
                        str(read['MAPQ_right']),
                        str(read['indel']['left']['indel_length']),
                        str(read['indel']['left']['indel_count']),
                        str(read['indel']['left']['span']),
                        str(read['indel']['right']['indel_length']),
                        str(read['indel']['right']['indel_count']),
                        str(read['indel']['right']['span']),
                    ]

                    writer.writerow(output)

                    if sv_str != '-':
                        arr = sv_str.split('_')
                        output = [arr[0], arr[1], arr[2], arr[3], arr[4]]
                        writer_bed2.writerow(output)

                    sv_str = '-'

    def gen_filtered_read_fasta(self):
        search_size = self.config['search_size']
        samtools = self.config['samtools']

        # gen bed from bed2
        bed = '%s.bed' % (self.output_prefix)
        temp_bam = '%s_temp.bam' % (self.output_prefix)

        cmd = "python bed2tobed.py -in_bed2 %s -out_bed %s -search_size %s" % (self.input_bed2, bed, search_size)
        #print('cmd', cmd)
        run_shell_cmd(cmd)

        # gen temp bam
        cmd = "%s view -@ 48 -L %s %s -b -M > %s && %s index -@ 48 %s" % (samtools, bed, self.orig_bam, temp_bam, samtools, temp_bam)
        #print('cmd', cmd)
        run_shell_cmd(cmd)

        # gen fastq
        cmd = "python bam2fasta.py -input_bam %s -input_fastq %s -output_fasta %s" % (temp_bam, self.fastq, self.filtered_read_fasta)
        #print('cmd', cmd)
        run_shell_cmd(cmd)

    def run(self):
        if self.input_bed2:
            self.load_bed2()

        if self.action in ['all', 'local_validate']:
            self.gen_filtered_read_fasta()
            self.gen_altref_file()

            if self.action == 'all':
                self.gen_ref_file()
            elif self.action == 'local_validate':
                self.gen_localref_file()

            self.align()
            self.gen_result()
            self.gen_filtered_bed2()
        else:
            self.gen_result()
            self.gen_filtered_bed2()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='run')
    parser.add_argument('-name', '--name', help='name', required=True)
    parser.add_argument('-outputprefix', '--output_prefix', help='output prefix', required=True)
    parser.add_argument('-origbam', '--orig_bam', help='orig bam', required=True)
    parser.add_argument('-inputbed2', '--input_bed2', help='input bed2 file', required=True)

    parser.add_argument('-fastq', '--fastq', help='fastq', required=False)
    parser.add_argument('-action', '--action', help='action [all]', required=False)

    init_logger()

    options = parser.parse_args()
    altref = Altref(options)
    altref.run()
