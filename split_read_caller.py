import csv
import pysam
import re
from multiprocessing import Pool
from argparse import ArgumentParser

from utility import (
    get_var,
    get_short_name,
    is_same_arm,
    load_cytobands,
    load_common_config,
    cigar_string_to_tuples,
    get_compact_cigar_string,
)

"""
class SplitReadSv finds SV candidates of all chromosome based on the cigar string and split reads.
"""


class SplitReadCallerOptions:
    def __init__(self, bam_file, chrom_list, min_sv_size, max_sv_size, target_sv_type, out_bed2, depth_bed2, term_threshold, lt_10k_bed2, fai_file):
        self.bam_file = bam_file
        self.chrom_list = chrom_list
        self.min_sv_size = int(min_sv_size)
        self.max_sv_size = int(max_sv_size)
        self.target_sv_type = target_sv_type
        self.out_bed2 = out_bed2
        self.depth_bed2 = depth_bed2
        self.term_threshold = term_threshold
        self.lt_10k_bed2 = lt_10k_bed2
        self.fai_file = fai_file


class SplitReadCaller:
    def __init__(self, options):
        self.load_config()

        self.options = options

        load_cytobands()

    def load_config(self):
        config = {}

        config['cigarDelFragmentThreshold'] = int(get_var('denovo_common', 'cigarDelFragmentThreshold'))
        config['cigarDelFragmentMinDist'] = int(get_var('denovo_common', 'cigarDelFragmentMinDist'))
        config['cigarDelThreshold'] = int(get_var('denovo_common', 'cigarDelThreshold'))

        config['cigarInsFragmentThreshold'] = int(get_var('denovo_common', 'cigarInsFragmentThreshold'))
        config['cigarInsFragmentMinDist'] = int(get_var('denovo_common', 'cigarInsFragmentMinDist'))
        config['cigarInsThreshold'] = int(get_var('denovo_common', 'cigarInsThreshold'))

        config['clipThreshold'] = int(get_var('denovo_common', 'clipThreshold'))
        config['split_gap_threshold'] = int(get_var('split_read_caller', 'split_gap_threshold'))
        #config['lt_10k_min_sv_size'] = int(get_var('split_read_caller', 'lt_10k_min_sv_size'))

        config['bp_buf_size'] = int(get_var('bp_tools', 'buf_size'))

        config.update(load_common_config())

        self.config = config

    def get_cigar_info(self, read, chrom):
        config = self.config

        info_list = []
        indel_list = []

        tuples = read.cigartuples
        ref_pos = read.reference_start
        query_pos = 0

        read_strand = '-' if read.is_reverse else '+'
        read_length = read.infer_read_length()

        for (type, length) in tuples:
            if config['cType'][type] == 'D' and length > config['cigarDelFragmentThreshold']:
                bp_query_pos = query_pos if read_strand == '+' else read_length - query_pos
                info_list.append({
                    'bp_chrom': chrom,
                    'bp_chrom2': chrom,
                    'bp_start': ref_pos,
                    'bp_end': ref_pos+length,
                    'read_name': read.query_name,
                    'supp_type': 'cigar',
                    'bpType': 'DEL',
                    'read_strand': read_strand,
                    'bp_query_pos': bp_query_pos,
                    'read_seq': read.query_sequence
                })

            if config['cType'][type] == 'I' and length > config['cigarInsFragmentThreshold']:
                bp_query_pos = query_pos if read_strand == '+' else read_length - query_pos
                info_list.append({
                    'bp_chrom': chrom,
                    'bp_chrom2': chrom,
                    'bp_start': ref_pos,
                    'bp_end': ref_pos+length,
                    'read_name': read.query_name,
                    'supp_type': 'cigar',
                    'bpType': 'INS',
                    'read_strand': read_strand,
                    'bp_query_pos': bp_query_pos,
                    'read_seq': read.query_sequence
                })

            if config['cType'][type] in config['cTypeConsumeRef']:
                ref_pos += length
            if config['cType'][type] in config['cTypeConsumeQuery']:
                query_pos += length

        #merge and output
        # assume there is only 1 target del/ins covered by a read (?)
        curr_bp_type = None
        curr_bp_start = None
        curr_bp_end = None
        for info in info_list:
            if not curr_bp_type:
                curr_bp_type = info['bpType']
            elif info['bpType'] != curr_bp_type:
                # maybe a better way to handle
                continue

            # check if 1st fragment
            if not curr_bp_start:
                curr_bp_start = info['bp_start']
                curr_bp_end = info['bp_end']
                bp_query_pos = info['bp_query_pos']
            else:
                if info['bp_start']-curr_bp_end < config['cigarDelFragmentMinDist']:
                    curr_bp_end = info['bp_end']
                    bp_query_pos = info['bp_query_pos']
                else:
                    # more than 1 target del/ins
                    # maybe a better way to handle
                    pass

        if curr_bp_type == 'DEL':
            if curr_bp_start and curr_bp_end-curr_bp_start > config['cigarDelThreshold']:
                indel_list.append({
                    'bp_chrom': chrom,
                    'bp_chrom2': chrom,
                    'bp_start': curr_bp_start,
                    'bp_end': curr_bp_end,
                    'read_name': read.query_name,
                    'supp_type': 'cigar',
                    'sv_type': curr_bp_type,
                    'read_strand': read_strand,
                    'bp_query_pos': bp_query_pos,
                    'read_seq': read.query_sequence
                })
        elif curr_bp_type == 'INS':
            if curr_bp_start and curr_bp_end-curr_bp_start > config['cigarDelThreshold']:
                indel_list.append({
                    'bp_chrom': chrom,
                    'bp_chrom2': chrom,
                    'bp_start': curr_bp_start,
                    'bp_end': curr_bp_end,
                    'read_name': read.query_name,
                    'supp_type': 'cigar',
                    'sv_type': curr_bp_type,
                    'read_strand': read_strand,
                    'bp_query_pos': bp_query_pos,
                    'read_seq': read.query_sequence
                })

        return indel_list

    def find_sv_region(self, chrom):
        sv_region_list = []
        covered_reads = {}

        bam = pysam.AlignmentFile(self.options.bam_file, "rb")
        for read in bam.fetch(chrom):
            if read.is_secondary:
                continue

            # if read.is_supplementary:
            #    continue

            if read.query_name in covered_reads:
                continue

            covered_reads[read.query_name] = 1

            bp_infoList = self.get_bp_info(read, chrom)

            for bp_info in bp_infoList:
                if 'ALL' not in self.options.target_sv_type and bp_info['sv_type'] not in self.options.target_sv_type:
                    continue

                if bp_info['sv_type'] not in ['TRA', 'TRA_INV']:
                    bp_length = bp_info['bp_end']-bp_info['bp_start']

                    if self.options.min_sv_size and bp_length + (bp_info['overlap'] if "overlap" in bp_info else 0) < self.options.min_sv_size:
                        continue

                    if self.options.max_sv_size and bp_length > self.options.max_sv_size:
                        continue

                sv_region_list.append(bp_info)

        return sv_region_list

    def get_query_pos(self, orig_tuples, strand, total_read_length):
        start = self.get_query_start_pos(orig_tuples, strand)
        end = self.get_query_end_pos(orig_tuples, strand, total_read_length)

        return {'query_start': start, 'query_end': end}

    def get_split_info_detail(self, chrom, read_name, read_length, supp_list, read):
        config = self.config
        supp_list.sort(key=lambda x: int(x.split(',')[1]))

        info = []
        split_gap_threshold = config['split_gap_threshold']
        bp_chrom2 = None

        prev = {}
        for suppStr in supp_list:
            bp_chrom2 = None

            arr = suppStr.split(',')
            supp_chrom = arr[0]
            supp_ref_start = int(arr[1])
            supp_ref_end = int(arr[2])
            supp_strand = arr[3]
            supp_cigar_string = arr[4]

            cigarRe = r"(\d+)([MIDNSHP=X])"
            cigar_tuple = []
            for m in re.finditer(cigarRe, supp_cigar_string):
                cigar_tuple.append((m.group(1), m.group(2)))
            supp_for_clipped = cigar_tuple[0][0] if cigar_tuple[0][1] == "S" else 0
            supp_back_clipped = cigar_tuple[-1][0] if cigar_tuple[-1][1] == "S" else 0
            supp_clipped = len(read.query_sequence)-max(int(supp_for_clipped), int(supp_back_clipped))

            if supp_chrom not in ['X', 'Y']:
                try:
                    int(supp_chrom)
                except ValueError:
                    continue

            tuples = cigar_string_to_tuples(supp_cigar_string)
            ref_length = 0
            for tuple in tuples:
                if tuple[0] in config['cTypeConsumeRef']:
                    ref_length += int(tuple[1])

            supp_ref_end = supp_ref_start + ref_length

            supp_query_pos = self.get_query_pos(tuples, supp_strand, read_length)

            curr = {
                'ref_chrom': supp_chrom,
                'ref_start': supp_ref_start,
                'ref_end': supp_ref_end,
                'query_start': supp_query_pos['query_start'],
                'query_end': supp_query_pos['query_end'],
                'strand': supp_strand,
                'clipped': supp_clipped
            }

            if not prev:
                prev = curr
                continue

            detect_type = None
            sv_type = None

            # TRA
            if curr['ref_chrom'] != prev['ref_chrom']:
                if prev['ref_chrom'] == chrom:
                    bp1, bp2 = prev, curr
                elif curr['ref_chrom'] == chrom:
                    bp1, bp2 = curr, prev
                else:
                    continue

                bp_start, bp_end, bp_start_query, bp_end_query = None, None, None, None

                if bp1['strand'] == bp2['strand']:
                    # simple
                    # case 1: (+)->| ... |->(+)
                    if bp1['strand'] == '+' and bp1['query_start'] < bp2['query_start']:
                        bp_start, bp_start_query = bp1['ref_end'], bp1['query_end']
                        bp_end, bp_end_query = bp2['ref_start'], bp2['query_start']
                        detect_type = 'TRA'

                    # case 2: (-)<-| ... |<-(-)
                    elif bp1['strand'] == '-' and bp1['query_start'] > bp2['query_start']:
                        bp_start, bp_start_query = bp1['ref_end'], bp1['query_start']
                        bp_end, bp_end_query = bp2['ref_start'], bp2['query_end']
                        detect_type = 'TRA'

                    # case 3: |<-(-) ... (-)<-|
                    elif bp1['strand'] == '-' and bp1['query_start'] < bp2['query_start']:
                        bp_start, bp_start_query = bp1['ref_start'], bp1['query_end']
                        bp_end, bp_end_query = bp2['ref_end'], bp2['query_start']
                        detect_type = 'TRA'

                    # case 4: |->(+) ... (+)->|
                    elif bp1['strand'] == '+' and bp1['query_start'] > bp2['query_start']:
                        bp_start, bp_start_query = bp1['ref_start'], bp1['query_start']
                        bp_end, bp_end_query = bp2['ref_end'], bp2['query_end']
                        detect_type = 'TRA'

                    if detect_type:
                        # check for arm
                        if not is_same_arm(bp1['ref_chrom'], bp_start, bp2['ref_chrom'], bp_end):
                            continue

                        query_diff = abs(bp_start_query - bp_end_query)

                        if query_diff < split_gap_threshold:
                            sv_type = 'TRA'
                            bp_chrom2 = bp2['ref_chrom']
                else:
                    bp_start, bp_end, bp_start_query, bp_end_query = None, None, None, None
                    detect_type = None

                    # TRA + INV
                    if not detect_type:
                        # case 1: (+)->| ... (-)<-|
                        if bp1['strand'] == '+' and bp1['query_start'] < bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_end'], bp1['query_end']
                            bp_end, bp_end_query = bp2['ref_end'], bp2['query_start']
                            detect_type = 'TRA_INV'

                        # case 2: (-)<-| ... (+)->|
                        elif bp1['strand'] == '-' and bp1['query_start'] > bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_end'], bp1['query_start']
                            bp_end, bp_end_query = bp2['ref_end'], bp2['query_end']
                            detect_type = 'TRA_INV'

                        # case 3:      |<-(-) ... |->(+)
                        elif bp1['strand'] == '-' and bp1['query_start'] < bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_start'], bp1['query_end']
                            bp_end, bp_end_query = bp2['ref_start'], bp2['query_start']
                            detect_type = 'TRA_INV'

                        # case 4:      |->(+) ... |<-(-)
                        elif bp1['strand'] == '+' and bp1['query_start'] > bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_start'], bp1['query_start']
                            bp_end, bp_end_query = bp2['ref_start'], bp2['query_end']
                            detect_type = 'TRA_INV'

                        if detect_type:
                            # check for arm
                            if is_same_arm(bp1['ref_chrom'], bp_start, bp2['ref_chrom'], bp_end):
                                continue

                            query_diff = abs(bp_start_query - bp_end_query)
                            # if query_diff < split_gap_threshold:
                            if True:
                                sv_type = 'TRA_INV'
                                bp_chrom2 = bp2['ref_chrom']

            # other types
            # else:
            elif prev['ref_chrom'] == chrom and curr['ref_chrom'] == chrom:
                if prev['ref_end'] < curr['ref_start']:
                    bp1, bp2 = prev, curr
                else:
                    bp1, bp2 = curr, prev

                if bp1['strand'] == bp2['strand']:
                    bp_start, bp_end, bp_start_query, bp_end_query = None, None, None, None
                    detect_type = None

                    # DUP
                    if not detect_type:
                        # case 1: |->(+) ... (+)->|
                        if bp1['strand'] == '+' and bp1['query_start'] > bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_start'], bp1['query_start']
                            bp_end, bp_end_query = bp2['ref_end'], bp2['query_end']
                            detect_type = 'DUP'

                        # case 2: |<-(-) ... (-)<-|
                        elif bp1['strand'] == '-' and bp1['query_start'] < bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_start'], bp1['query_end']
                            bp_end, bp_end_query = bp2['ref_end'], bp2['query_start']
                            detect_type = 'DUP'

                        if detect_type:
                            query_diff = abs(bp_start_query - bp_end_query)
                            if query_diff < split_gap_threshold:
                                sv_type = 'DUP'

                    # INS/DEL
                    if not detect_type:
                        # case 1: (+)->| ... |->(+)
                        if bp1['strand'] == '+' and bp1['query_start'] < bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_end'], bp1['query_end']
                            #bp_chrom2, bp_end, bp_end_query = bp2['ref_chrom'], bp2['ref_start'], bp2['query_start']
                            bp_end, bp_end_query = bp2['ref_start'], bp2['query_start']
                            detect_type = 'INDEL'
                        # case 1: (-)<-| ... |<-(-)
                        elif bp1['strand'] == '-' and bp1['query_start'] > bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_end'], bp1['query_start']
                            #bp_chrom2, bp_end, bp_end_query = bp2['ref_chrom'], bp2['ref_start'], bp2['query_end']
                            bp_end, bp_end_query = bp2['ref_start'], bp2['query_end']
                            detect_type = 'INDEL'

                        if detect_type:
                            # INS
                            ref_diff = abs(bp_start - bp_end)
                            if ref_diff < split_gap_threshold:
                                sv_type = 'INS'

                            # DEL
                            query_diff = abs(bp_start_query - bp_end_query)

                            if query_diff < split_gap_threshold:
                                # if True:
                                sv_type = 'DEL'
                                #bp_end = bp_end + query_diff

                elif bp1['strand'] != bp2['strand']:
                    bp_start, bp_end, bp_start_query, bp_end_query = None, None, None, None
                    detect_type = None

                    # INV
                    if not detect_type:
                        # case 1: (+)->| ... (-)<-|
                        if bp1['strand'] == '+' and bp1['query_start'] < bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_end'], bp1['query_end']
                            bp_end, bp_end_query = bp2['ref_end'], bp2['query_start']
                            detect_type = 'INV'

                        # case 2: (-)<-| ... (+)->|
                        elif bp1['strand'] == '-' and bp1['query_start'] > bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_end'], bp1['query_start']
                            bp_end, bp_end_query = bp2['ref_end'], bp2['query_end']
                            detect_type = 'INV'

                        # case 3:      |<-(-) ... |->(+)
                        elif bp1['strand'] == '-' and bp1['query_start'] < bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_start'], bp1['query_end']
                            bp_end, bp_end_query = bp2['ref_start'], bp2['query_start']
                            detect_type = 'INV'

                        # case 4:      |->(+) ... |<-(-)
                        elif bp1['strand'] == '+' and bp1['query_start'] > bp2['query_start']:
                            bp_start, bp_start_query = bp1['ref_start'], bp1['query_start']
                            bp_end, bp_end_query = bp2['ref_start'], bp2['query_end']
                            detect_type = 'INV'

                        if detect_type:
                            query_diff = abs(bp_start_query - bp_end_query)
                            if query_diff < split_gap_threshold:
                                sv_type = 'INV'

            if sv_type:
                if not bp_chrom2:
                    bp_chrom2 = chrom

                if bp1['query_start'] < bp2['query_start']:
                    read_strand = bp1['strand']
                    bp_query_pos = bp1['query_end']
                else:
                    read_strand = bp2['strand']
                    bp_query_pos = bp2['query_end']

                info.append({'bp_chrom': chrom, 'bp_start': bp_start, 'bp_end': bp_end, 'read_name': read_name,
                             'supp_type': 'splitRead', 'sv_type': sv_type, 'bp_chrom2': bp_chrom2, 'read_strand': read_strand, 'bp_query_pos': bp_query_pos, 'read_seq': read.query_sequence, 'overlap': max(0, curr['clipped']+prev['clipped']-len(read.query_sequence))})

            prev = curr

        return info

    def get_split_info(self, read, read_chrom):
        split_info = []
        supp_list = []

        read_name = read.query_name
        read_length = read.infer_read_length()

        if read.has_tag("SA"):
            sa_list = read.get_tag("SA").strip(';').split(';')
            for sa in sa_list:
                arr = sa.split(',')
                chrom = arr[0]
                ref_start = int(arr[1])
                strand = arr[2]
                cigar_string = arr[3]

                pos_info = self.get_pos_info(cigar_string, strand, ref_start)
                ref_end = pos_info['ref_end']
                info = '%s,%s,%s,%s,%s,%s,%s' % (chrom, ref_start, ref_end, strand, cigar_string,
                                                 pos_info['query_start'], pos_info['query_end'])
                supp_list.append(info)

            # add current read alignment
            cigar_string = get_compact_cigar_string(read.cigarstring, 0, 0)
            strand = '-' if read.is_reverse else '+'
            ref_start = read.reference_start

            pos_info = self.get_pos_info(cigar_string, strand, ref_start)
            ref_end = pos_info['ref_end']
            info = '%s,%s,%s,%s,%s,%s,%s' % (read_chrom, ref_start, ref_end, strand, cigar_string,
                                             pos_info['query_start'], pos_info['query_end'])
            supp_list.append(info)

            split_info = self.get_split_info_detail(read_chrom, read_name, read_length, supp_list, read)

        return split_info

    def get_query_start_pos(self, orig_tuples, strand):
        # equiv. to count hardclip/softclip at the beginning

        if strand == '-':
            tuples = orig_tuples[::-1]
        else:
            tuples = orig_tuples

        pos = 0
        for (type, length) in tuples:
            # if cType[type] in 'HS':
            if type in 'HS':
                pos += int(length)
            else:
                break

        return pos

    def get_query_end_pos(self, orig_tuples, strand, total_read_length):
        # equiv. to count querystring up to ending softclip/hardclip

        if strand == '-':
            tuples = orig_tuples[::-1]
        else:
            tuples = orig_tuples

        end_clip_count = 0
        rev_tuples = tuples[::-1]
        for (type, length) in rev_tuples:
            if type in 'HS':
                end_clip_count += int(length)
            else:
                break

        return total_read_length - end_clip_count - 1

    def get_pos_info(self, cigar_string, strand, ref_start):
        config = self.config

        read_length = 0
        ref_length = 0

        tuples = cigar_string_to_tuples(cigar_string)
        for (type, length) in tuples:
            if type in config['cTypeConsumeRef']:
                ref_length += length
            # if type in cTypeConsumeQuery:
            if type in config['cTypeConsumeRead']:
                read_length += length

        ref_end = ref_start + ref_length
        query_start = self.get_query_start_pos(tuples, strand)
        query_end = self.get_query_end_pos(tuples, strand, read_length)

        output = {
            'read_length': read_length,
            'ref_length': ref_length,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'query_start': query_start,
            'query_end': query_end,
        }

        return output

    def get_bp_info(self, read, chrom):
        return self.get_cigar_info(read, chrom) + self.get_split_info(read, chrom)

    def __call__(self, chrom):
        return {
            'sv_region_list': self.find_sv_region(chrom)
        }

    def gen_bed2(self, results):
        options = self.options
        config = self.config

        output = {}
        output_lt_10k = {}
        lt_10k_max_sv_size = 10000

        with open(self.options.out_bed2, 'w') as f, open(self.options.lt_10k_bed2, 'w') as f_lt_10k:
            for result in results:
                for sv_region in result['sv_region_list']:
                    chrom = sv_region['bp_chrom']
                    start = str(sv_region['bp_start'])
                    chrom2 = sv_region['bp_chrom2']
                    end = str(sv_region['bp_end'])
                    sv_type = sv_region['sv_type']
                    overlap = sv_region['overlap'] if 'overlap' in sv_region else 0

                    read_name = get_short_name(sv_region['read_name'])
                    # read_seq = sv_region['read_seq']
                    read_strand = sv_region['read_strand']
                    bp_query_pos = str(sv_region['bp_query_pos'])

                    supp_type = 'unknown'
                    if 'supp_type' in sv_region:
                        supp_type = sv_region['supp_type']

                    #line = [chrom, start, chrom2, end, sv_type, supp_type, read_name, read_seq]
                    line = [chrom, start, chrom2, end, sv_type, supp_type, read_name, read_strand, bp_query_pos]

                    # if sv_type in ['DEL','DUP']  and int(end) - int(start) < options.min_sv_size:
                    #    continue

                    if sv_type in ['DEL', 'DUP'] and int(end) - int(start) < lt_10k_max_sv_size:
                        output_lt_10k['|'.join(line)] = 1

                        if int(end) - int(start) + overlap >= options.min_sv_size - config['bp_buf_size']:
                            output['|'.join(line)] = 1
                    else:
                        output['|'.join(line)] = 1

            # sort output
            toIntStr = lambda text: text if text.isdigit() else '23' if text == 'X' else '24'
            output = sorted(output, key=lambda k: toIntStr(k.split('|')[0]).rjust(2) + k.split('|')[1].rjust(9))
            output_lt_10k = sorted(output_lt_10k, key=lambda k: toIntStr(k.split('|')[0]).rjust(2) + k.split('|')[1].rjust(9))

            writer = csv.writer(f, delimiter='\t')
            for line in output:
                writer.writerow(line.split('|'))

            writer_lt_10k = csv.writer(f_lt_10k, delimiter='\t')
            for line in output_lt_10k:
                writer_lt_10k.writerow(line.split('|'))

        print('gen_bed2: >10k, count=%d' % len(output))
        print('gen_bed2: <=10k, count=%d' % len(output_lt_10k))

    def load_depth_region(self):
        options = self.options

        #depth_buf_size = options.depth_buf_size
        depth_buf_size = 600000
        depth_region = {'DEL': [], 'DUP': []}
        with open(options.depth_bed2) as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                arr = line.strip().split()
                chrom, start, end, region_type = arr[0], int(float(arr[1])), int(float(arr[2])), arr[4]

                # consider DEL and DUP regions only
                if region_type not in ['DEL', 'DUP']:
                    continue

                chrom = 'X' if chrom == '23' else 'Y' if chrom == '24' else chrom
                start1, start2, end1, end2 = start - depth_buf_size, start + depth_buf_size, end - depth_buf_size, end + depth_buf_size
                start1 = 1 if start1 <= 0 else start1
                end1 = 1 if end1 <= 0 else end1

                depth_region[region_type].append({'chrom':chrom,'start':start,'end':end,'start1':start1,'start2':start2,'end1':end1,'end2':end2})

        self.depth_region = depth_region

        #print('depth_region', self.depth_region)

    def is_sa_in_region(self, sa, region_list):
        for region in region_list:
            if (
                sa['chrom'] == region['chrom'] and
                (
                    region['start1'] <= sa['start'] <= region['start2'] or
                    region['end1'] <= sa['start'] <= region['end2']
                )
            ):
                return True

        return False

    def get_sa_ref_end(self, sa):
        config = self.config

        sa_end = sa['start']
        sa_cigar_tuples = cigar_string_to_tuples(sa['cigar_string'])
        for tuple in sa_cigar_tuples:
            if tuple[0] in config['cTypeConsumeRef']:
                sa_end += int(tuple[1])

        return sa_end

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

    def call_indel(self):
        config = self.config
        options = self.options
        indel_dict = {}

        print('call_indel_working')

        for region in self.depth_region['DEL']:
            #target_bp_regions = [[region['start1'],region['start2']],[region['end1'],region['end2']]]
            #target_bp_regions = [[region['start1'],region['start2']]]

            target_bp_regions = []

            if region['start1'] > options.term_threshold:
                target_bp_regions.append([region['start1'], region['start2']])

            chrom_size = self.chrom_info[region['chrom']]
            if region['end2'] < chrom_size - options.term_threshold:
                target_bp_regions.append([region['end1'], region['end2']])

            for bp in target_bp_regions:
                bam = pysam.AlignmentFile(options.bam_file, 'rb')
                for read in bam.fetch(region['chrom'], bp[0], bp[1]):
                    if read.is_secondary:
                        continue

                    if not read.has_tag('SA'):
                        continue

                    # del_side = ''
                    # if bp[0] == region['start1']:
                    #     del_side = 'START'
                    # elif bp[0] == region['end1']:
                    #     del_side = 'END'

                    read_name = get_short_name(read.query_name)
                    read_length = read.infer_read_length()
                    read_strand = '-' if read.is_reverse else '+'
                    cigar_tuples = cigar_string_to_tuples(read.cigarstring)
                    query_pos1 = self.get_query_pos(cigar_tuples, read_strand, read_length)

                    for sa_str in read.get_tag('SA').strip(';').split(';'):
                        arr = sa_str.split(',')
                        sa = {'chrom': arr[0], 'start': int(arr[1]), 'strand': arr[2], 'cigar_string': arr[3]}

                        #ins_info = self.get_ins_info(sa, self.depth_region['DUP'])
                        # if ins_info['ins_side']:
                        if self.is_sa_in_region(sa, self.depth_region['DUP']):
                            sa['end'] = self.get_sa_ref_end(sa)

                            query_pos2 = self.get_query_pos(cigar_string_to_tuples(sa['cigar_string']), sa['strand'], read_length)

                            del_chrom, ins_chrom = region['chrom'], sa['chrom']
                            del_from_pos, del_to_pos = 0, 0
                            ins_from_pos, ins_to_pos = 0, 0

                            if query_pos1['query_start'] < query_pos2['query_start']:
                                del_read_part, ins_read_part = 1, 2
                            else:
                                del_read_part, ins_read_part = 2, 1

                            if read_strand == '+':
                                if del_read_part == 1:
                                    del_from_pos = read.reference_end
                                else:
                                    del_to_pos = read.reference_start
                            else:
                                if del_read_part == 1:
                                    del_to_pos = read.reference_start
                                else:
                                    del_from_pos = read.reference_end

                            if sa['strand'] == '+':
                                if ins_read_part == 1:
                                    ins_to_pos = sa['end']
                                else:
                                    ins_from_pos = sa['start']
                            else:
                                if ins_read_part == 1:
                                    ins_from_pos = sa['start']
                                else:
                                    ins_to_pos = sa['end']

                            # if read_strand == sa['strand']:
                            #    sv_type = 'DELINS'
                            # else:
                            #    sv_type = 'DELINS2'

                            if query_pos1['query_start'] < query_pos2['query_start']:
                                if abs(query_pos1['query_end'] - query_pos2['query_start']) > config['split_gap_threshold']:
                                    continue
                            else:
                                if abs(query_pos2['query_end'] - query_pos1['query_start']) > config['split_gap_threshold']:
                                    continue

                            sv_type = ''
                            supp_type = 'adhoc'
                            if del_from_pos:
                                if ins_from_pos:
                                    sv_type = 'DEL-FROM-INS-FROM'
                                else:
                                    sv_type = 'DEL-FROM-INS-TO'
                            else:
                                if ins_from_pos:
                                    sv_type = 'DEL-TO-INS-FROM'
                                else:
                                    sv_type = 'DEL-TO-INS-TO'

                            if query_pos1['query_start'] < query_pos2['query_start']:
                                bp_query_pos = query_pos1['query_end']
                            else:
                                bp_query_pos = query_pos2['query_end']

                            bp_start = del_from_pos if del_from_pos else del_to_pos
                            bp_end = ins_from_pos if ins_from_pos else ins_to_pos

                            if del_chrom == ins_chrom and abs(bp_start - bp_end) < 1000000:
                                continue

                            bp_read_strand = '+'
                            if 'DEL-FROM' in sv_type:
                                bp_read_strand = '+' if read_strand == '+' else '-'
                            else:
                                bp_read_strand = '+' if read_strand == '-' else '-'

                            if read_name not in indel_dict:
                                indel_dict[read_name] = {
                                    'bp_chrom': del_chrom,
                                    'bp_start': bp_start,
                                    'bp_chrom2': ins_chrom,
                                    'bp_end': bp_end,
                                    'sv_type': sv_type,
                                    'read_name': read_name,
                                    'read_strand': bp_read_strand,
                                    'bp_query_pos': bp_query_pos,
                                    'read_seq': '',
                                    'supp_type': supp_type,
                                    'overlap': 0
                                }

        return {'sv_region_list': indel_dict.values()}

    def run(self):
        self.load_chrom_info()
        self.load_depth_region()

        chrom_list = self.options.chrom_list.split(',')

        results = []
        pool = Pool(len(chrom_list))
        results = pool.map(self, chrom_list)
        pool.close()
        pool.join()

        indel_result = self.call_indel()
        results.append(indel_result)

        if self.options.out_bed2:
            self.gen_bed2(results)


def get_split_read_caller_options(options):
    bam_file = options.bam_file
    chrom_list = options.chrom_list
    min_sv_size = options.min_sv_size
    max_sv_size = options.max_sv_size
    target_sv_type = options.target_sv_type
    out_bed2 = options.out_bed2
    depth_bed2 = options.depth_bed2
    term_threshold = options.term_threshold
    lt_10k_bed2 = options.lt_10k_bed2
    fai_file = options.fai_file

    return SplitReadCallerOptions(
        bam_file,
        chrom_list,
        min_sv_size,
        max_sv_size,
        target_sv_type,
        out_bed2,
        depth_bed2,
        term_threshold,
        lt_10k_bed2,
        fai_file,
    )


if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-bam', '--bam_file', help='bam file', required=True)
    parser.add_argument('-chromlist', '--chrom_list', help='chrom list', required=True)
    parser.add_argument('-minsize', '--min_sv_size', help='min sv size', required=True, type=int)
    parser.add_argument('-maxsize', '--max_sv_size', help='max sv size', required=True, type=int)
    parser.add_argument('-svtype', '--target_sv_type', help='target sv type', required=True)
    parser.add_argument('-outbed2', '--out_bed2', help='output bed2', required=True)

    parser.add_argument('-depth_bed2', '--depth_bed2', help='depth bed', required=False)
    parser.add_argument('-term_threshold', '--term_threshold', help='term threshold', required=False, type=int)
    parser.add_argument('-lt_10k_bed2', '--lt_10k_bed2', help='less-than-10k bed2', required=True)
    parser.add_argument('-fai_file', '--fai_file', help='ref. index file', required=True)

    options = parser.parse_args()

    split_read_caller_options = get_split_read_caller_options(options)

    split_read_caller = SplitReadCaller(split_read_caller_options)
    split_read_caller.run()
