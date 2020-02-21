from __future__ import print_function
import pysam
import sys
import os
import re
from random import randint
from time import gmtime, strftime
import subprocess
from shutil import copyfile
import glob
import string
import configparser
import csv
from multiprocessing import Pool

ASSEMBLE_SCRIPT = '%s/assemble.sh' % os.path.dirname(os.path.realpath(__file__))
GET_SEQ_FROM_FASTQ_SCRIPT = '%s/get_seq_by_name.sh' % os.path.dirname(os.path.realpath(__file__))

config = None

cytobands = {}

# working: start

cType = 'MIDNSHP=X'
cTypeConsumeRef = 'MDN=X'
cTypeConsumeQuery = 'MIS=X'
cTypeConsumeRead = 'MIS=XH'

def has_hard_clip(read):
    tuples = read.cigartuples

    if cType[tuples[0][0]] in 'H' or cType[tuples[-1][0]] in 'H':
        return True

    return False

def get_seq(read, fastq_prefix, ext='.fastq.gz', return_all=False):
    result = has_hard_clip(read)

    if has_hard_clip(read):
        seq = get_seq_from_fastq(read.query_name, fastq_prefix, ext, return_all)
    else:
        if read.is_reverse:
            seq = rev_comp(read.query_sequence)
        else:
            seq = read.query_sequence

    return seq

# working: end

def load_common_config():
    config = {}

    config['cType'] = get_var('denovo_common', 'cType')
    config['cTypeConsumeRef'] = get_var('denovo_common', 'cTypeConsumeRef')
    config['cTypeConsumeQuery'] = get_var('denovo_common', 'cTypeConsumeQuery')
    config['cTypeConsumeRead'] = get_var('denovo_common', 'cTypeConsumeRead')

    config['minimap2'] = get_var('common', 'minimap2')
    config['samtools'] = get_var('common', 'samtools')
    config['ref_ver'] = get_var('common', 'ref_ver')
    config['ref'] = get_var('common', 'ref_%s' % (config['ref_ver']))

    return config


def get_var(group, var):
    global config

    if not config:
        load_config()

    return config[group][var]


def load_config():
    global config

    config = configparser.ConfigParser()
    config.read('config.ini')


def get_ref(chr, start_orig, end):
    ref_ver = get_var('common', 'ref_ver')

    ref = pysam.FastaFile(get_var('common', 'ref_%s' % ref_ver))

    padding = ''
    if start_orig < 1:
        padding = 'N' * (-start_orig)
        start = 1
    else:
        start = start_orig

    nts = padding + ref.fetch(str(chr), start-1, end)

    return nts.upper()


def rev_comp(seq):
    trans = string.maketrans('ACGT', 'TGCA')
    return seq.translate(trans)[::-1]


def get_read_mapq(read):
    mapq = read.mapping_quality
    return mapq


def get_read_qual(read):
    mapq = get_read_mapq(read)
    qual = 1 - 10 ** (-mapq / 10.0)
    return qual


def get_seq_from_fastq(read_name, fastq_prefix, ext='.fastq.gz', return_all=False):
    seq = ''

    files = glob.glob(fastq_prefix + ext)
    isFound = 0
    for filename in files:
        cmd = '%s %s %s' % (GET_SEQ_FROM_FASTQ_SCRIPT, read_name, filename)
        output = subprocess.check_output(cmd, shell=True).splitlines()

        if len(output) == 4 and output[0].startswith('@' + read_name):
            if return_all:
                return output
            else:
                seq = output[1]
                isFound = 1
                break

    if not seq:
        print("ERROR: get_seq_from_fastq: read %s not found. fastq_prefix=%s" % (read_name, fastq_prefix))
        sys.exit(0)

    return seq


def cigar_string_to_tuples(cigar_string):
    tuples = re.findall(r'(\d+)(\w)', cigar_string)

    cigar_tuples = []
    for tuple in tuples:
        cigar_tuples.append([tuple[1], int(tuple[0])])

    return cigar_tuples

def get_compact_cigar_string(cigar_string, start_clip_count, end_clip_count):
    cond_cigar = {'startS': 0, 'M': 0, 'D': 0, 'I': 0, 'endS': 0}

    cigar_tuples = cigar_string_to_tuples(cigar_string)

    if cigar_tuples[0][0] in 'HS':
        cond_cigar['startS'] += int(cigar_tuples[0][1]) + start_clip_count

    if cigar_tuples[-1][0] in 'HS':
        cond_cigar['endS'] += int(cigar_tuples[-1][1]) + end_clip_count

    for (type, length) in cigar_tuples:
        if type == 'M':
            cond_cigar['M'] += int(length)
        elif type == 'I':
            cond_cigar['I'] += int(length)
        elif type == 'D':
            cond_cigar['D'] += int(length)

    cond_cigar_string = ''
    for key in ['startS', 'M', 'I', 'D', 'endS']:
        if cond_cigar[key] >= 0:
            if key in ['startS', 'endS']:
                type = 'S'
            else:
                type = key
            cond_cigar_string += '%d%s' % (cond_cigar[key], type)

    return cond_cigar_string


def run_shell_cmd(cmd, stdinInput=None):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    if stdinInput:
        out, err = p.communicate(stdinInput)
    else:
        out, err = p.communicate()

    #if err:
    #    print(err)

    return out


def load_cytobands():
    global cytobands

    cytobands = {}
    with open('data_file/cytoband') as f:
        for line in f:
            arr = line.strip().split()
            if arr[0][0] == '#':
                continue

            chrom = arr[0].replace('chr', '')
            if not chrom.isdigit() and chrom not in ['X', 'Y']:
                continue

            chromStart = int(arr[1])
            chromEnd = int(arr[2])
            name = arr[3]
            #gieStain = arr[4]

            if chrom not in cytobands:
                cytobands[chrom] = []

            cytobands[chrom].append({'chromStart': chromStart, 'chromEnd': chromEnd, 'name': name})


def get_cytoband(chrom, pos_str):
    global cytobands

    ##not thread-safe ?
    #if not cytobands:
    #    load_cytobands()

    cytoband = ''
    pos = int(pos_str)

    list = cytobands[chrom]
    found = 0
    for entry in list:
        if pos >= entry['chromStart'] and pos < entry['chromEnd']:
            cytoband = entry['name']
            break

    return cytoband


def is_same_arm(chrom1, pos1, chrom2, pos2):
    cytoband1 = get_cytoband(chrom1, pos1)
    cytoband2 = get_cytoband(chrom2, pos2)

    if cytoband1 and cytoband2 and cytoband1[0] == cytoband2[0]:
        return True
    else:
        return False


def align(refInfoList, queryInfoList, file_prefix):
    # gen fasta
    fasta_name = file_prefix + '.fasta'
    fasta = open(fasta_name, 'w')

    for queryInfo in queryInfoList:
        query_seq_name = queryInfo['query_seq_name']
        query_seq = queryInfo['query_seq']

        print('>%s' % query_seq_name, file=fasta)
        print('%s' % query_seq, file=fasta)

    fasta.close()

    # gen fa
    ref_name = file_prefix + '.fa'
    ref = open(ref_name, 'w')

    for refInfo in refInfoList:
        ref_seq_name = refInfo['ref_seq_name']
        ref_seq = refInfo['ref_seq']

        print('>%s' % ref_seq_name, file=ref)
        seqArr = [ref_seq[i: i + 70] for i in range(0, len(ref_seq), 70)]
        for seq in seqArr:
            print('%s' % seq, file=ref)

    ref.close()

    # realign
    samtools = get_var('common', 'samtools')
    cmd = ''
    #cmd = 'timeout 120 '
    if get_var('common', 'aligner') == 'minimap2':
        cmd += "%s -t 48 -a %s.fa %s.fasta | %s sort -@ 24 -o %s.bam - && %s index -@ 24 %s.bam &> %s.log" % \
            (get_var('common', 'minimap2'), file_prefix, file_prefix, samtools, file_prefix, samtools, file_prefix, file_prefix)
    else:
        cmd += "%s -t 24 -r %s.fa -q %s.fasta -x ont | %s sort -@ 24 -o %s.bam - && %s index %s.bam &> %s.log" % \
            (get_var('common', 'ngmlr'), file_prefix, file_prefix, samtools, file_prefix, samtools, file_prefix, file_prefix)

    run_shell_cmd(cmd)


def gen_altref_align_seq(sv_str, working_path, fastq_prefix, query_read_name_list, config):
    arr = sv_str.split('_')
    bp_chrom, bp_start, bp_chrom2, bp_end, sv_type = arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4]

    # ref pos
    ref_start = bp_start - config['ref_buf_size']
    if ref_start < 1:
        ref_start = 1
    ref_end = bp_end + config['ref_buf_size']

    # file_prefix
    working_path = '%s/%s_%d_%s_%d_%s/validate' % (working_path, bp_chrom, bp_start, bp_chrom2, bp_end, sv_type)
    try:
        os.makedirs(working_path)
    except:
        pass
    file_prefix = '%s/%s_%d_%s_%d_%s' % (working_path, bp_chrom, bp_start, bp_chrom2, bp_end, sv_type)

    ref_seq = ''
    bpSeq = ''
    if sv_type == 'DEL':
        bpSeq = ''
    elif sv_type == 'INV':
        bpSeq += rev_comp(get_ref(bp_chrom, bp_start, bp_end))
    elif sv_type == 'DUP':
        bpSeq += get_ref(bp_chrom, bp_start, bp_end)
        bpSeq += get_ref(bp_chrom, bp_start, bp_end)
    elif sv_type == 'INS':
        bpSeq = 'XXX'
    elif sv_type in ['TRA', 'TRA_INV']:
        bpSeq = ''

    ref_info_list = []
    # altref seq
    if sv_type == 'TRA':
        altref_seq = get_ref(bp_chrom, bp_start - config['ref_buf_size'], bp_start - 1)
        altref_seq += get_ref(bp_chrom2, bp_end + 1, bp_end + config['ref_buf_size'])

        altref_seq += get_ref(bp_chrom2, bp_end - config['ref_buf_size'], bp_end)
        altref_seq += get_ref(bp_chrom, bp_start, bp_start + config['ref_buf_size'])
    elif sv_type == 'TRA_INV':
        altref_seq = get_ref(bp_chrom, bp_start - config['ref_buf_size'], bp_start - 1)
        altref_seq += rev_comp(get_ref(bp_chrom2, bp_end - config['ref_buf_size'], bp_end))

        altref_seq += rev_comp(get_ref(bp_chrom, bp_start, bp_start + config['ref_buf_size']))
        altref_seq += get_ref(bp_chrom2, bp_end + 1, bp_end + config['ref_buf_size'])
    else:
        altref_seq = get_ref(bp_chrom, ref_start, bp_start - 1)
        altref_seq += bpSeq
        altref_seq += get_ref(bp_chrom2, bp_end + 1, ref_end)

    altref_seq_name_list = ['altref']
    ref_info_list.append({'ref_seq_name': 'altref', 'ref_seq': altref_seq})

    # ref seq
    if sv_type in ['TRA', 'TRA_INV']:
        ref_seq = get_ref(bp_chrom, bp_start - config['ref_buf_size'], bp_start + config['ref_buf_size'])
        ref_seq += get_ref(bp_chrom2, bp_end - config['ref_buf_size'], bp_end + config['ref_buf_size'])
    else:
        #ref_seq = get_ref(bp_chrom, ref_start, ref_end)
        max_chrom_length = int(get_var('sv_validator', 'max_chrom_length'))
        ref_seq = get_ref(bp_chrom, 1, max_chrom_length)

    ref_seq_name_list = ['ref']
    ref_info_list.append({'ref_seq_name': 'ref', 'ref_seq': ref_seq})

    query_info_list = []
    for query_seq_name in query_read_name_list:
        query_seq = get_seq_from_fastq(query_seq_name, fastq_prefix)
        query_info_list.append({'query_seq_name': query_seq_name, 'query_seq': query_seq})

    return {'file_prefix': file_prefix, 'ref_seq_name_list': ref_seq_name_list, 'altref_seq_name_list': altref_seq_name_list,
            'ref_info_list': ref_info_list, 'query_info_list': query_info_list}

def get_clip_seq(read, fastq_prefix):
    config = load_common_config()

    clip = {'start': {'length': 0, 'seq': ''}, 'end': {'length': 0, 'seq': ''}}

    cigar_tuples = read.cigartuples
    reference_start = read.reference_start+1
    query_seq = read.query_sequence

    # start clip
    tuple = cigar_tuples[0]
    clip['start']['length'] = tuple[1]
    if config['cType'][tuple[0]] in 'HS' and clip['start']['length'] > int(get_var('sv_realigner', 'split_read_clip_threshold')):
        if config['cType'][tuple[0]] == 'S':
            clip['start']['seq'] = query_seq[:tuple[1]]
            clip['start']['query_pos'] = 0
        else:
            clip['start']['seq'] = get_hard_clip_seq(read, 0, 0, clip['start']['length'], fastq_prefix)
            clip['start']['query_pos'] = 0

    # end clip
    tuple = cigar_tuples[len(cigar_tuples)-1]
    clip['end']['length'] = tuple[1]
    if config['cType'][tuple[0]] in 'HS' and clip['end']['length'] > int(get_var('sv_realigner', 'split_read_clip_threshold')):
        query_pos = 0
        ref_pos = reference_start

        for (type, length) in cigar_tuples[:-1]:
            if config['cType'][type] in config['cTypeConsumeRef']:
                ref_pos += length
            if config['cType'][type] in config['cTypeConsumeQuery']:
                query_pos += length

        if config['cType'][tuple[0]] == 'S':
            clip['end']['seq'] = query_seq[query_pos:]
            clip['end']['query_pos'] = query_pos
        else:
            start_hard_clip_count = 0
            if config['cType'][cigar_tuples[0][0]] == 'H':
                start_hard_clip_count = cigar_tuples[0][1]
            clip['end']['seq'] = get_hard_clip_seq(read, start_hard_clip_count, query_pos, clip['end']['length'], fastq_prefix)
            clip['end']['query_pos'] = start_hard_clip_count + query_pos

    return clip

def get_hard_clip_seq(read, start_hard_clip_count, query_start_pos, length, fastq_prefix):
    seq = get_seq_from_fastq(read.query_name, fastq_prefix)

    if seq:
        if read.is_reverse:
            seq = rev_comp(seq)

        if query_start_pos == 0:
            if length == -1:
                return seq
            else:
                return seq[:length]
        else:
            start = start_hard_clip_count + query_start_pos
            if length == -1:
                return seq[start:]
            else:
                end = start + length
                return seq[start:end]
    else:
        print('cannot find %s' % read.query_name)
        return ''

def gen_altref_seq(sv_str, buf_size):
    bp_chrom, bp_start, bp_chrom2, bp_end, bp_type = sv_str.split('_')
    bp_start, bp_end = int(bp_start), int(bp_end)

    ref_info = []

    if bp_type not in ['TRA', 'TRAINV']:
        desc = ''
        seq = get_ref(bp_chrom, bp_start-1-buf_size, bp_start-1)

        if bp_type == 'DEL':
            pass
        elif bp_type == 'INV':
            seq += rev_comp(get_ref(bp_chrom, bp_start, bp_end))
        elif bp_type == 'DUP':
            seq += get_ref(bp_chrom, bp_start, bp_end)
            seq += get_ref(bp_chrom, bp_start, bp_end)
        elif bp_type == 'INS':
            seq = 'N' * (bp_end - bp_start)

        seq += get_ref(bp_chrom, bp_end+1, bp_end+1+buf_size)
        ref_info.append({'ref_seq_name': 'altref_%s' % sv_str, 'ref_seq': seq, 'desc': desc})
    else:
        desc = ''
        desc2 = ''
        if bp_type == 'TRA':
            seq = get_ref(bp_chrom, 1, bp_start-1)
            desc = 'bp:%s' % len(seq)
            seq += get_ref(bp_chrom2, bp_end, config['chrom_max_size'])
            seq2 = get_ref(bp_chrom2, 1, bp_end)
            desc2 = 'bp:%s' % len(seq2)
            seq2 += get_ref(bp_chrom, bp_start+1, config['chrom_max_size'])
        elif bp_type == 'TRAINV':
            seq = get_ref(bp_chrom, 1, bp_start-1)
            desc = 'bp:%s' % len(seq)
            seq += rev_comp(get_ref(bp_chrom2, 1, bp_end-1))
            seq2 = rev_comp(get_ref(bp_chrom, bp_start, config['chrom_max_size']))
            desc2 = 'bp:%s' % len(seq2)
            seq2 += get_ref(bp_chrom2, bp_end+1, config['chrom_max_size'])

        ref_info.append({'ref_seq_name': 'altref_%s' % sv_str, 'ref_seq': seq, 'desc': desc})
        ref_info.append({'ref_seq_name': 'altref_%s_2' % sv_str, 'ref_seq': seq2, 'desc': desc2})

    return ref_info

"""
def get_normal_avg_depth(sv_str):
    bam_file = get_var('chain_caller', 'normal_male_bam')

    arr = sv_str.split('_')
    chrom, start, end = arr[0], int(arr[1]), int(arr[2])

    max_depth = 8
    cmd = "samtools depth -a -a %s -r %s:%d-%d | \
          awk '{ max=%d; depth = ($3<max?$3:max); sum += depth } END { if (NR > 0) print sum / NR }'" % \
          (bam_file, chrom, start, end, max_depth)

    output = run_shell_cmd(cmd).split('\n')

    return float(output[0])

def filter_depth_file(gender, orig_depth_file, new_depth_file):
    # load depth file
    depth_regions = []
    with open(orig_depth_file) as f:
        for line in f:
            if not line.strip() or line[0] == '#':
                continue

            arr = line.strip().split()
            arr[0] = 'X' if arr[0] == '23' else 'Y' if arr[0] == '24' else arr[0]

            chrom, start, end, score, region_type = arr[0], int(float(arr[1])), int(float(arr[2])), int(float(arr[3])), arr[4]

            sv_str = '%s_%d_%d' % (chrom, start, end)
            normal_avg_depth = get_normal_avg_depth(sv_str)

            is_skip = 0
            if normal_avg_depth < 0.5:
                is_skip = 1
            elif gender == 'f' or chrom not in ['X','Y']:
                if gender == 'f' and chrom == 'Y':
                    is_skip = 1

            if is_skip:
                print('gender = %s, normal avg depth = %f, depth region skipped: %s_%s_%s' % (gender, normal_avg_depth, chrom, start, end))
                continue

            depth_regions.append({'chrom': chrom, 'start': start, 'end': end, 'score': score, 'region_type': region_type})

    # save depth file
    with open(new_depth_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for region in depth_regions:
            row = [region['chrom'], str(region['start']), str(region['end']), str(region['score']), region['region_type']]
            writer.writerow(row)
"""

"""
def get_normal_avg_depth(sv_str):

    bam_file = get_var('chain_caller', 'normal_male_bam')

    arr = sv_str.split('_')
    chrom, start, end = arr[0], int(arr[1]), int(arr[2])
    
    max_depth = 8
    cmd = "samtools depth -a -a %s -r %s:%d-%d | \
          awk '{ max=%d; depth = ($3<max?$3:max); sum += depth } END { if (NR > 0) print sum / NR }'" % \
          (bam_file, chrom, start, end, max_depth)

    output = run_shell_cmd(cmd).split('\n')

    return {'sv_str':sv_str, 'normal_avg_depth':float(output[0])}
"""
def load_depth_list():

    f = open("depth/depth_list", 'r')
    l = {}

    for line in f:
        chrom, start, end, depth = line.strip().split("\t")
        start = int(start)
        end = int(end)
        depth = float(depth)

        if chrom in l:
            l[chrom].append((start, end, depth))
        else:
            l[chrom] = []
            l[chrom].append((start, end, depth))

    return l

def get_normal_avg_depth(sv_str):

    depth_list = load_depth_list()

    arr = sv_str.split('_')
    chrom, start, end = arr[0], int(arr[1]), int(arr[2])

    depth_sum = 0
    depth_count = 0

    if chrom in depth_list:
        for depth in depth_list[chrom]:
            if depth[0]>=start and depth[1]<=end:
                depth_sum += min(8,depth[2])
                depth_count += 1

    if depth_count != 0:
        return {'sv_str':sv_str, 'normal_avg_depth':1.0*depth_sum/depth_count}
    else:
        return 0

        

def filter_depth_file(gender, orig_depth_file, new_depth_file):
    # load depth file
    depth_regions = []

    sv_str_list = []
    with open(orig_depth_file) as f:
        for line in f:
            if not line.strip() or line[0] == '#':
                continue

            arr = line.strip().split()
            arr[0] = 'X' if arr[0] == '23' else 'Y' if arr[0] == '24' else arr[0]

            chrom, start, end, score, region_type = arr[0], int(float(arr[1])), int(float(arr[2])), int(float(arr[3])), arr[4]

            sv_str = '%s_%d_%d_%d_%s' % (chrom, start, end, score, region_type)
            sv_str_list.append(sv_str)

    pool = Pool(100)
    results = pool.map(get_normal_avg_depth, sv_str_list)
    pool.close()
    pool.join()

    for result in results:
        is_skip = 0

        sv_str = result['sv_str']
        normal_avg_depth = result['normal_avg_depth']

        if result['normal_avg_depth'] < 0.5:
            is_skip = 1
        elif gender == 'f' or chrom not in ['X','Y']:
            if gender == 'f' and chrom == 'Y':
                is_skip = 1

        if is_skip:
            print('gender = %s, normal avg depth = %f, depth region skipped: %s_%s_%s' % (gender, normal_avg_depth, chrom, start, end))
            continue

        arr = sv_str.split('_')
        chrom, start, end, score, region_type = arr[0], arr[1], arr[2], arr[3], arr[4]
        depth_regions.append({'chrom': chrom, 'start': start, 'end': end, 'score': score, 'region_type': region_type})

    # save depth file
    with open(new_depth_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for region in depth_regions:
            row = [region['chrom'], str(region['start']), str(region['end']), str(region['score']), region['region_type']]
            writer.writerow(row)

def get_sorted_sv_str_list(sv_str_dict):
    toIntStr = lambda text: text if text.isdigit() else '23' if text == 'X' else '24'
    return sorted(sv_str_dict, key=lambda k: k.split('_')[4].rjust(17) + \
        toIntStr(k.split('_')[0]).rjust(2) + k.split('_')[1].rjust(9) + \
        toIntStr(k.split('_')[2]).rjust(2) + k.split('_')[3].rjust(9))

def get_short_name(read_name):
    return read_name[0:12]
