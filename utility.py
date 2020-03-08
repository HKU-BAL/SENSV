from os import path
from multiprocessing import Pool
from collections import defaultdict

# ASSEMBLE_SCRIPT = '%s/assemble.sh' % path.dirname(path.realpath(__file__))

config = None

cytobands = {}


def init_logger():
    import logging
    from sys import stdout

    formatter = logging.Formatter('[%(asctime)s - %(levelname)s] - %(message)s')
    handler = logging.StreamHandler(stdout)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(handler)


def log(message):
    import logging
    from time import time
    from functools import wraps

    def inner_function(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            logging.info(f'{message}: begin')
            start_time = time()

            function(*args, **kwargs)

            end_time = time()
            logging.info(f'{message}: end. time_elapsed: {int(end_time - start_time)} sec')

        return wrapper

    return inner_function


def load_common_config():
    return {
        'cType': get_var('denovo_common', 'cType'),
        'cTypeConsumeRef': get_var('denovo_common', 'cTypeConsumeRef'),
        'cTypeConsumeQuery': get_var('denovo_common', 'cTypeConsumeQuery'),
        'cTypeConsumeRead': get_var('denovo_common', 'cTypeConsumeRead'),

        'minimap2': get_var('common', 'minimap2'),
        'samtools': get_var('common', 'samtools'),
        'ref_ver': get_var('common', 'ref_ver'),
        'ref': get_var('common', 'ref_%s' % (get_var('common', 'ref_ver'))),
    }


def optional_string_to_optional_int(string):
    if string is None:
        return None
    return int(string)


def get_var(group, var):
    global config
    if not config:
        load_config()

    try:
        return config[group][var]
    except:
        return None


def set_var(group, var, value):
    global config
    if not config:
        load_config()

    config[group][var] = value


def exit_on_not_found(file_path, message=None):
    from pathlib import Path
    if Path(file_path).is_file():
        return

    if message is not None:
        import logging
        logging.error(f'{message}. exit.')
    exit(0)


def base_directory():
    from pathlib import Path
    return Path(path.dirname(__file__))


def load_config():
    from configparser import ConfigParser
    global config

    config_file = base_directory() / "config.ini"

    config = ConfigParser()
    config.read(config_file)


def get_ref(chr, start_orig, end):
    import pysam

    ref_ver = get_var('common', 'ref_ver')

    padding_N = 'N' * (-start_orig)
    start = max(start_orig, 1) - 1

    ref = pysam.FastaFile(get_var('common', 'ref_%s' % ref_ver))

    return (padding_N + ref.fetch(str(chr), start, end)).upper()


def rev_comp(seq):
    trans = str.maketrans('ACGT', 'TGCA')
    return seq.translate(trans)[::-1]


def get_seq_from_fastq(read_name, fastq_prefix, ext='.fastq.gz', return_all=False):
    import subprocess
    from sys import exit
    from glob import glob

    GET_SEQ_FROM_FASTQ_SCRIPT = '%s/get_seq_by_name.sh' % path.dirname(path.realpath(__file__))

    seq = ''
    for filename in glob(fastq_prefix + ext):
        cmd = f'{GET_SEQ_FROM_FASTQ_SCRIPT} {read_name} {filename}'
        output = subprocess.check_output(cmd, shell=True, universal_newlines=True).splitlines()

        if not (len(output) == 4 and output[0].startswith('@' + read_name)):
            continue

        if return_all:
            return output
        else:
            seq = output[1]
            break

    if not seq:
        print("ERROR: get_seq_from_fastq: read %s not found. fastq_prefix=%s" % (read_name, fastq_prefix))
        exit(0)

    return seq


def cigar_string_to_tuples(cigar_string):
    from re import findall
    cigar_tuples = findall(r'(\d+)(\w)', cigar_string)

    return [(cigar_tuple[1], int(cigar_tuple[0])) for cigar_tuple in cigar_tuples]


def get_compact_cigar_string(cigar_string, start_clip_count, end_clip_count):
    cond_cigar = {'startS': 0, 'M': 0, 'D': 0, 'I': 0, 'endS': 0}

    cigar_tuples = cigar_string_to_tuples(cigar_string)
    first_cigar_tuple, last_cigar_tuple = cigar_tuples[0], cigar_tuples[-1]

    if first_cigar_tuple[0] in 'HS':
        cond_cigar['startS'] += first_cigar_tuple[1] + start_clip_count

    if last_cigar_tuple[0] in 'HS':
        cond_cigar['endS'] += last_cigar_tuple[1] + end_clip_count

    for (operation_type, length) in cigar_tuples:
        if operation_type in 'MID':
            cond_cigar[operation_type] += length

    cond_cigar_string = ''
    for key in ['startS', 'M', 'I', 'D', 'endS']:
        if cond_cigar[key] < 0:
            continue

        if key in ['startS', 'endS']:
            operation_type = 'S'
        else:
            operation_type = key

        cond_cigar_string += '%d%s' % (cond_cigar[key], operation_type)

    return cond_cigar_string


def run_shell_cmd(cmd, stdin_input=None):
    from subprocess import Popen, PIPE

    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE, universal_newlines=True)
    out, _err = p.communicate(stdin_input) if stdin_input else p.communicate()

    return out


def load_cytobands():
    global cytobands

    cytoband_file_path = base_directory() / 'data_file' / 'cytoband'

    cytobands = defaultdict(list)
    with open(cytoband_file_path, 'r') as f:
        for line in f:
            arr = line.strip().split()
            if arr[0][0] == '#':
                continue

            chr = arr[0].replace('chr', '')
            if not chr.isdigit() and chr not in ['X', 'Y']:
                continue

            chr_start, chr_end = int(arr[1]), int(arr[2])
            name = arr[3]
            #gieStain = arr[4]

            cytobands[chr].append({'chromStart': chr_start, 'chromEnd': chr_end, 'name': name})


def get_cytoband(chr, pos_str):
    global cytobands

    # not thread-safe ?
    # if not cytobands:
    #    load_cytobands()

    for entry in cytobands[chr]:
        if entry['chromStart'] <= int(pos_str) < entry['chromEnd']:
            return entry['name']

    return ''


def is_same_arm(chrom1, pos1, chrom2, pos2):
    cytoband1 = get_cytoband(chrom1, pos1)
    cytoband2 = get_cytoband(chrom2, pos2)

    return (cytoband1 and cytoband2 and cytoband1[0] == cytoband2[0])


def align(ref_info_list, query_info_list, file_prefix):
    # gen fasta
    with open(file_prefix + '.fasta', 'w') as fasta:
        for query_info in query_info_list:
            query_seq_name = query_info['query_seq_name']
            query_seq = query_info['query_seq']

            print('>%s' % query_seq_name, file=fasta)
            print('%s' % query_seq, file=fasta)

    # gen fa
    with open(file_prefix + '.fa', 'w') as ref:
        for ref_info in ref_info_list:
            ref_seq_name = ref_info['ref_seq_name']
            ref_seq = ref_info['ref_seq']

            print('>%s' % ref_seq_name, file=ref)
            sequence_list = [ref_seq[i: i + 70] for i in range(0, len(ref_seq), 70)]
            for sequence in sequence_list:
                print('%s' % sequence, file=ref)

    # realign
    # cmd = 'timeout 120 '
    cmd = ''
    if get_var('common', 'aligner') == 'minimap2':
        minimap2 = get_var('common', 'minimap2')
        cmd += "%s -t 48 -a %s.fa %s.fasta" % (minimap2, file_prefix, file_prefix)
    else:
        ngmlr = get_var('common', 'ngmlr')
        cmd += "%s -t 24 -r %s.fa -q %s.fasta -x ont" % (ngmlr, file_prefix, file_prefix)

    samtools = get_var('common', 'samtools')
    cmd += " | %s sort -@ 24 -o %s.bam - && %s index %s.bam &> %s.log" % (
        samtools, file_prefix, samtools, file_prefix, file_prefix)

    run_shell_cmd(cmd)


def gen_altref_align_seq(sv_str, working_path, fastq_prefix, query_read_name_list, config):
    from os import makedirs

    arr = sv_str.split('_')
    bp_chrom, bp_start, bp_chrom2, bp_end, sv_type = arr[0], int(arr[1]), arr[2], int(arr[3]), arr[4]

    # ref pos
    ref_start = max(bp_start - config['ref_buf_size'], 1)
    ref_end = bp_end + config['ref_buf_size']

    # file_prefix
    working_path = '%s/%s_%d_%s_%d_%s/validate' % (working_path, bp_chrom, bp_start, bp_chrom2, bp_end, sv_type)
    try:
        makedirs(working_path)
    except:
        pass
    file_prefix = '%s/%s_%d_%s_%d_%s' % (working_path, bp_chrom, bp_start, bp_chrom2, bp_end, sv_type)

    ref_seq = ''
    bp_seq = ''
    if sv_type == 'DEL':
        bp_seq = ''
    elif sv_type == 'INV':
        bp_seq += rev_comp(get_ref(bp_chrom, bp_start, bp_end))
    elif sv_type == 'DUP':
        bp_seq += get_ref(bp_chrom, bp_start, bp_end)
        bp_seq += get_ref(bp_chrom, bp_start, bp_end)
    elif sv_type == 'INS':
        bp_seq = 'XXX'
    elif sv_type in ['TRA', 'TRA_INV']:
        bp_seq = ''

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
        altref_seq += bp_seq
        altref_seq += get_ref(bp_chrom2, bp_end + 1, ref_end)

    altref_seq_name_list = ['altref']
    ref_info_list = [{'ref_seq_name': 'altref', 'ref_seq': altref_seq}]

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

    return {
        'file_prefix': file_prefix,
        'ref_seq_name_list': ref_seq_name_list,
        'altref_seq_name_list': altref_seq_name_list,
        'ref_info_list': ref_info_list,
        'query_info_list': query_info_list,
    }


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
            clip['end']['seq'] = get_hard_clip_seq(read, start_hard_clip_count,
                                                   query_pos, clip['end']['length'], fastq_prefix)
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


def load_depth_list():
    depth_list_path = base_directory() / 'depth' / 'depth_list'

    l = defaultdict(list)
    with open(depth_list_path, 'r') as f:
        for line in f:
            chrom, start, end, depth = line.strip().split("\t")
            start = int(start)
            end = int(end)
            depth = float(depth)

            l[chrom].append((start, end, depth))

    return l


def get_normal_avg_depth(sv_str):
    arr = sv_str.split('_')
    chrom, start, end = arr[0], int(arr[1]), int(arr[2])

    depth_sum = 0
    depth_count = 0

    depth_list = load_depth_list()
    if chrom in depth_list:
        for depth in depth_list[chrom]:
            if depth[0] >= start and depth[1] <= end:
                depth_sum += depth[2]
                depth_count += 1

    if depth_count != 0:
        return {
            'sv_str': sv_str,
            'normal_avg_depth': 1.0*depth_sum/depth_count
        }
    else:
        return 0


def filter_depth_file(gender, orig_depth_file, new_depth_file, nprocs):
    # load depth file
    depth_regions = []

    sv_str_list = []
    with open(orig_depth_file) as f:
        for line in f:
            if not line.strip() or line[0] == '#':
                continue

            arr = line.strip().split()
            arr[0] = 'X' if arr[0] == '23' else 'Y' if arr[0] == '24' else arr[0]

            chrom, start, end, score, region_type = \
                arr[0], int(float(arr[1])), int(float(arr[2])), int(float(arr[3])), arr[4]

            sv_str = '%s_%d_%d_%d_%s' % (chrom, start, end, score, region_type)
            sv_str_list.append(sv_str)

    pool = Pool(processes=nprocs)
    results = pool.map(get_normal_avg_depth, sv_str_list)
    pool.close()
    pool.join()

    for result in results:
        sv_str = result['sv_str']
        normal_avg_depth = result['normal_avg_depth']

        is_skip = 0
        if result['normal_avg_depth'] < 0.5:
            is_skip = 1
        elif gender == 'f' or chrom not in ['X', 'Y']:
            if gender == 'f' and chrom == 'Y':
                is_skip = 1

        if is_skip:
            print('gender = %s, normal avg depth = %f, depth region skipped: %s_%s_%s' %
                  (gender, normal_avg_depth, chrom, start, end))
            continue

        arr = sv_str.split('_')
        chrom, start, end, score, region_type = arr[0], arr[1], arr[2], arr[3], arr[4]
        depth_regions.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'score': score,
            'region_type': region_type,
        })

    # save depth file
    import csv
    with open(new_depth_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for region in depth_regions:
            row = [
                region['chrom'],
                str(region['start']),
                str(region['end']),
                str(region['score']),
                region['region_type'],
            ]
            writer.writerow(row)


def get_sorted_sv_str_list(sv_str_dict):
    def toIntStr(text): return text if text.isdigit() else '23' if text == 'X' else '24'
    return sorted(
        sv_str_dict,
        key=lambda k: (
            k.split('_')[4].rjust(17) +
            toIntStr(k.split('_')[0]).rjust(2) +
            k.split('_')[1].rjust(9) +
            toIntStr(k.split('_')[2]).rjust(2) +
            k.split('_')[3].rjust(9)
        )
    )


def get_short_name(read_name):
    return read_name[0:12]
