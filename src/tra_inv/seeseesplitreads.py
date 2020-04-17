import sys
import os
import shlex
import tempfile
from threading import Thread
from subprocess import PIPE
from pathlib import Path
from itertools import chain, repeat
from multiprocessing import Pool
from argparse import ArgumentParser
from collections import defaultdict
from intervaltree import IntervalTree

from sv import SV, SV_Utilities
from dp import DP, get_dp_result
from read import Read
from utils import get_contig_list, get_reference_from_fasta, interval_from, subprocess_popen, major_contigs_set
from samtoolsViewViewer import SamtoolsViewViewer


def reads_from_samtools(samtools, bam, ctgname, ctgrange, mapq):
    viewer = SamtoolsViewViewer(
        samtools=samtools,
        bam=bam,
        has_tag="SA",
        ctg_list=[ctgname] if ctgname else get_contig_list(bam, samtools),
        ctg_range=ctgrange if ctgrange else None,
        mapq_filter=mapq,
    )

    return viewer.reads


def breakdowned_reads_from(reads, processes):
    with Pool(processes=processes) as pool:
        breakdowned_reads = pool.map(Read.breakdowned_reads_from, reads)

    return chain.from_iterable(breakdowned_reads)


def filtered_reads_with(reads, mapq, min_sv_size, min_read_length):
    def get_SV_size(read):
        pos1, pos2 = read.split_position
        return abs(pos1 - pos2)

    filtered_reads = {}
    for read in reads:
        read2 = read.SA_reads[0]
        if read.RNAME not in major_contigs_set or read2.RNAME not in major_contigs_set:
            continue
        if int(read.MAPQ) < mapq or int(read2.MAPQ) < mapq:
            continue
        if read.SV_type not in ["INV", "TRA"]:
            continue
        if len(read.SEQ) < min_read_length:
            continue
        if read.CigarString.first_cigar_tuple[1] == "H" or read.CigarString.last_cigar_tuple[1] == "H":
            continue
        if read.SV_type != "TRA" and get_SV_size(read) < min_sv_size:
            continue
        if read.QNAME in filtered_reads and len(filtered_reads[read.QNAME].SEQ) >= len(read.SEQ):
            continue

        filtered_reads[read.QNAME + f'{read.split_no}'] = read

    return filtered_reads.values()


def filtered_reads_before_dp(reads):
    interval_trees = defaultdict(IntervalTree)

    for read in reads:
        contig1, contig2 = read.RNAME, read.SA_reads[0].RNAME
        breakpoint_position1, breakpoint_position2 = read.split_position
        if contig1 > contig2 or breakpoint_position1 > breakpoint_position2:
            contig1, contig2 = contig2, contig1
            breakpoint_position1, breakpoint_position2 = breakpoint_position2, breakpoint_position1

        interval1 = interval_from(breakpoint_position1, flanking_bases=1000)
        interval2 = interval_from(breakpoint_position2, flanking_bases=1000)

        # node_data is an interval trees storing ending intervals together with the supported read
        node_data = defaultdict(IntervalTree)
        node_data[read.SV_type + "#" + contig2].addi(interval2.begin, interval2.end, [read])

        interval_trees[contig1].addi(interval1.begin, interval1.end, node_data)

    # merge the overlapping starting intervals (node_data)
    def merge_interval_trees(current_reduced_interval_trees, new_interval_trees):
        for contig in new_interval_trees.keys():
            for interval_obj in new_interval_trees[contig]:
                current_reduced_interval_trees[contig].addi(*(interval_obj))
        return current_reduced_interval_trees
    for interval_tree_begin in interval_trees.values():
        interval_tree_begin.merge_overlaps(data_reducer=merge_interval_trees, strict=False)

    # merge the overlapping ending intervals
    def merge_node_data(current_reduced_node_data, new_node_data):
        return current_reduced_node_data + new_node_data
    for interval_tree_begin in interval_trees.values():
        for _, _, node_data in interval_tree_begin:
            for interval_tree_end in node_data.values():
                interval_tree_end.merge_overlaps(data_reducer=merge_node_data, strict=False)

    filtered_reads = []
    for interval_tree_begin in interval_trees.values():
        for _, _, node_data in interval_tree_begin:
            for interval_tree_end in node_data.values():
                for _, _, reads in interval_tree_end:
                    if len(reads) <= 1:
                        continue
                    filtered_reads += reads

    return filtered_reads


def dp_results_from(reads, ref, processes):
    print("... Performing dp on %d reads ..." % len(reads))
    with Pool(processes=processes) as pool:
        result = pool.starmap(get_dp_result, zip(reads, repeat(ref)))

    return result


def ngmlr_realign(ngmlr, ref, reads):
    """
    Realigns all output split reads with ngmlr, returns
    a dict with QNAME as key and pos as value.
    """

    def pump_input(pipe, reads):

        aligned_reads = {}

        with pipe:
            for read in reads:
                if read.QNAME in aligned_reads:
                    continue

                aligned_reads[read.QNAME] = True

                string = "@" + read.QNAME + "\n"
                string += read.SEQ + "\n"
                string += "+\n"
                string += "~"*len(read.SEQ) + "\n"

                pipe.write(string)

    def get_read_pos_info(line):
        read = Read.from_str(line)

        pos_array = [(read.RNAME, read.POS)]

        if "SA" in read.tags:
            for SA_read in read.SA_reads:
                pos_array.append((SA_read.RNAME, SA_read.POS))

        return read.QNAME, pos_array

    #aligned_reads = {}

    # with tempfile.NamedTemporaryFile() as temp_file:
    #    for read in reads:
    #        if read.QNAME in aligned_reads:
    #            continue

    #        aligned_reads[read.QNAME] = True

    #        string = "@" + read.QNAME + "\n"
    #        string += read.SEQ + "\n"
    #        string += "+\n"
    #        string += "~"*len(read.SEQ) + "\n"

    #        temp_file.write(string.encode("utf-8"))
    #        temp_file.flush()

    #        print(string)

    with subprocess_popen(shlex.split(f"{ngmlr} -t 48 -r {ref} -x ont"), stdin=PIPE, stdout=PIPE) as p:

        Thread(target=pump_input, args=[p.stdin, reads]).start()

        remapped_reads = {}

        with p.stdout:
            for line in iter(p.stdout.readline, b''):

                if len(line) == 0:
                    break

                if line[0] == "@":
                    continue

                QNAME, pos_array = get_read_pos_info(line)

                remapped_reads[QNAME] = pos_array

    return remapped_reads


def filtered_reads_with_remapping(reads, pos_dict):

    filtered_reads = []

    for read in reads:
        distance = 10001
        if read.QNAME not in pos_dict:
            continue

        for pos_info in pos_dict[read.QNAME]:
            if read.RNAME != pos_info[0]:
                continue

            if abs(int(read.POS)-int(pos_info[1])) < distance:
                distance = abs(int(read.POS)-int(pos_info[1]))

        if distance <= 10000:
            filtered_reads.append(read)

    return filtered_reads


def sv_candidates_from(reads, dp_results):
    interval_trees = defaultdict(IntervalTree)
    for (read, dp_result) in zip(reads, dp_results):
        output, similar_score = dp_result
        if similar_score < 1.0:
            continue
        # print(output, similar_score)

        sv_candidate = SV(
            read=read,
            breakpoint_position1=output[2],
            breakpoint_position2=output[3],
            dp_result=dp_result,
        )

        contig1, contig2 = sv_candidate.breakpoint_chromosomes
        interval1, interval2 = sv_candidate.breakpoint_intervals

        node_data = defaultdict(IntervalTree)
        node_data[sv_candidate.type + "#" + contig2].addi(interval2.begin, interval2.end, sv_candidate)

        interval_trees[contig1].addi(interval1.begin, interval1.end, node_data)

    # merge the overlapping starting intervals (node_data)
    def merge_interval_trees(current_reduced_interval_trees, new_interval_trees):
        for contig in new_interval_trees.keys():
            for interval_obj in new_interval_trees[contig]:
                current_reduced_interval_trees[contig].addi(*(interval_obj))
        return current_reduced_interval_trees
    for interval_tree_begin in interval_trees.values():
        interval_tree_begin.merge_overlaps(data_reducer=merge_interval_trees, strict=False)

    # merge the overlapping ending intervals
    def merge_sv_candidate(current_reduced_sv_candidate, new_sv_candidate):
        current_reduced_sv_candidate.merge(new_sv_candidate)
        return current_reduced_sv_candidate
    for interval_tree_begin in interval_trees.values():
        for _, _, node_data in interval_tree_begin:
            for interval_tree_end in node_data.values():
                interval_tree_end.merge_overlaps(data_reducer=merge_sv_candidate, strict=False)

    sv_candidates = []
    for interval_tree_begin in interval_trees.values():
        for _, _, node_data in interval_tree_begin:
            for interval_tree_end in node_data.values():
                for _, _, sv_candidate in interval_tree_end:
                    sv_candidates.append(sv_candidate)

    return sv_candidates


def is_valid_sv_with_args(sv, ref, allele_frequencies_from, min_af):
    if len(sv.supported_reads) <= 1:
        return False

    for chr, breakpoint_interval in zip(sv.breakpoint_chromosomes, sv.breakpoint_intervals):
        if "N" in get_reference_from_fasta(ref, chr, breakpoint_interval.begin, breakpoint_interval.end):
            return False

    for af in allele_frequencies_from(sv):
        if af < min_af:
            return False

    return True


def filtered_tra_sv_from(sv_candidates, min_sv_size):
    filtered_sv_candidate = []
    for sv_candidate in sv_candidates:
        if sv_candidate.type == "INV":

            start, end = sv_candidate.breakpoint1.position, sv_candidate.breakpoint2.position

            if abs(end-start) < min_sv_size:
                continue

            filtered_sv_candidate.append(sv_candidate)
            continue

        if sv_candidate.type == "TRA":
            read_exist = [False, False]
            for read in sv_candidate.supported_reads:
                if read.RNAME == sv_candidate.breakpoint1.chr:
                    pos1 = int(read.POS)
                    pos2 = int(read.POS) + read.CigarString.ref_length
                else:
                    pos1 = int(read.SA_reads[0].POS)
                    pos2 = int(read.SA_reads[0].POS) + read.SA_reads[0].CigarString.ref_length

                if abs(pos1-sv_candidate.breakpoint1.position) < abs(pos2-sv_candidate.breakpoint1.position):
                    read_exist[0] = True
                else:
                    read_exist[1] = True

            if read_exist[0] and read_exist[1]:
                filtered_sv_candidate.append(sv_candidate)

    return filtered_sv_candidate


def sv_list_from(sv_candidates, ref, allele_frequencies_from, min_af, processes):
    with Pool(processes=processes) as pool:
        is_valid_sv_list = pool.starmap(
            is_valid_sv_with_args,
            zip(sv_candidates, repeat(ref), repeat(allele_frequencies_from), repeat(min_af)),
        )

    for is_valid_sv, sv_candidate in zip(is_valid_sv_list, sv_candidates):
        if is_valid_sv:
            yield sv_candidate


def get_vcf_rows_from(sv_list, vcf_rows_from, processes):
    with Pool(processes=processes) as pool:
        vcf_rows = pool.map(vcf_rows_from, sv_list)

    return chain.from_iterable(vcf_rows)


def output_vcf(vcf_rows, output_file_path, header_file_path):
    with open(output_file_path, "w") as vcf:
        with open(header_file_path) as header:
            for line in header.read().splitlines():
                print(line, file=vcf)

        for vcf_row in vcf_rows:
            print(vcf_row, file=vcf)


def test(args):
    reads = reads_from_samtools(args.samtools, args.bam, args.ctgname, args.ctgrange, args.mapq)
    reads = breakdowned_reads_from(reads, args.processes)
    reads = filtered_reads_with(reads, args.mapq, args.min_sv_size, args.min_read_length)
    reads = filtered_reads_before_dp(reads)

    dp_results = dp_results_from(reads, ref=args.ref, processes=args.processes)
    #remap_pos_dict = ngmlr_realign("ngmlr", args.ref, reads)
    #reads = filtered_reads_with_remapping(reads, remap_pos_dict)

    sv_candidates = sv_candidates_from(reads, dp_results)
    sv_candidates = filtered_tra_sv_from(sv_candidates, args.min_sv_size)

    sv_utils = SV_Utilities(args.bam, args.samtools)
    sv_list = sv_list_from(sv_candidates, args.ref, sv_utils.allele_frequencies_from, args.min_af, args.processes)
    vcf_rows = get_vcf_rows_from(sv_list, vcf_rows_from=sv_utils.vcf_rows_from, processes=args.processes)
    output_vcf(vcf_rows, output_file_path=args.output, header_file_path=args.header_file_path)


if __name__ == "__main__":
    parser = ArgumentParser(description='A script to see see what split reads look like >.0')

    parser.add_argument('--bam', help='The bam file to be read :D', required=True)
    parser.add_argument('--ref', help='The corresponding fasta file >.0', required=True)
    parser.add_argument('--ctgname', help='contig name 0.0', required=False)
    parser.add_argument('--ctgrange', help='contig range 0v0', required=False)
    parser.add_argument('--min_sv_size', help='minimum SV ize to be detected .v. (default: 10000)',
                        required=False, type=int, default=10000)
    parser.add_argument('--output', help='output file path 0.<', required=False, default='output.vcf')
    parser.add_argument('--mapq', help='MAPQ filtering', required=False, type=int, default=10)
    parser.add_argument('--min_read_length', help='min read length for filtering (exclusive)',
                        required=False, type=int, default=0)
    parser.add_argument('--min_af', help='min af nya', required=False, type=float, default=0.2)
    parser.add_argument('--samtools', help='samtools executable', required=False, type=str, default='samtools')
    parser.add_argument('--header_file_path', help='header file path', required=True, type=str, default='')
    parser.add_argument('--processes', help='# of processes', required=False, type=int, default=40)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    test(args=parser.parse_args())
