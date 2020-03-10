import sys
import os
from multiprocessing import Pool
from argparse import ArgumentParser
from collections import defaultdict
from intervaltree import IntervalTree

from sv import SV, DP_Result, interval_from
from dp import DP
from utils import get_contig_list, get_reference_from_fasta, depth_stat_from
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


def filtered_reads_with(reads, mapq, min_sv_size, min_read_length):
    def get_SV_size(read, min_sv_size):
        chr1 = read.RNAME

        read2 = read.SA_reads[0]
        chr2 = read2.RNAME
        if chr1 != chr2:
            return min_sv_size

        pos1, pos2 = read.split_position
        return abs(pos1 - pos2)

    filtered_reads = {}
    for read in reads:
        if not read.SV_type_with_mapq(mapq_filter=mapq) in ["INV", "TRA"]:
            continue
        if len(read.SEQ) < min_read_length:
            continue
        if read.CigarString.first_cigar_tuple[1] == "H" or read.CigarString.last_cigar_tuple[1] == "H":
            continue
        if get_SV_size(read, min_sv_size) < args.min_sv_size:
            continue
        if read.QNAME in filtered_reads and len(filtered_reads[read.QNAME].SEQ) >= len(read.SEQ):
            continue

        filtered_reads[read.QNAME] = read

    return list(filtered_reads.values())


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


def get_dp_result(read):
    read2 = read.SA_reads[0]

    cigar = read.CigarString
    cigar_first = int(cigar.cigar_tuples[0][0]) if cigar.cigar_tuples[0][1] in "SH" else 0
    cigar_end = int(cigar.cigar_tuples[-1][0]) if cigar.cigar_tuples[-1][1] in "SH" else 0
    pos1, pos2 = read.split_position
    interval1, interval2 = interval_from(pos1, flanking_bases=20000), interval_from(pos2, flanking_bases=20000)

    ref1 = get_reference_from_fasta(args.ref, read.RNAME, interval1.begin, interval1.end)
    ref2 = get_reference_from_fasta(args.ref, read2.RNAME, interval2.begin, interval2.end)

    dp = DP(read.SEQ, ref1, ref2, interval1.begin, interval2.begin,
            read.is_forward_strand != read2.is_forward_strand, cigar_end < cigar_first)

    dp_result = dp.get_read_breakpoint()
    dp_similar_score = dp.get_similarity_score(dp_result[0])

    return (read, dp_result, dp_similar_score)


def dp_result_from(reads, processes):
    print("... Performing dp on %d reads ..." % len(reads))
    pool = Pool(processes=processes)
    result = pool.map(get_dp_result, reads)
    pool.close()

    return result


def sv_candidates_from(dp_result):
    interval_trees = defaultdict(IntervalTree)
    for (read, output, similar_score) in dp_result:
        if similar_score < 1.0:
            continue
        # print(output, similar_score)

        sv_candidate = SV(
            read=read,
            breakpoint_position1=output[2],
            breakpoint_position2=output[3],
            dp_result=DP_Result(output, similar_score),
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


def is_valid_sv_with_args(sv, ref, bam, samtools):
    if len(sv.supported_reads) <= 1:
        return False

    for chr, breakpoint_interval in zip(sv.breakpoint_chromosomes, sv.breakpoint_intervals):
        if "N" in get_reference_from_fasta(ref, chr, breakpoint_interval.begin, breakpoint_interval.end):
            return False

    for breakpoint in [sv.breakpoint1, sv.breakpoint2]:
        interval = interval_from(breakpoint.position, flanking_bases=50)
        mean_left, mean_right = depth_stat_from(bam, samtools, breakpoint, interval.begin, interval.end)
        af = len(sv.supported_reads) / max(max(mean_left, mean_right), 1)
        if af < 0.2:
            return False

    return True


def final_sv_candidates_from(sv_candidates, ref, bam, samtools, processes):
    pool = Pool(processes=processes)
    is_valid_sv_list = pool.starmap(
        is_valid_sv_with_args,
        [(sv_candidate, ref, bam, samtools) for sv_candidate in sv_candidates]
    )
    pool.close()

    final_sv_candidates = []
    for is_valid_sv, sv_candidate in zip(is_valid_sv_list, sv_candidates):
        if is_valid_sv:
            final_sv_candidates.append(sv_candidate)
    return final_sv_candidates


def output(sv_list):
    from pathlib import Path
    repo_directory = Path(os.path.dirname(__file__)).parent.parent

    with open(args.output, 'w') as vcf:
        header = open(repo_directory / 'data' / 'output' / 'header').read().splitlines()
        for line in header:
            print(line, file=vcf)

        for sv in sv_list:
            # log
            # print(sv)
            # for dp_result in sv.dp_results:
            #     print(dp_result.output, dp_result.similar_score)

            # TODO - find representitive read
            read = sv.first_supported_read
            QUAL = str(int(sv.first_similar_score * 10))
            TYPE = read.SV_type
            chr1, chr2 = sv.breakpoint_chromosomes
            start, end = sv.breakpoint1.position, sv.breakpoint2.position

            if TYPE == "INV":
                print("%s\t%d\t.\tN\t<%s>\t%s\tPASS\tSVLEN=-%d;SVTYPE=%s;END=%d\tGT\t./." %
                      (chr1, start, TYPE, QUAL, end - start, TYPE, end), file=vcf)
            elif TYPE == "TRA":
                read2 = read.SA_reads[0]
                cigar = read2.CigarString
                cigar_first = int(cigar.first_cigar_tuple[0]) if cigar.first_cigar_tuple[1] in "SH" else 0
                cigar_end = int(cigar.last_cigar_tuple[0]) if cigar.last_cigar_tuple[1] in "SH" else 0

                if cigar_first < cigar_end:
                    TYPE = f'N[{chr2}:{end}['
                else:
                    TYPE = f']{chr2}:{end}]N'

                print("%s\t%d\t.\tN\t%s\t%s\tPASS\tSVLEN=0;SVTYPE=BND\tGT\t./." % (chr1, start, TYPE, QUAL), file=vcf)


def test(args):
    reads = reads_from_samtools(args.samtools, args.bam, args.ctgname, args.ctgrange, args.mapq)
    reads = filtered_reads_with(reads, args.mapq, args.min_sv_size, args.min_read_length)
    reads = filtered_reads_before_dp(reads)

    result = dp_result_from(reads, args.processes)

    sv_candidates = sv_candidates_from(result)
    final_sv_candidates = final_sv_candidates_from(sv_candidates, args.ref, args.bam, args.samtools, args.processes)

    output(final_sv_candidates)


if __name__ == "__main__":

    parser = ArgumentParser(description='A script to see see what split reads look like >.0')

    parser.add_argument('--bam', help='The bam file to be read :D', required=True)
    parser.add_argument('--ref', help='The corresponding fasta file >.0', required=True)
    parser.add_argument('--ctgname', help='contig name 0.0', required=False)
    parser.add_argument('--ctgrange', help='contig range 0v0', required=False)
    parser.add_argument('--min_sv_size', help='minimum SV ize to be detected .v. (default: 10000)',
                        required=False, type=int, default=10000)
    parser.add_argument('--output', help='output file name 0.<', required=False, default='output.vcf')
    parser.add_argument('--mapq', help='MAPQ filtering', required=False, type=int, default=10)
    parser.add_argument('--min_read_length', help='min read length for filtering (exclusive)',
                        required=False, type=int, default=0)
    parser.add_argument('--samtools', help='samtools executable', required=False, type=str, default='samtools')
    parser.add_argument('--processes', help='# of processes', required=False, type=int, default=40)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    test(args)
