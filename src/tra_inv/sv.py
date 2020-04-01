from collections import namedtuple
from utils import interval_from, merge_two_intervals, depth_stat_from

Breakpoint = namedtuple('Breakpoint', ['chr', 'position'])

class SV:
    def __init__(self, read, breakpoint_position1, breakpoint_position2, dp_result=None, flanking_bases=50):
        """
            SV Init.
            read:
                - one Read (supported read) from samtoolsViewViewer
            breakpoint_position1, breakpoint_position2:
                - two breakpoint position w.r.t. the reference
        """
        read1, read2 = read, read.SA_reads[0]

        self.breakpoint1, self.breakpoint2 = SV.sort_breakpoints(
            Breakpoint(read1.RNAME, breakpoint_position1), Breakpoint(read2.RNAME, breakpoint_position2)
        )
        self.type = read.SV_type
        self.supported_reads = [read]
        self.dp_results = [dp_result] if dp_result is not None else []

        self.breakpoint_intervals = [
            interval_from(self.breakpoint1.position, flanking_bases),
            interval_from(self.breakpoint2.position, flanking_bases),
        ]

    def __str__(self):
        chr1, pos1 = self.breakpoint1.chr, self.breakpoint1.position
        chr2, pos2 = self.breakpoint2.chr, self.breakpoint2.position

        return (
            "Point 1: %s:%d, Point 2: %s:%d. Type: %s. Number of reads supporting: %d" % (
                chr1, pos1, chr2, pos2, self.type, len(self.supported_reads)
            )
        )

    @property
    def first_supported_read(self):
        return self.supported_reads[0]

    @property
    def first_similar_score(self):
        return self.dp_results[0].similar_score

    @property
    def breakpoint_chromosomes(self):
        return self.breakpoint1.chr, self.breakpoint2.chr

    @classmethod
    def sort_breakpoints(cls, breakpoint1, breakpoint2):
        chr1, pos1 = breakpoint1.chr, breakpoint1.position
        chr2, pos2 = breakpoint2.chr, breakpoint2.position

        if chr2 < chr1 or (chr2 == chr1 and pos2 < pos1):
            return breakpoint2, breakpoint1
        else:
            return breakpoint1, breakpoint2

    def merge(self, sv):
        self.dp_results += sv.dp_results
        self.supported_reads = self.supported_reads + sv.supported_reads
        self.breakpoint_intervals = [
            merge_two_intervals(*(intervals)) for intervals in zip(self.breakpoint_intervals, sv.breakpoint_intervals)
        ]


class SV_Utilities:

    def __init__(self, bam, samtools):
        self.bam = bam
        self.samtools = samtools

    def allele_frequencies_from(self, sv):
        """get af from two breakpoints of sv"""
        for breakpoint in [sv.breakpoint1, sv.breakpoint2]:
            interval = interval_from(breakpoint.position, flanking_bases=50)
            mean_left, mean_right = depth_stat_from(self.bam, self.samtools, breakpoint, interval.begin, interval.end)
            af = len(sv.supported_reads) / max(max(mean_left, mean_right), 1)
            yield af

    def output_info_from(self, sv):
        """get QUAL and SCORE from sv"""
        af = sum(self.allele_frequencies_from(sv))
        score = 0
        for dp_result in sv.dp_results:
            score += dp_result.similar_score
        score = 1.0 * score / len(sv.dp_results)

        qual = str(int(score * 10) + int(min(af * 40 - 10, 30)))
        score = str(int(sv.first_similar_score * 10) + int(min(af * 40 - 10, 30)))

        return qual, score

    def vcf_rows_from(self, sv):
        """
        Get vcf-formatted rows of a sv, each array item is a string contains vcf rows for output

        If not inversion / translocation, return empty array []
        """
        # TODO - find representitive read
        read = sv.first_supported_read

        QUAL, _SCORE = self.output_info_from(sv)

        TYPE = read.SV_type
        chr1, chr2 = sv.breakpoint_chromosomes
        start, end = sv.breakpoint1.position, sv.breakpoint2.position

        if TYPE == "INV":
            return [f"{chr1}\t{start}\t.\tN\t<{TYPE}>\t{QUAL}\tPASS\tSVLEN=-{end - start};SVTYPE={TYPE};END={end}\tGT\t./."]

        elif TYPE == "TRA":
            read2 = read.SA_reads[0]
            cigar = read2.CigarString
            cigar_first = int(cigar.first_cigar_tuple[0]) if cigar.first_cigar_tuple[1] in "SH" else 0
            cigar_end = int(cigar.last_cigar_tuple[0]) if cigar.last_cigar_tuple[1] in "SH" else 0

            if cigar_first < cigar_end:
                if read2.is_forward_strand == read.is_forward_strand:
                    TYPE = "N[%s:%d[" % (chr2, end)
                    TYPE2 = "N[%s:%d[" % (chr1, start+1)
                    start2 = end-1
                else:
                    TYPE = "N]%s:%d]" % (chr2, end)
                    TYPE2 = "[%s:%d[N" % (chr1, start+1)
                    start2 = end+1
            else:
                if read2.is_forward_strand == read.is_forward_strand:
                    TYPE = "]%s:%d]N" % (chr2, end)
                    TYPE2 = "]%s:%d]N" % (chr1, start-1)
                    start2 = end+1
                else:
                    TYPE = "[%s:%d[N" % (chr2, end)
                    TYPE2 = "N]%s:%d]" % (chr1, start-1)
                    start2 = end-1

            return [
                f"{chr1}\t{start}\t.\tN\t{TYPE}\t{QUAL}\tPASS\tSVLEN=0;SVTYPE=BND\tGT\t./.",
                f"{chr2}\t{start2}\t.\tN\t{TYPE2}\t{QUAL}\tPASS\tSVLEN=0;SVTYPE=BND\tGT\t./.",
            ]

        return []
