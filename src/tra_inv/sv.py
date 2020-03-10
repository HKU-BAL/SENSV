from collections import defaultdict, namedtuple

Breakpoint = namedtuple('Breakpoint', ['chr', 'position'])
Interval = namedtuple('Interval', ['begin', 'end'])
DP_Result = namedtuple('DP_Result', ['output', 'similar_score'])


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
    def sort_breakpoints(self, breakpoint1, breakpoint2):
        chr1, pos1 = breakpoint1.chr, breakpoint1.position
        chr2, pos2 = breakpoint2.chr, breakpoint2.position

        if chr2 < chr1 or (chr2 == chr1 and pos2 < pos1):
            return breakpoint2, breakpoint1
        else:
            return breakpoint1, breakpoint2

    @classmethod
    def merge_two_intervals(self, interval1, interval2):
        return Interval(min(interval1.begin, interval2.begin), max(interval1.end, interval2.end))

    def merge(self, sv):
        self.dp_results += sv.dp_results
        self.supported_reads = self.supported_reads + sv.supported_reads
        self.breakpoint_intervals = [
            SV.merge_two_intervals(*(intervals)) for intervals in zip(self.breakpoint_intervals, sv.breakpoint_intervals)
        ]


def interval_from(position, flanking_bases):
    # intervals are [start, end)
    return Interval(max(position - flanking_bases, 0), position + flanking_bases + 1)
