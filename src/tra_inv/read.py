from collections import namedtuple

from cigar_parser import CigarString
from utils import major_contigs_set

SA_Read = namedtuple('SA_Read', [
    'RNAME',
    'POS',
    'is_forward_strand',
    'MAPQ',
    'CigarString',
])


class Read:

    def __init__(self, read):
        self.QNAME, self.FLAG, self.FLAGS, self.RNAME, self.POS, self.MAPQ, CIGAR, \
            self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, _QUAL, self.tags = Read.parse(read)

        self._SA_reads = None
        self._SV_type = None
        self.CigarString = CigarString(CIGAR)

    def __str__(self):
        string = (
            "QNAME: %s\nPOS: %s:%s\nFLAG: %s\nTAG:\n" % (self.QNAME, self.RNAME, self.POS, self.FLAG) +
            "".join(["%s: %s\n" % (key, value) for (key, value) in self.tags.items()])
        )

        return string

    @classmethod
    def parse(cls, read):
        columns = read.strip().split()

        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL = columns[:11]
        FLAGS = Read.parse_flags_from(FLAG)
        TAGS = Read.parse_tags_from(columns)

        return QNAME, FLAG, FLAGS, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, TAGS

    @classmethod
    def parse_flags_from(cls, flag_str):
        flag_int_value = int(flag_str)
        return [a == a & flag_int_value for a in [pow(2, n) for n in range(12)]]

    @classmethod
    def parse_tags_from(cls, columns):
        tags = {}
        for tag in columns[11:]:
            tag_list = tag.split(":")
            tags[tag_list[0]] = tag_list[1:]

        return tags

    @property
    def SA_reads(self):
        if self._SA_reads is not None:
            return self._SA_reads

        """
            SA: Other canonical alignments in a chimeric alignment

            SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
            Other canonical alignments in a chimeric alignment, formatted as a semicolon-delimited list.
            Each element in the list represents a part of the chimeric alignment. Conventionally, at a
            supplementary line, the first element points to the primary line. Strand is either ‘+’ or ‘-’,
            indicating forward/reverse strand, corresponding to FLAG bit 0x10. Pos is a 1-based coordinate
        """
        SA_reads = []
        for SA_str in self.tags["SA"][1].strip().split(";"):
            SA_values = SA_str.split(',')
            if len(SA_values) != 6:
                continue

            RNAME, POS, STRAND, CIGAR, MAPQ, _NM = tuple(SA_values)

            SA_reads.append(SA_Read(
                RNAME=RNAME,
                POS=POS,
                is_forward_strand=STRAND == "+",
                MAPQ=int(MAPQ),
                CigarString=CigarString(CIGAR),
            ))

        # only use the first SA read
        self._SA_reads = [SA_reads[0]]
        return self._SA_reads

    @property
    def is_forward_strand(self):
        return self.FLAGS[4] is not True

    @property
    def split_position(self):
        read1 = self
        read2 = self.SA_reads[0]

        read1_direction = read1.CigarString.clipped_length() < read1.CigarString.clipped_length(backward=True)
        read2_relative_strand = (read1.is_forward_strand != read2.is_forward_strand) == read1_direction

        return (
            int(read1.POS) + (read1.CigarString.ref_length) * read1_direction,
            int(read2.POS) + (read2.CigarString.ref_length) * read2_relative_strand
        )

    @property
    def SV_type(self):
        if self._SV_type is not None:
            return self._SV_type

        read1 = self
        read2 = self.SA_reads[0]

        read1_ref = int(read1.POS)
        read2_ref = int(read2.POS)

        read1_strand = read1.is_forward_strand
        read2_strand = read2.is_forward_strand

        read1_query_adv = read1.CigarString.clipped_length(backward=read1_strand)
        read2_query_adv = read2.CigarString.clipped_length(backward=read2_strand)

        if read1.RNAME not in major_contigs_set or read2.RNAME not in major_contigs_set:
            self._SV_type = ""
        elif read1.RNAME == read2.RNAME:
            arr = ["DEL", "INV", "INV", "DUP", "DUP", "INV", "INV", "DEL"]
            if read1_query_adv < read2_query_adv:
                self._SV_type = arr[read1_strand + 2*read2_strand + 4*(read1_ref > read2_ref)]
            else:
                self._SV_type = arr[read2_strand + 2*read1_strand + 4*(read2_ref > read1_ref)]
        else:
            self._SV_type = "TRA"

        return self._SV_type

    def SV_type_with_mapq(self, mapq_filter):
        read2 = self.SA_reads[0]
        if read2.MAPQ < mapq_filter:
            return ""
        return self.SV_type
