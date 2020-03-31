
CIGAR_CONSUME_QUERY_OPERATORS = set("MIS=X")
CIGAR_CONSUME_REFERENCE_OPERATORS = set("MDN=X")


class CigarString:

    def __init__(self, cigar=""):
        # self.cigar = cigar

        self.cigar_tuples, self.query_length, self.ref_length = CigarString.parse(cigar)

    @classmethod
    def parse(self, cigar):
        cigar_tuples = []

        count_query = 0
        count_reference = 0

        advance = 0
        for c in str(cigar):
            if c.isdigit():
                advance = advance * 10 + int(c)
                continue

            cigar_tuples.append((advance, c))

            if c in CIGAR_CONSUME_QUERY_OPERATORS:
                count_query += advance

            if c in CIGAR_CONSUME_REFERENCE_OPERATORS:
                count_reference += advance

            # reset advance
            advance = 0

        # only first and last cigar tuple is used
        cigar_tuples = [cigar_tuples[0], cigar_tuples[-1]]

        return cigar_tuples, count_query, count_reference

    @property
    def first_cigar_tuple(self):
        if len(self.cigar_tuples) <= 0:
            return (0, "")
        return self.cigar_tuples[0]

    @property
    def last_cigar_tuple(self):
        if len(self.cigar_tuples) <= 0:
            return (0, "")
        return self.cigar_tuples[-1]

    def clipped_length(self, backward=False):
        cigar_tuple_count, cigar_tuple_operation = self.first_cigar_tuple if not backward else self.last_cigar_tuple

        return cigar_tuple_count if cigar_tuple_operation in "SH" else 0
