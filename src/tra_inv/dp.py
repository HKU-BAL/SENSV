from dp_module import get_max_array

verbose = False


REVERSE_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


def reverse_complement_of(sequence):
    return "".join(reversed([REVERSE_COMPLEMENT[base] for base in sequence]))


class DP:

    def __init__(self, query, ref1, ref2, ref1_start=0, ref2_start=0, reverse=False, reverse_ref=False):
        if reverse_ref:
            self.ref1 = ref2
            self.ref2 = ref1
        else:
            self.ref1 = ref1
            self.ref2 = ref2

        self.query = query.upper()
        self.ref1_start = ref1_start
        self.ref2_start = ref2_start
        self.reverse_ref = reverse_ref
        self.reverse = reverse
        self.max_array1 = []
        self.max_array2 = []

    def first_alignment_arguments(self):
        if not self.reverse_ref:
            return self.query, self.ref1
        else:
            if self.reverse:
                return reverse_complement_of(self.query)[::-1], self.ref1[::-1]
            else:
                return self.query, self.ref1

    def second_alignment_arguments(self):
        if not self.reverse_ref:
            if self.reverse:
                return reverse_complement_of(self.query), self.ref2
            else:
                return self.query[::-1], self.ref2[::-1]
        else:
            return self.query[::-1], self.ref2[::-1]

    def get_read_breakpoint(self):
        first_alignment_arguments = self.first_alignment_arguments()
        second_alignment_arguments = self.second_alignment_arguments()

        if verbose:
            print(*(first_alignment_arguments))
            print(*(second_alignment_arguments))

        first_align = get_max_array(*(first_alignment_arguments))
        second_align = get_max_array(*(second_alignment_arguments))

        if len(first_align) != len(second_align):
            print("... dp: Error in getting breakpoint ...")
            return None

        second_align = second_align[::-1]

        self.max_array1 = first_align
        self.max_array2 = second_align

        combined = []
        for first_align_tuple, second_align_tuple in zip(first_align, second_align):
            if verbose:
                print(first_align_tuple, second_align_tuple, first_align_tuple[0] + second_align_tuple[0])
            combined.append(first_align_tuple[0] + second_align_tuple[0])

        max_index = combined.index(max(combined))
        bp1 = first_align[max_index][1]
        bp2 = second_align[max_index][1]

        if self.reverse:
            if not self.reverse_ref:
                bp2 = len(self.ref2) + 1 - bp2
            else:
                bp1 = len(self.ref1) + 1 - bp1

        if self.reverse_ref:
            bp1, bp2 = bp2, bp1

        self.bp1, self.bp2 = bp1, bp2

        if verbose:
            print(
                "ref1", self.ref1_start, "ref2", self.ref2_start,
                "len", len(self.ref1), "max", max_index, "bp", bp1, bp2
            )

        if self.reverse_ref:
            return (max_index, max(combined), self.ref1_start+len(self.ref1)-bp1, self.ref2_start+bp2)
        else:
            return (max_index, max(combined), self.ref1_start+bp1, self.ref2_start+len(self.ref2)-bp2)

    def get_similarity_score(self, max_index):
        ref1_align_score = self.max_array1[max_index][0]
        ref2_align_score = self.max_array2[max_index][0]

        ref2_align_to_ref1 = self.max_array1[min(max_index+20000, len(self.max_array1)-1)][0] - ref1_align_score
        ref1_align_to_ref2 = self.max_array2[max(max_index-20000, 0)][0] - ref2_align_score

        # ref1_score = 1.0 * (ref2_align_score - ref2_align_to_ref1) / (min(20000, len(self.max_array1)-max_index)+1)
        # ref2_score = 1.0 * (ref1_align_score - ref1_align_to_ref2) / (min(20000, max_index)+1)

        score = 1.0 * (ref2_align_score - ref2_align_to_ref1 + ref1_align_score - ref1_align_to_ref2) / \
            (min(20000, len(self.max_array1)-max_index) + min(20000, max_index) + 1)

        # if self.reverse_ref:
        #     print(
        #         self.ref1_start+len(self.ref1)-self.bp1, self.ref2_start+self.bp2, ref1_align_score,
        #         ref2_align_to_ref1, ref1_score, ref2_align_score, ref1_align_to_ref2, ref2_score, max_index
        #     )
        # else:
        #     print(
        #         self.ref1_start+self.bp1, self.ref2_start+len(self.ref2)-self.bp2, ref1_align_score,
        #         ref2_align_to_ref1, ref1_score, ref2_align_score, ref1_align_to_ref2, ref2_score, max_index
        #     )

        return score
