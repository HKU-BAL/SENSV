#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dp.h"

#define BUFSIZE 500000

align_result result;

int main(int argc, char *argv[]) {
    if (argc != 1) {
        printf("usage: %s\n", argv[0]);
        return -1;
    }

    char buf[BUFSIZE+1];
    char *query_seq;
    char *ref_seq;

    fgets(buf, BUFSIZE, stdin);

    query_seq = strtok(buf, " ");
    ref_seq = strtok(NULL, " ");

    if (ref_seq[strlen(ref_seq)-1] == '\n') {
        ref_seq[strlen(ref_seq)-1] = 0;
    }

    dp(query_seq, ref_seq, &result, 1);
    fprintf(stdout, "%d %d %d %d %d %d %d %d %d %d %d %d", \
        result.score, result.query_start, result.query_end, result.ref_start, result.ref_end, \
        result.gap_start, result.gap_end, result.query_pos, result.global_ref_start, result.global_ref_end, \
        result.global_gap_start, result.global_gap_end);
}
