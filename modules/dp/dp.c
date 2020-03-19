#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dp.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void dp(char *query, char *ref, align_result *result, int quota);

typedef struct s_cost {
    int gap_start;
    int gap_end;
    int query_pos;
    float cost;
} Cost;

void dp(char *query, char *ref, align_result *result, int quota) {
    const int isDebug = 0;
    //const int quota = 1; // how many times 'free gap penalty' can be used
    const float NEGINF = -999999.0;

    if (isDebug) {
        fprintf(stderr, "[INFO] DP Called\n");
    }

    int i, j, k, l, z, v, gapStart, gapEnd;
    float O1 = -4.0, O2 = -24.0, E1 = -2.0, E2 = -1.0;
    int rlen, qlen;

    Cost ***M;
    Cost ***D;
    Cost ***I;
    Cost ***FD;
    Cost new_cost;

    float ***D_len;
    float ***I_len;
    float score1, score2;

    char refAlign[63000];
    char queryAlign[63000];

    int MA = 2;
    int MM = -4;
    int GO = -4;
    int GE = -2;

    float maxM = 0;
    int maxI = 0;
    int maxJ = 0;
    int maxZ = 0;
    int maxM_gap_start = -1;
    int maxM_gap_end = -1;
    int maxM_query_pos = -1;

    float deletion_cost;
    float insertion_cost;
    float free_deletion_cost;
    float match_mismatch_cost;

    rlen = strlen(ref);
    qlen = strlen(query);

    // alloc
    // if (isDebug) {
    //     printf("[INFO] Before Alloc\n");
    // }

    size_t SIZE_I_COST, SIZE_J_COST, SIZE_Z_COST;
    size_t SIZE_J_FLOAT, SIZE_I_FLOAT, SIZE_Z_FLOAT;

    SIZE_Z_COST = (1+quota) * sizeof(Cost **);
    SIZE_Z_FLOAT = (1+quota) * sizeof(float **);
    SIZE_I_COST = (1+1) * sizeof(Cost *);
    SIZE_I_FLOAT = (1+1) * sizeof(float *);
    SIZE_J_COST = (1+rlen) * sizeof(Cost);
    SIZE_J_FLOAT = (1+rlen) * sizeof(float);

    M = (Cost ***)malloc(SIZE_Z_COST);
    D = (Cost ***)malloc(SIZE_Z_COST);
    I = (Cost ***)malloc(SIZE_Z_COST);
    FD = (Cost ***)malloc(SIZE_Z_COST);
    D_len = (float ***)malloc(SIZE_Z_FLOAT);
    I_len = (float ***)malloc(SIZE_Z_FLOAT);

    for (z=0; z < 1+quota; ++z) {
        M[z] = (Cost **)malloc(SIZE_I_COST);
        D[z] = (Cost **)malloc(SIZE_I_COST);
        I[z] = (Cost **)malloc(SIZE_I_COST);
        FD[z] = (Cost **)malloc(SIZE_I_COST);
        D_len[z] = (float **)malloc(SIZE_I_FLOAT);
        I_len[z] = (float **)malloc(SIZE_I_FLOAT);

        for (i=0; i < 1+1; ++i) {
            M[z][i] = (Cost *)malloc(SIZE_J_COST);
            D[z][i] = (Cost *)malloc(SIZE_J_COST);
            I[z][i] = (Cost *)malloc(SIZE_J_COST);
            FD[z][i] = (Cost *)malloc(SIZE_J_COST);
            D_len[z][i] = (float *)malloc(SIZE_J_FLOAT);
            I_len[z][i] = (float *)malloc(SIZE_J_FLOAT);

            if (z == 0 && i == 0) {
                for (j=0; j<=rlen; ++j) {
                    // D[z][i][j] = (Cost) { .gap_end = -1, .cost = NEGINF + 0.0 };
                    // M[z][i][j] = (Cost) { .gap_end = -1, .cost = 0.0 };
                    D[z][i][j] = (Cost) { .gap_start = -1, .gap_end = -1, .query_pos = -1, .cost = NEGINF + 0.0 };
                    M[z][i][j] = (Cost) { .gap_start = -1, .gap_end = -1, .query_pos = -1, .cost = 0.0 };
                    D_len[z][i][j] = NEGINF + 0.0;
                }
            } else {
                memcpy(D[z][i], D[0][0], SIZE_J_COST);
                memcpy(M[z][i], M[0][0], SIZE_J_COST);
                memcpy(D_len[z][i], D_len[0][0], SIZE_J_FLOAT);
            }

            memcpy(I[z][i], D[z][i], SIZE_J_COST);
            memcpy(FD[z][i], D[z][i], SIZE_J_COST);
            memcpy(I_len[z][i], D_len[z][i], SIZE_J_FLOAT);
        }
    }

    // if (isDebug) {
    //     printf("[INFO] After Alloc\n");
    // }

    int prev_i, curr_i;

    for (i=1; i<=qlen; ++i) {
        prev_i = (i + 1) % 2;
        curr_i = (i % 2);

        for (j=1; j<=rlen; ++j) {

            for (z=0; z<=quota; ++z) {
                // min(O1 + k * E1, O2 + k * E2) (length of the gap)
                // for O1 = 4, O2 = 24, E1 = 2, E2 = 1
                // k <= 20 will choose O1 + k * E1
                // k > 20 will choose O2 + k * E2

                // insertion
                score1 = M[z][prev_i][j].cost + GO + E1;
                if (I_len[z][prev_i][j] + 1 <= 20) {
                    score2 = I[z][prev_i][j].cost + E1;
                } else {
                    score2 = I[z][prev_i][j].cost + E2;
                }
                if (score1 >= score2) {
                    I_len[z][curr_i][j] = 1;
                    I[z][curr_i][j].cost = score1;
                } else {
                    I_len[z][curr_i][j] = I_len[z][prev_i][j] + 1;
                    I[z][curr_i][j].cost = score2;
                }
                I[z][curr_i][j].gap_start = M[z][prev_i][j].gap_start;
                I[z][curr_i][j].gap_end = M[z][prev_i][j].gap_end;
                I[z][curr_i][j].query_pos = M[z][prev_i][j].query_pos;

                // deletion
                score1 = M[z][curr_i][j-1].cost + GO + E1;
                if (D_len[z][curr_i][j-1] + 1 <= 20) {
                    score2 = D[z][curr_i][j-1].cost + E1;
                } else {
                    score2 = D[z][curr_i][j-1].cost + E2;
                }
                if (score1 >= score2) {
                    D_len[z][curr_i][j] = 1;
                    D[z][curr_i][j].cost = score1;
                } else {
                    D_len[z][curr_i][j] = D_len[z][curr_i][j-1] + 1;
                    D[z][curr_i][j].cost = score2;
                }
                D[z][curr_i][j].gap_start = M[z][curr_i][j-1].gap_start;
                D[z][curr_i][j].gap_end = M[z][curr_i][j-1].gap_end;
                D[z][curr_i][j].query_pos = M[z][curr_i][j-1].query_pos;

                // free deletion
                if (z >= 1) {
                    // score1: use free deletion at this position
                    // score2: continue to use the free deletion at this position
                    score1 = M[z-1][curr_i][j-1].cost;
                    score2 = FD[z][curr_i][j-1].cost;
                    if (score1 >= score2) {
                        FD[z][curr_i][j].cost = score1;
                        FD[z][curr_i][j].gap_start = j;
                        FD[z][curr_i][j].query_pos = i;
                    } else {
                        FD[z][curr_i][j].cost = score2;
                        FD[z][curr_i][j].gap_start = FD[z][curr_i][j-1].gap_start;
                        FD[z][curr_i][j].query_pos = FD[z][curr_i][j-1].query_pos;
                    }
                    FD[z][curr_i][j].gap_end = j;
                    // FD[z][curr_i][j] = MAX(M[z-1][curr_i][j-1], FD[z][curr_i][j-1]);
                }

                deletion_cost = D[z][curr_i][j].cost;
                insertion_cost = I[z][curr_i][j].cost;
                free_deletion_cost = FD[z][curr_i][j].cost;
                match_mismatch_cost = M[z][prev_i][j-1].cost + (query[i-1] == ref[j-1] ? MA : MM);

                new_cost.cost = 0;
                new_cost.gap_start = -1;
                new_cost.gap_end = -1;
                if (new_cost.cost <= match_mismatch_cost) {
                    new_cost.cost = match_mismatch_cost;
                    new_cost.gap_start = M[z][prev_i][j-1].gap_start;
                    new_cost.gap_end = M[z][prev_i][j-1].gap_end;
                    new_cost.query_pos = M[z][prev_i][j-1].query_pos;
                }
                if (new_cost.cost <= free_deletion_cost) {
                    new_cost.cost = free_deletion_cost;
                    new_cost.gap_start = FD[z][curr_i][j].gap_start;
                    new_cost.gap_end = FD[z][curr_i][j].gap_end;
                    new_cost.query_pos = FD[z][curr_i][j].query_pos;
                }
                if (new_cost.cost <= insertion_cost) {
                    new_cost.cost = insertion_cost;
                    new_cost.gap_start = I[z][curr_i][j].gap_start;
                    new_cost.gap_end = I[z][curr_i][j].gap_end;
                    new_cost.query_pos = I[z][curr_i][j].query_pos;
                }
                if (new_cost.cost <= deletion_cost) {
                    new_cost.cost = deletion_cost;
                    new_cost.gap_start = D[z][curr_i][j].gap_start;
                    new_cost.gap_end = D[z][curr_i][j].gap_end;
                    new_cost.query_pos = D[z][curr_i][j].query_pos;
                }

                M[z][curr_i][j].cost = new_cost.cost;
                M[z][curr_i][j].gap_start = new_cost.gap_start;
                M[z][curr_i][j].gap_end = new_cost.gap_end;
                M[z][curr_i][j].query_pos = new_cost.query_pos;

                if (new_cost.cost > maxM) {
                    maxI = i;
                    maxJ = j;
                    maxZ = z;
                    maxM = new_cost.cost;

                    maxM_gap_start = new_cost.gap_start;
                    maxM_gap_end = new_cost.gap_end;
                    maxM_query_pos = new_cost.query_pos;
                }
            }
        }
    }

    if (isDebug) {
        fprintf(stderr,
            "maxM=%.0f, maxI=%d, maxJ=%d, maxZ=%d\ngap_start: %d, gap_end: %d, query_pos: %d\n",
            maxM, maxI, maxJ, maxZ,
            maxM_gap_start, maxM_gap_end, maxM_query_pos
        );
    }

    // printf("%s,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s\n", \
    //     queryName, maxM, i, qlen-maxI, j, gapStart, gapEnd, \
    //     query, ref, queryAlign, refAlign);

    // free
    for (z=0; z < 1+quota; ++z) {
        for (i=0; i < 1+1; ++i) {
            free(D[z][i]);
            free(FD[z][i]);
            free(I[z][i]);
            free(M[z][i]);
            free(D_len[z][i]);
            free(I_len[z][i]);
        }
        free(D[z]);
        free(FD[z]);
        free(I[z]);
        free(M[z]);
        free(D_len[z]);
        free(I_len[z]);
    }
    free(D);
    free(FD);
    free(I);
    free(M);
    free(D_len);
    free(I_len);

    result->score = maxM;
    result->query_start = -1;
    result->query_end = maxI;
    result->ref_start = -1;
    result->ref_end = maxJ;
    result->gap_start = maxM_gap_start;
    result->gap_end = maxM_gap_end;
    result->query_pos = maxM_query_pos;
}
