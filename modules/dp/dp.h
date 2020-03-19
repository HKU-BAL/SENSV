#include <stdint.h>

typedef struct {
  int score;
  int query_start;
  int query_end;
  int ref_start;
  int ref_end;
  int gap_start;
  int gap_end;
  int query_pos;
  int global_ref_start;
  int global_ref_end;
  int global_gap_start;
  int global_gap_end;

  int n_cigar;
  uint32_t * cigar; // "MIDN"
} align_result;

void dp(char *query, char *ref, align_result *result, int quota);
