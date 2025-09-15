#ifndef _MOTIF_H_
#define _MOTIF_H_

#include "ic_config.h"
#include "ic_common.h"
#include "ic_util.h"
#include "ic_type.h"

#define NUM_AA_SCORE_MATRIX 24
extern const char SCORE_MATRIX_AAS[];
extern const char BLOSUM62[];
extern const int BLOSUM62_SCORE_MATRIX[NUM_AA_SCORE_MATRIX][NUM_AA_SCORE_MATRIX];

int max3(int i, int j, int k);
void convert_pseq_to_i_vector(const char* pseq, i_vector *vec);
char* convert_i_vector_to_pseq(i_vector *vec);
int calculate_global_alignment_score(const char* seq1, const char* seq2);
void sort_cluster(si_vector *vec);

#endif
