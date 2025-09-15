#ifndef _ALL_H_
#define _ALL_H_

#include "ic_config.h"
#include "ic_common.h"
#include "ic_util.h"
#include "ic_type.h"

void count_pattern(parameter_t *global_options, t_vector *tvec, c_vector *rvec, int cdr3_len, const char* chain, unsigned long pattern_mask, av_hash *pvh);
bool match_all_pattern(char* pattern, char* cdr3);
void discover_enriched_pattern(parameter_t *global_options, const char* chain);
void do_all_clustering(parameter_t* global_options);

#endif
