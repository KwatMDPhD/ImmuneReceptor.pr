#ifndef _UTIL_H_
#define _UTIL_H_

#include "ic_common.h"
#include "ic_type.h"
#include "ic_io.h"
#include "ic_config.h"
#include "jansson.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

// #define MAX2(a,b)       ((a < b) ?  (b) : (a))
// #define MAX3(a,b,c)     ((MAX2(a,b) < c) ?  (c) : (MAX2(a,b)))
// #define MAX4(a,b,c,d)   ((MAX3(a,b, c) < d) ?  (d) : (MAX3(a,b,c)))
// #define MIN2(a,b)       ((a > b) ?  (b) : (a))
// #define MIN3(a,b,c)     ((MIN2(a,b) > c) ?  (c) : (MIN2(a,b)))
// #define MIN4(a,b,c,d)   ((MIN3(a,b, c) > d) ?  (d) : (MIN3(a,b,c)))


#define safe_free(p) safer_free((void**)&(p))
#define icalloc(n, type)    (type *)_icalloc(n, sizeof(type))
#define dbg_print(fmt, ...) \
	do { if (DEBUG) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, \
		                __LINE__, __func__, __VA_ARGS__); } while (0)



void safer_free(void **pp);
void *_icalloc(long number, int size);
void die(const char* format, ...);
void debug(char *format, ...);
void system_call(char *cmd);

//compute fisher exact test
double gammln(double xx);
double factln(int n);
double fisher_exact_test(int n11, int n12, int n21, int n22);

double compute_enrichment_pvalue(int population, int sub_population_with_feature, int sample_population, int sample_population_with_feature);
double hypergeometric_pmf(int x, int m, int n, int k);
double log_binomial_coefficient(int population, int sample);

//string manuplication
void ic_to_upper(char s[]);
void ic_to_lower(char s[]);
bool ic_starts_with(const char* str, const char* pattern);
bool ic_ends_with(const char *str, const char *pattern);
void ic_split_string(const char* istring, char delim, s_vector *array);
bool ic_all_space(const char *str);
char* ic_strmove(char *tgt, char *src);
char* ic_strtriml(char *str);
char* ic_strtrimr(char *str);
char* ic_strip(char* str);
char* ic_generate_random_string(int size);

void read_hla(parameter_t *global_options);
void read_target_cdr3(parameter_t *global_options, t_vector *vector);
void read_refer_cdr3(const char* ifile, c_vector *vector);
char* purge_cdr3(parameter_t *global_options, const char* chain, const char* cdr3);
void read_cluster_cdr3(parameter_t *global_options, const char* ifile, const char* chain, s_set *set);
int count_unique_target_cdr3(parameter_t *global_options, const char* chain);
int count_unique_refer_cdr3(parameter_t *global_options, const char* chain);
bool match_pattern(const char* cdr3, const char* pattern, int ignored_left, int ignored_rght);

void output_hla_association(parameter_t *global_options, const char* chain);
void output_cluster(parameter_t *global_options, const char* chain);
void prepare_data(parameter_t *global_options);
void format_input_cdr3_file(parameter_t *global_options);

void sort_cluster_on_entfrc(const char* filename);
void combinations(i_vector *iterable, int r, iv_vector *rslt);
void get_v_gene(const char* pattern, char* vgene);
bool within(i_vector *small, i_vector *large);
int issubset(i_vector *v1, i_vector *v2);
float entropy_fraction(s_vector *ivec);
int purge_cluster(parameter_t * global_options, s_vector *ivec, ii_hash *dict);
void create_overlap_matrix(parameter_t *global_options, p_vector *pvec, vi_vector *matrix);
int overlap_cluster(si_hash *h1, pattern_t *p2);
void connected_components(vi_vector *matrix, vi_vector *rslt);
void output_one_cluster(parameter_t *global_options, svi_hash *cdr3toindex, const char* chain, pattern_t *p, int lnum, FILE *ofp);
void output_html(parameter_t *global_options, const char* chain, p_vector *pvec, vi_vector *viv);

#endif
