#ifndef _TYPE_H_
#define _TYPE_H_

#include "ic_common.h"
#include "ic_type.h"

// vector containers
typedef kvec_t (int) i_vector;
void destruct_i_vector(i_vector *array);
void sort_i_vector(i_vector *array);
int int_compar(const void* a, const void* b);

typedef kvec_t (i_vector *) iv_vector;
void destruct_iv_vector(iv_vector *array);

typedef kvec_t (float) f_vector;
void destruct_f_vector(f_vector *array);
typedef kvec_t (char *) s_vector;
void destruct_s_vector(s_vector *array);
void destruct_s_vector_nd(s_vector *array);
void sort_s_vector(s_vector *array);
int str_compar(const void* a, const void* b);
void sort_s_vector_by_length(s_vector *array);
int str_length_compar(const void* a, const void* b);

typedef kvec_t (i_vector *) vi_vector;
void destruct_vi_vector(vi_vector *array);

// map containers
KHASH_SET_INIT_STR(set)
typedef khash_t (set) s_set;
void append_to_s_set_nd(s_set *sset, const char* key);
void append_to_s_set(s_set *sset, const char* key);
void destruct_s_set_nd(s_set *sset);
void destruct_s_set(s_set *sset);
bool key_exist_s_set(s_set *hash, const char *key);
char* get_value_s_set(s_set *hash, const char *key);

KHASH_MAP_INIT_STR(sim, int)
typedef khash_t (sim) si_hash;
void destruct_si_hash(si_hash *hash);
void add_to_si_hash(si_hash *hash, char* key);
void destruct_si_hash_nd(si_hash *hash);
void add_to_si_hash_nd(si_hash *hash, char* key);
void append_to_si_hash(si_hash *hash, const char *key, int val);
void append_to_si_hash_nd(si_hash *hash, const char *key, int val);
bool key_exist_si_hash(si_hash *hash, const char *key);
int get_value_si_hash(si_hash *hash, const char* key, const int default_value);

KHASH_MAP_INIT_STR(sfm, float)
typedef khash_t (sfm) sf_hash;
void destruct_sf_hash(sf_hash *hash);
void append_to_sf_hash(sf_hash *hash, const char *key, float val);
bool key_exist_sf_hash(sf_hash *hash, const char *key);
float get_value_sf_hash(sf_hash *hash, const char *key, float default_value);

KHASH_MAP_INIT_INT(ifm, float)
typedef khash_t (ifm) if_hash;
void destruct_if_hash(if_hash *hash);
bool key_exist_if_hash(if_hash *hash, int key);

KHASH_MAP_INIT_INT(iim, int)
typedef khash_t (iim) ii_hash;
void append_to_ii_hash(ii_hash *hash, int key, int val);
void add_to_ii_hash(ii_hash *hash, int key, int val);
bool key_exist_ii_hash(ii_hash *hash, int key);
void destruct_ii_hash(ii_hash *hash);

KHASH_MAP_INIT_INT64(lim, int)
typedef khash_t (lim) li_hash;
void append_to_li_hash(li_hash *hash, unsigned long key, int val);
void add_to_li_hash(li_hash *hash, unsigned long key, int val);
bool key_exist_li_hash(li_hash *hash, unsigned long key);
void destruct_li_hash(li_hash *hash);


KHASH_MAP_INIT_STR(svsm, s_vector*)
typedef khash_t (svsm) svs_hash;
void append_to_svs_hash(svs_hash *hash, const char *key, const char* val);
void destruct_svs_hash(svs_hash *hash);
bool key_exist_svs_hash(svs_hash *hash, const char *key);
s_vector* get_value_svs_hash(svs_hash* hash, const char* key);
void append_to_svs_hash_nd(svs_hash *hash, const char *key, const char* val);
void destruct_svs_hash_nd(svs_hash *hash);

KHASH_MAP_INIT_INT(ivim, i_vector*)
typedef khash_t (ivim) ivi_hash;
void append_to_ivi_hash(ivi_hash *hash,  int key, int val);
void destruct_ivi_hash(ivi_hash *hash);
bool key_exist_ivi_hash(ivi_hash *hash, int key);
i_vector* get_value_ivi_hash(ivi_hash* hash, int key);



KHASH_MAP_INIT_INT(ivsm, s_vector*)
typedef khash_t (ivsm) ivs_hash;
void append_to_ivs_hash(ivs_hash *hash, const int key, const char* val);
void destruct_ivs_hash(ivs_hash *hash);
void append_to_ivs_hash_nd(ivs_hash *hash, const int key, const char* val);
void destruct_ivs_hash_nd(ivs_hash *hash);
bool key_exist_ivs_hash(ivs_hash *hash, const int key);
s_vector* get_value_ivs_hash(ivs_hash* hash, const int key);

KHASH_MAP_INIT_STR(svim, i_vector*)
typedef khash_t (svim) svi_hash;
void append_to_svi_hash(svi_hash *hash, const char *key, int val);
void destruct_svi_hash(svi_hash *hash);
i_vector* get_value_svi_hash(svi_hash* hash, const char* key);

KHASH_MAP_INIT_STR(ssim, si_hash*)
typedef khash_t (ssim) ssi_hash;
void add_to_ssi_hash(ssi_hash *hash, const char *key1, const char *key2);
void destruct_ssi_hash(ssi_hash *hsh);
si_hash* get_value_ssi_hash(ssi_hash *hash, const char *key);
int size_content_for_key_ssi_hash(ssi_hash *hash, const char *key);


// tcr entry and container
typedef struct {
	char *va;
	char *ja;
	char *vb;
	char *jb;
	char *cdr3a;
	char *cdr3b;
	double freq;
	char *sid;
	// condition holds either tissue, timepoint, et al or combination of those
	char *condition;
} tcr_t;
tcr_t* create_tcr_entry();
void destruct_tcr_entry(tcr_t* entry);
typedef kvec_t (tcr_t *) t_vector;
void destruct_t_vector(t_vector *array);

typedef struct {
	char *v;
	char *j;
	char *cdr3;
} cdr_t;
cdr_t* create_cdr_entry();
void destruct_cdr_entry(cdr_t* entry);
typedef kvec_t (cdr_t *) c_vector;
void destruct_c_vector(c_vector *array);


typedef struct {
	char* pattern;
	float pvalue;
	float entfrc;
	s_vector *members;
}pattern_t;

pattern_t* init_pattern();
void destruct_pattern(pattern_t *p);
typedef kvec_t (pattern_t*) p_vector;
void destruct_p_vector(p_vector *array);
int pattern_compare(const void *a, const void *b);
void sort_p_vector(p_vector* array);
kstring_t* pattern_to_kstring(pattern_t *p);
pattern_t* line_to_pattern(const char* line);

typedef struct {
	char *key;
	float val;
} sf_t;
int sf_compare_by_value(const void *s1, const void *s2);
typedef kvec_t (sf_t*) sf_vector;
void destruct_sf_vector(sf_vector *array);
void sort_sf_vector(sf_vector* array);


typedef struct {
	char *key;
	int val;
} si_t;
int si_compare_by_key(const void *s1, const void *s2);

typedef kvec_t (si_t*) si_vector;
void destruct_si_vector(si_vector *array);
void sort_si_vector(si_vector* array);

typedef struct {
	int idx1;
	int idx2;
	int dist;
}i3_t;
typedef kvec_t (i3_t*) i3_vector;
int i3_compare_by_dist(const void *a, const void *b);
void destruct_i3_vector(i3_vector *array);
void sort_i3_vector(i3_vector* array);

typedef struct {
	int idx;
	int sze;
}is_t;

typedef kvec_t (is_t*) is_vector;
int is_compare_by_size(const void *a, const void *b);
void destruct_is_vector(is_vector *array);
void sort_is_vector(is_vector* array);


typedef struct _node {
	struct _node *left;
	struct _node *right;
	int idx;
} node_t;
node_t* create_node();
void destruct_node(node_t *node);

KHASH_MAP_INIT_INT(itm, node_t*)
typedef khash_t (itm) it_hash;
void append_to_it_hash(it_hash *hash, int key, node_t* val);
node_t* get_value_it_hash(it_hash *hash, int key);
bool key_exist_it_hash(it_hash *hash, int key);
void destruct_it_hash(it_hash *hash);

void find_all_leaf_index(node_t* root, i_vector *vec);

typedef struct {
	char *pattern;
	int target_count;
	int target_total;
	int refer_count;
	int refer_total;
	int ove;
	double p_value;
	int number_dot;
	float entropy_fraction;
} all_t;
int all_compar(const void *a, const void *b);

typedef kvec_t (all_t*) all_vector;
void destruct_all_vector(all_vector *array);
void sort_all_vector(all_vector* array);
KHASH_MAP_INIT_STR(avm, all_vector*)
typedef khash_t (avm) av_hash;
void append_to_av_hash(av_hash* hash, char* key, all_t *pat);
void destruct_av_hash(av_hash *hash);

#endif
