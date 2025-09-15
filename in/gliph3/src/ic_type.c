#include "ic_type.h"
#include "ic_util.h"

void destruct_s_vector(s_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		safe_free(kv_A(*array, i));
	}
	kv_destroy(*array);
}

void destruct_vi_vector(vi_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		i_vector *vec = kv_A(*array, i);
		destruct_i_vector(vec);
		safe_free(vec);
	}
	kv_destroy(*array);
}

void destruct_s_vector_nd(s_vector *array)
{
	kv_destroy(*array);
}

void destruct_i_vector(i_vector *array)
{
	kv_destroy(*array);
}

void destruct_iv_vector(iv_vector *array)
{
	for(int i = 0; i < kv_size(*array); i++)
	{
		i_vector *vec = kv_A(*array, i);
		destruct_i_vector(vec);
		safe_free(vec);
	}
	kv_destroy(*array);
}


void destruct_f_vector(f_vector *array)
{
	kv_destroy(*array);
}

void destruct_sf_vector(sf_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		sf_t *item = kv_A(*array, i);
		safe_free(item->key);
		safe_free(item);
	}
	kv_destroy(*array);
}

void sort_sf_vector(sf_vector* array)
{
	int length = kv_size(*array);
	sf_t **buf = icalloc(length, sf_t*);
	for (int i = 0; i < kv_size(*array); ++i)
	{
		buf[i]=kv_A(*array, i);
	}
	qsort(buf, length, sizeof(sf_t *), sf_compare_by_value);
	for(int i = 0; i < length; i++) {
		kv_A(*array, i) = buf[i];
	}
	safe_free(buf);
}

void destruct_si_vector(si_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		si_t *item = kv_A(*array, i);
		safe_free(item->key);
		safe_free(item);
	}
	kv_destroy(*array);
}

void sort_si_vector(si_vector* array)
{
	int length = kv_size(*array);
	si_t **buf = icalloc(length, si_t*);
	for (int i = 0; i < kv_size(*array); ++i)
	{
		buf[i]=kv_A(*array, i);
	}
	qsort(buf, length, sizeof(si_t *), si_compare_by_key);
	for(int i = 0; i < length; i++) {
		kv_A(*array, i) = buf[i];
	}
	safe_free(buf);
}

void append_to_s_set_nd(s_set *hash, const char* key)
{
	int ret;
	khiter_t k;

	k = kh_get(set, hash, key);
	if (k == kh_end(hash))
	{
		k = kh_put(set, hash, key, &ret);
		if(!ret)
		{
			kh_del(set, hash, k);
		}
	}
}

void destruct_s_set_nd(s_set *sset)
{
	kh_destroy(set, sset);
}

void append_to_s_set(s_set *sset, const char* key)
{
	int ret;
	khiter_t k;

	k = kh_get(set, sset, key);
	if (k == kh_end(sset))
	{
		char* ckey = strdup(key);
		k = kh_put(set, sset, ckey, &ret);
		if(!ret)
		{
			safe_free(ckey);
			kh_del(set, sset, k);
		}
	}
}

void destruct_s_set(s_set *sset)
{
	for (khiter_t k = kh_begin(sset); k != kh_end(sset); k++)
	{
		if (kh_exist(sset, k))
		{
			char* key = (char*)kh_key(sset, k);
			safe_free(key);
		}
	}
	kh_destroy(set, sset);
}

bool key_exist_s_set(s_set *hash, const char *key)
{
	khiter_t k = kh_get(set, hash, key);
	if(k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

char* get_value_s_set(s_set *hash, const char *key)
{
	khiter_t k = kh_get(set, hash, key);
	if(k != kh_end(hash))
	{
		char* key =(char*)kh_key(hash, k);
		return key;
	}
	else
	{
		return(0);
	}
}

void destruct_si_hash(si_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			safe_free(kh_key(hash, k));
		}
	}
	kh_destroy(sim, hash);
}

void destruct_si_hash_nd(si_hash *hash)
{
	kh_destroy(sim, hash);
}

void add_to_si_hash(si_hash *hash, char* key)
{
	int ret;
	khiter_t k;

	k = kh_get(sim, hash, key);
	if (k == kh_end(hash))
	{
		char* ckey = strdup(key);
		k = kh_put(sim, hash, ckey, &ret);
		if(!ret)
		{
			safe_free(ckey);
			kh_del(sim, hash, k);
		}
		kh_val(hash, k) = 1;
	}
	else
	{
		kh_val(hash, k) += 1;
	}
}

void add_to_si_hash_nd(si_hash *hash, char* key)
{
	int ret;
	khiter_t k;

	k = kh_get(sim, hash, key);
	if (k == kh_end(hash))
	{
		k = kh_put(sim, hash, key, &ret);
		kh_val(hash, k) = 1;
	}
	else
	{
		kh_val(hash, k) += 1;
	}
}

void append_to_si_hash(si_hash *hash, const char *key, int val)
{
	int ret;
	khiter_t k;

	k = kh_get(sim, hash, key);
	if (k == kh_end(hash))
	{
		char* ckey = strdup(key);
		k = kh_put(sim, hash, ckey, &ret);
		if(!ret)
		{
			safe_free(ckey);
			kh_del(sim, hash, k);
			return;
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) = val;
	}
}

void append_to_si_hash_nd(si_hash *hash, const char *key, int val)
{
	int ret;
	khiter_t k;

	k = kh_get(sim, hash, key);
	if (k == kh_end(hash))
	{
		k = kh_put(sim, hash, key, &ret);
		if(!ret)
		{
			kh_del(sim, hash, k);
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) = val;
	}
}

void append_to_sf_hash(sf_hash *hash, const char *key, float val)
{
	int ret;
	khiter_t k;

	k = kh_get(sfm, hash, key);
	if (k == kh_end(hash))
	{
		char* ckey = strdup(key);
		k = kh_put(sfm, hash, ckey, &ret);
		if(!ret)
		{
			safe_free(ckey);
			kh_del(sfm, hash, k);
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) += val;
	}
}

float get_value_sf_hash(sf_hash *hash, const char *key, float default_value)
{
	khiter_t k = kh_get(sfm, hash, key);
	if (k == kh_end(hash))
	{
		return default_value;
	}
	else
	{
		return kh_val(hash, k);
	}
}

bool key_exist_si_hash(si_hash *hash, const char *key)
{
	khiter_t k = kh_get(sim, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

int get_value_si_hash(si_hash *hash, const char* key, const int default_value)
{
	khiter_t k = kh_get(sim, hash, key);

	if (k == kh_end(hash))
	{
		return(default_value);
	}
	else
	{
		return(kh_val(hash, k));
	}
}

i_vector* get_value_svi_hash(svi_hash* hash, const char* key)
{
	khiter_t k = kh_get(svim, hash, key);

	if (k == kh_end(hash))
	{
		return(0);
	}
	else
	{
		return(kh_val(hash, k));
	}
}

bool key_exist_sf_hash(sf_hash *hash, const char *key)
{
	khiter_t k = kh_get(sfm, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

void destruct_sf_hash(sf_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			safe_free(kh_key(hash, k));
		}
	}
	kh_destroy(sfm, hash);
}

void append_to_li_hash(li_hash *hash, unsigned long key, int val)
{
	int ret;
	khiter_t k = kh_get(lim, hash, key);

	if (k == kh_end(hash))
	{
		k = kh_put(lim, hash, key, &ret);
		if(!ret)
		{
			kh_del(lim, hash, k);
			return;
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) = val;
	}
}
void add_to_li_hash(li_hash *hash, unsigned long key, int val)
{
	int ret;
	khiter_t k = kh_get(lim, hash, key);

	if (k == kh_end(hash))
	{
		k = kh_put(lim, hash, key, &ret);
		if(!ret)
		{
			kh_del(lim, hash, k);
			return;
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) += val;
	}
}
bool key_exist_li_hash(li_hash *hash, unsigned long key)
{
	khiter_t k = kh_get(lim, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}

}
void destruct_li_hash(li_hash *hash)
{
	kh_destroy(lim, hash);
}

void append_to_ii_hash(ii_hash *hash, int key, int val)
{
	int ret;
	khiter_t k = kh_get(iim, hash, key);

	if (k == kh_end(hash))
	{
		k = kh_put(iim, hash, key, &ret);
		if(!ret)
		{
			kh_del(iim, hash, k);
			return;
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) = val;
	}
}

void add_to_ii_hash(ii_hash *hash, int key, int val)
{
	int ret;
	khiter_t k = kh_get(iim, hash, key);

	if (k == kh_end(hash))
	{
		k = kh_put(iim, hash, key, &ret);
		if(!ret)
		{
			kh_del(iim, hash, k);
			return;
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) += val;
	}
}

bool key_exist_ii_hash(ii_hash *hash, int key)
{
	khiter_t k = kh_get(iim, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

void append_to_ivi_hash(ivi_hash *hash,  int key, int val)
{
	int ret;
	khiter_t k = kh_get(ivim, hash, key);
	i_vector *vector;
	if (k == kh_end(hash))
	{
		k = kh_put(ivim, hash, key, &ret);
		if(!ret)
		{
			kh_del(ivim, hash, k);
			return;
		}
		vector = icalloc(1, i_vector);
		kv_init(*vector);
		kv_push(int, *vector, val);
		kh_val(hash, k) = vector;
	}
	else
	{
		vector = kh_val(hash, k);
		kv_push(int, *vector, val);
	}
}

void destruct_ivi_hash(ivi_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			i_vector *vector = kh_val(hash, k);
			destruct_i_vector(vector);
			safe_free(vector);
		}
	}
	kh_destroy(ivim, hash);
}
bool key_exist_ivi_hash(ivi_hash *hash, int key)
{
	khiter_t k = kh_get(ivim, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

i_vector* get_value_ivi_hash(ivi_hash* hash, int key)
{
	khiter_t k = kh_get(ivim, hash, key);

	if (k == kh_end(hash))
	{
		return(0);
	}
	else
	{
		return(kh_val(hash, k));
	}
}


bool key_exist_svs_hash(svs_hash *hash, const char *key)
{
	khiter_t k = kh_get(svsm, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

s_vector* get_value_svs_hash(svs_hash* hash, const char* key)
{
	khiter_t k = kh_get(svsm, hash, key);

	if (k == kh_end(hash))
	{
		return(0);
	}
	else
	{
		return(kh_val(hash, k));
	}
}

s_vector* get_value_ivs_hash(ivs_hash* hash, const int key)
{
	khiter_t k = kh_get(ivsm, hash, key);

	if (k == kh_end(hash))
	{
		return(0);
	}
	else
	{
		return(kh_val(hash, k));
	}
}

bool key_exist_if_hash(if_hash *hash, int key)
{
	khiter_t k = kh_get(ifm, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

void destruct_ii_hash(ii_hash *hash)
{
	kh_destroy(iim, hash);
}

void destruct_if_hash(if_hash *hash)
{
	kh_destroy(ifm, hash);
}

void append_to_svs_hash_nd(svs_hash *hash, const char *key, const char* val)
{
	int ret;
	khiter_t k = kh_get(svsm, hash, key);
	s_vector *vector;
	if (k == kh_end(hash))
	{
		k = kh_put(svsm, hash, key, &ret);
		if(!ret)
		{
			kh_del(svsm, hash, k);
			return;
		}
		vector = icalloc(1, s_vector);
		kv_init(*vector);
		kv_push(char *, *vector, (char*)val);
		kh_val(hash, k) = vector;
	}
	else
	{
		vector = kh_val(hash, k);
		kv_push(char *, *vector, (char*)val);
	}
}

void destruct_svs_hash_nd(svs_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			s_vector *vector = kh_val(hash, k);
			destruct_s_vector_nd(vector);
			safe_free(vector);
		}
	}
	kh_destroy(svsm, hash);
}

void append_to_svs_hash(svs_hash *hash, const char *key, const char* val)
{
	int ret;
	khiter_t k = kh_get(svsm, hash, key);
	s_vector *vector;
	char *ckey;
	if (k == kh_end(hash))
	{
		ckey = strdup(key);
		k = kh_put(svsm, hash, ckey, &ret);
		if(!ret)
		{
			safe_free(ckey);
			kh_del(svsm, hash, k);
			return;
		}
		vector = icalloc(1, s_vector);
		kv_init(*vector);
		kv_push(char *, *vector, strdup(val));
		kh_val(hash, k) = vector;
	}
	else
	{
		vector = kh_val(hash, k);
		kv_push(char *, *vector, strdup(val));
	}
}

void destruct_svs_hash(svs_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			s_vector *vector = kh_val(hash, k);
			destruct_s_vector(vector);
			safe_free(vector);
			char *key = (char *)kh_key(hash, k);
			safe_free(key);
		}
	}
	kh_destroy(svsm, hash);
}

void append_to_ivs_hash(ivs_hash *hash, const int key, const char* val)
{
	int ret;
	khiter_t k = kh_get(ivsm, hash, key);
	s_vector *vector;
	if (k == kh_end(hash))
	{
		k = kh_put(ivsm, hash, key, &ret);
		if(!ret)
		{
			kh_del(ivsm, hash, k);
			return;
		}
		vector = icalloc(1, s_vector);
		kv_init(*vector);
		kv_push(char *, *vector, strdup(val));
		kh_val(hash, k) = vector;
	}
	else
	{
		vector = kh_val(hash, k);
		kv_push(char *, *vector, strdup(val));
	}
}

void destruct_ivs_hash(ivs_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			s_vector *vector = kh_val(hash, k);
			destruct_s_vector(vector);
			safe_free(vector);
		}
	}
	kh_destroy(ivsm, hash);
}

//not to duplicate insert, just store pointer
void append_to_ivs_hash_nd(ivs_hash *hash, const int key, const char* val)
{
	int ret;
	khiter_t k = kh_get(ivsm, hash, key);
	s_vector *vector;
	if (k == kh_end(hash))
	{
		k = kh_put(ivsm, hash, key, &ret);
		if(!ret)
		{
			kh_del(ivsm, hash, k);
			return;
		}
		vector = icalloc(1, s_vector);
		kv_init(*vector);
		kv_push(char *, *vector, (char*)val);
		kh_val(hash, k) = vector;
	}
	else
	{
		vector = kh_val(hash, k);
		kv_push(char *, *vector, (char*)val);
	}
}

//not destruct the content, just destroy container
void destruct_ivs_hash_nd(ivs_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			s_vector *vector = kh_val(hash, k);
			kv_destroy(*vector);
			safe_free(vector);
		}
	}
	kh_destroy(ivsm, hash);
}

bool key_exist_ivs_hash(ivs_hash *hash, const int key)
{
	khiter_t k = kh_get(ivsm, hash, key);

	if (k == kh_end(hash))
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

void append_to_svi_hash(svi_hash *hash, const char *key, int val)
{
	int ret;
	khiter_t k = kh_get(svim, hash, key);
	i_vector *vector;
	char *ckey;
	if (k == kh_end(hash))
	{
		ckey = strdup(key);
		k = kh_put(svim, hash, ckey, &ret);
		if(!ret)
		{
			safe_free(ckey);
			kh_del(svim, hash, k);
			return;
		}
		vector = icalloc(1, i_vector);
		kv_init(*vector);
		kv_push(int, *vector, val);
		kh_val(hash, k) = vector;
	}
	else
	{
		vector = kh_val(hash, k);
		kv_push(int, *vector, val);
	}
}

void destruct_svi_hash(svi_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			i_vector *vector = kh_val(hash, k);
			destruct_i_vector(vector);
			safe_free(vector);
			char *key = (char *)kh_key(hash, k);
			safe_free(key);
		}
	}
	kh_destroy(svim, hash);
}

tcr_t* create_tcr_entry()
{
	tcr_t *entry = icalloc(1, tcr_t);
	entry->va = NULL;
	entry->ja = NULL;
	entry->vb = NULL;
	entry->jb = NULL;
	entry->cdr3a = NULL;
	entry->cdr3b = NULL;
	entry->freq = 0;
	entry->sid = NULL;
	entry->condition = NULL;
	return entry;
}

void destruct_tcr_entry(tcr_t* entry)
{
	if (entry == NULL)
	{
		return;
	}

	if(entry->va != NULL)
	{
		safe_free(entry->va);
	}

	if(entry->ja != NULL)
	{
		safe_free(entry->ja);
	}

	if(entry->vb != NULL)
	{
		safe_free(entry->vb);
	}

	if(entry->jb != NULL)
	{
		safe_free(entry->jb);
	}

	if(entry->cdr3a != NULL)
	{
		safe_free(entry->cdr3a);
	}

	if(entry->cdr3b != NULL)
	{
		safe_free(entry->cdr3b);
	}

	if(entry->sid != NULL)
	{
		safe_free(entry->sid);
	}

	if(entry->condition != NULL)
	{
		safe_free(entry->condition);
	}

	safe_free(entry);
}

void destruct_t_vector(t_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		tcr_t *entry = (tcr_t*)kv_A(*array, i);
		destruct_tcr_entry(entry);
	}
	kv_destroy(*array);
}

cdr_t* create_cdr_entry()
{
	cdr_t *entry = icalloc(1, cdr_t);
	entry->v = NULL;
	entry->j = NULL;
	entry->cdr3 = NULL;
	return entry;
}

void destruct_cdr_entry(cdr_t* entry)
{
	if (entry == NULL)
	{
		return;
	}

	if(entry->v != NULL)
	{
		safe_free(entry->v);
	}

	if(entry->j != NULL)
	{
		safe_free(entry->j);
	}

	if(entry->cdr3 != NULL)
	{
		safe_free(entry->cdr3);
	}

	safe_free(entry);
}

void destruct_c_vector(c_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		cdr_t *entry = (cdr_t*)kv_A(*array, i);
		destruct_cdr_entry(entry);
	}
	kv_destroy(*array);
}

int sf_compare_by_value(const void *a, const void *b)
{
	sf_t *oa = *(sf_t *const *)a;
	sf_t *ob = *(sf_t *const *)b;

	if (oa->val > ob->val)
	{
		return(1);
	}
	else
	{
		return(-1);
	}
}

int si_compare_by_key(const void *a, const void *b)
{
	si_t *oa = *(si_t *const *)a;
	si_t *ob = *(si_t *const *)b;

	return(strcmp(oa->key, ob->key));
}

int is_compare_by_size(const void *a, const void *b)
{
	is_t *oa = *(is_t *const *)a;
	is_t *ob = *(is_t *const *)b;

	return(ob->sze - oa->sze);
}

void destruct_is_vector(is_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		is_t *obj = (is_t*)kv_A(*array, i);
		safe_free(obj);
	}
	kv_destroy(*array);
}

void sort_is_vector(is_vector* array)
{
	int length = kv_size(*array);
	is_t **buf = icalloc(length, is_t*);
	for (int i = 0; i < kv_size(*array); ++i)
	{
		buf[i]=kv_A(*array, i);
	}
	qsort(buf, length, sizeof(is_t *), is_compare_by_size);
	for(int i = 0; i < length; i++) {
		kv_A(*array, i) = buf[i];
	}
	safe_free(buf);
}

int i3_compare_by_dist(const void *a, const void *b)
{
	i3_t *oa = *(i3_t *const *)a;
	i3_t *ob = *(i3_t *const *)b;

	return(ob->dist - oa->dist);
}

void destruct_i3_vector(i3_vector *array)
{
	for (int i = 0; i < kv_size(*array); ++i)
	{
		i3_t *obj = (i3_t*)kv_A(*array, i);
		safe_free(obj);
	}
	kv_destroy(*array);
}

void sort_i3_vector(i3_vector* array)
{
	int length = kv_size(*array);
	i3_t **buf = icalloc(length, i3_t*);
	for (int i = 0; i < kv_size(*array); ++i)
	{
		buf[i]=kv_A(*array, i);
	}
	qsort(buf, length, sizeof(i3_t *), i3_compare_by_dist);
	for(int i = 0; i < length; i++) {
		kv_A(*array, i) = buf[i];
	}
	safe_free(buf);
}

void append_to_it_hash(it_hash *hash, int key, node_t* val)
{
	int ret;
	khiter_t k = kh_get(itm, hash, key);
	if (k == kh_end(hash))
	{
		k = kh_put(itm, hash, key, &ret);
		if(!ret)
		{
			kh_del(itm, hash, k);
			return;
		}
		kh_val(hash, k) = val;
	}
	else
	{
		kh_val(hash, k) = val;
	}
}

bool key_exist_it_hash(it_hash *hash, int key)
{
	khiter_t k = kh_get(itm, hash, key);
	if (k == kh_end(hash))
	{
		return(false);
	}
	else{
		return(true);
	}
}

node_t* get_value_it_hash(it_hash *hash, int key){
	khiter_t k = kh_get(itm, hash, key);
	if (k == kh_end(hash))
	{
		return(NULL);
	}
	else{
		node_t* node = kh_val(hash, k);
		return node;
	}
}

void destruct_it_hash(it_hash *hash)
{
	kh_destroy(itm, hash);
}

node_t* create_node()
{
	node_t* node = icalloc(1, node_t);
	node->left = 0;
	node->right = 0;
	node->idx = -1;
	return node;
}

void destruct_node(node_t *node)
{
	if(node == 0)
	{
		return;
	}

	if(node->left != 0)
	{
		destruct_node(node->left);
	}

	if(node->right != 0)
	{
		destruct_node(node->right);
	}

	safe_free(node);
}

void find_all_leaf_index(node_t* root, i_vector *vec)
{
	if(root == 0)
	{
		return;
	}

	if(root->left != 0)
	{
		find_all_leaf_index(root->left, vec);
	}
	if(root->idx >= 0)
	{
		kv_push(int, *vec, root->idx);
	}
	if(root->right != 0)
	{
		find_all_leaf_index(root->right, vec);
	}
}

pattern_t* init_pattern()
{
	pattern_t *p = icalloc(1, pattern_t);
	p->pattern = 0;
	p->pvalue = 0.0;
	p->members = icalloc(1, s_vector);
	kv_init(*(p->members));
	return p;
}

void destruct_pattern(pattern_t *p)
{
	// debug("0.3");
	destruct_s_vector(p->members);
	safe_free(p->members);
	// debug("0.4");
	if (p->pattern != 0)
	{
		safe_free(p->pattern);
	}
	// debug("0.5");
	safe_free(p);
}

void destruct_p_vector(p_vector *array)
{
	// debug("0.1");
	for(int i = 0; i < kv_size(*array); i++)
	{
		pattern_t *p = (pattern_t*)kv_A(*array, i);
		destruct_pattern(p);
	}
	// debug("0.2");
	kv_destroy(*array);
}

void sort_p_vector(p_vector* array)
{
	int length = kv_size(*array);
	pattern_t **buf = icalloc(length, pattern_t*);
	for (int i = 0; i < kv_size(*array); ++i)
	{
		buf[i]=kv_A(*array, i);
	}
	qsort(buf, length, sizeof(pattern_t *), pattern_compare);
	for(int i = 0; i < length; i++) {
		kv_A(*array, i) = buf[i];
	}
	safe_free(buf);
}

int pattern_compare(const void *a, const void *b)
{
	pattern_t *oa = *(pattern_t *const *)a;
	pattern_t *ob = *(pattern_t *const *)b;

	if(fabsf(oa->entfrc - ob->entfrc) < EPSILON)
	{
		if(fabsf(oa->pvalue - ob->pvalue) < EPSILON)
		{
			return(0);
		}else if (oa->pvalue > ob->pvalue)
		{
			return(1);
		}else{
			return(-1);
		}
	}
	else if (oa->entfrc > ob->entfrc)
	{
		return(1);
	}
	else
	{
		return(-1);
	}
}

kstring_t* pattern_to_kstring(pattern_t *p)
{
	assert(p != 0);
	kstring_t *kstr  = icalloc(1, kstring_t);
	ksprintf(kstr, "%s %.0e %.2f", p->pattern, p->pvalue, p->entfrc);
	for(int i = 0; i < kv_size(*(p->members)); i++)
	{
		ksprintf(kstr, " %s", (char*)kv_A(*(p->members), i));
	}
	ksprintf(kstr, "\n");
	return kstr;
}

pattern_t* line_to_pattern(const char* line)
{
	s_vector fields;
	kv_init(fields);
	ic_split_string(line, ' ', &fields);
	pattern_t *p = init_pattern();
	p->pattern = strdup((char*)kv_A(fields, 0));
	p->pvalue = atof((char*)kv_A(fields, 1));
	p->entfrc = atof((char*)kv_A(fields, 2));
	for(int i = 3; i < kv_size(fields); i++)
	{
		kv_push(char*, *(p->members), strdup((char*)kv_A(fields, i)));
	}
	destruct_s_vector(&fields);
	return p;
}

void add_to_ssi_hash(ssi_hash *hash, const char *key1, const char *key2)
{
	khiter_t k = kh_get(ssim, hash, key1);
	int ret;

	if (k == kh_end(hash))
	{
		k = kh_put(ssim, hash, strdup(key1), &ret);
		si_hash *hsh = kh_init(sim);
		add_to_si_hash(hsh, key2);
		kh_val(hash, k) = hsh;
	}
	else
	{
		si_hash *hsh = kh_val(hash, k);
		add_to_si_hash(hsh, key2);
	}
}

int size_content_for_key_ssi_hash(ssi_hash *hash, const char*key)
{
	khiter_t k = kh_get(ssim, hash, key);

	if (k == kh_end(hash))
	{
		return 0;
	}
	else
	{
		si_hash *hsh = kh_val(hash, k);
		return kh_size(hsh);
	}
}

void destruct_ssi_hash(ssi_hash *hsh)
{
	for(khiter_t k = kh_begin(hsh); k != kh_end(hsh); k++)
	{
		if(kh_exist(hsh, k))
		{
			char* key = (char*) kh_key(hsh, k);
			si_hash *hash = kh_val(hsh, k);
			safe_free(key);
			destruct_si_hash(hash);
		}
	}
	kh_destroy(ssim, hsh);
}

si_hash* get_value_ssi_hash(ssi_hash *hash, const char *key)
{
	khiter_t k = kh_get(ssim, hash, key);

	if (k == kh_end(hash))
	{
		return 0;
	}
	else
	{
		si_hash *hsh = kh_val(hash, k);
		return hsh;
	}
}

void sort_i_vector(i_vector *array)
{
	int size  = kv_size(*array);
	int *buf = icalloc(size, int);

	for (int i = 0; i < size; ++i)
	{
		buf[i] = kv_pop(*array);
	}
	qsort(buf, size, sizeof(int), int_compar);
	for (int i = 0; i < size; ++i)
	{
		kv_push(int, *array, buf[i]);
	}
	safe_free(buf);
}
int int_compar(const void* a, const void* b)
{
	int aa = *(const int *)a;
	int bb = *(const int *)b;
	return (aa-bb);
}

void sort_s_vector_by_length(s_vector *array)
{
	int size  = kv_size(*array);
	char **buf = icalloc(size, char*);

	for (int i = 0; i < size; ++i)
	{
		buf[i] = kv_pop(*array);
	}
	qsort(buf, size, sizeof(char *), str_length_compar);
	for (int i = 0; i < size; ++i)
	{
		kv_push(char*, *array, buf[i]);
	}
	safe_free(buf);
}
int str_length_compar(const void* a, const void* b)
{
	const char* aa = *(const char**)a;
	const char* bb = *(const char**)b;
	return (strlen(aa) - strlen(bb));
}


int str_compar(const void* a, const void* b)
{
	const char* aa = *(const char**)a;
	const char* bb = *(const char**)b;
	return strcmp(aa,bb);
}

void sort_s_vector(s_vector *array)
{
	int size  = kv_size(*array);
	char **buf = icalloc(size, char*);

	for (int i = 0; i < size; ++i)
	{
		buf[i] = kv_pop(*array);
	}
	qsort(buf, size, sizeof(char *), str_compar);
	for (int i = 0; i < size; ++i)
	{
		kv_push(char*, *array, buf[i]);
	}
	safe_free(buf);
}

int all_compar(const void *a, const void *b)
{
	const all_t *aa = *(all_t* const *)a;
	const all_t *bb = *(all_t* const *)b;

	if(aa->entropy_fraction * 1.05 < bb->entropy_fraction)
	{
		return -1;
	}else if (aa->entropy_fraction > 1.05 * bb->entropy_fraction)
	{
		return 1;
	}else{
		if(aa->p_value * 1.05 < bb->p_value)
		{
			return -1;
		}
		else if(aa->p_value > 1.05 * bb->p_value)
		{
			return 1;
		}
		else{
			// pvalue is the same
			if (aa->ove > bb->ove)
			{
				return -1;
			}
			else if (aa->ove < bb->ove)
			{
				return 1;
			}else{
				// p_value and ove are the same
				if(aa->number_dot < bb->number_dot)
				{
					return -1;
				}else if (aa->number_dot > bb->number_dot)
				{
					return 1;
				}else{
					return 0;
				}
			}
		}
	}
}

void destruct_all_vector(all_vector *array)
{
	for(int i = 0; i < kv_size(*array); i++)
	{
		all_t *entry = kv_A(*array, i);
		safe_free(entry->pattern);
		safe_free(entry);
	}
	kv_destroy(*array);
}

void sort_all_vector(all_vector* array)
{
	int size  = kv_size(*array);
	all_t **buf = icalloc(size, all_t*);

	for (int i = 0; i < size; ++i)
	{
		buf[i] = kv_pop(*array);
	}
	qsort(buf, size, sizeof(all_t *), all_compar);
	for (int i = 0; i < size; ++i)
	{
		kv_push(all_t*, *array, buf[i]);
	}
	safe_free(buf);
}

void append_to_av_hash(av_hash* hash, char* key, all_t *pat)
{
	int ret;
	khiter_t k = kh_get(avm, hash, key);
	all_vector *vector;
	char *ckey;
	if (k == kh_end(hash))
	{
		ckey = strdup(key);
		k = kh_put(avm, hash, ckey, &ret);
		if(!ret)
		{
			safe_free(ckey);
			kh_del(avm, hash, k);
			return;
		}
		vector = icalloc(1, all_vector);
		kv_init(*vector);
		kv_push(all_t*, *vector, pat);
		kh_val(hash, k) = vector;
	}
	else
	{
		vector = kh_val(hash, k);
		kv_push(all_t*, *vector, pat);
	}
}

void destruct_av_hash(av_hash *hash)
{
	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			all_vector *vector = kh_val(hash, k);
			destruct_all_vector(vector);
			safe_free(vector);
			char *key = (char *)kh_key(hash, k);
			safe_free(key);
		}
	}
	kh_destroy(avm, hash);
}
