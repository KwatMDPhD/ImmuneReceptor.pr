#include "ic_util.h"
#include "ic_clust.h"
#include "heap.h"
#include "queue.h"

void safer_free(void **pp)
{
	if (pp != NULL && *pp != NULL)
	{
		free(*pp);
		*pp = NULL;
	}
}

inline void *_icalloc(long number, int size)
{
	void *p = (void *)calloc(number, size);
	if (p == NULL)
	{
		die("icalloc failure requesting %d of size %d bytes", number, size);
	}
	return(p);
}

void die(const char* format, ...)
{
	va_list vargs;
	va_start (vargs, format);
	fprintf (stderr, "%d: ", __LINE__);
	vfprintf (stderr, format, vargs);
	fprintf (stderr, ".\n");
	va_end (vargs);
	exit (EXIT_FAILURE);
}

void debug(char *format, ...)
{
	time_t T  = time(NULL);
	struct tm tm = *localtime(&T);
	va_list args;

	va_start(args, format);
	fprintf(stderr, "DEBUG @ %02d:%02d:%02d: ", tm.tm_hour, tm.tm_min, tm.tm_sec);
	vfprintf(stderr, format, args);
	fprintf(stderr, "\n");
	va_end(args);
}

void system_call(char *cmd)
{
	// debug(cmd);
	if (system(cmd) != 0)
	{
		die(cmd);
	}
}

char* ic_generate_random_string(int size)
{
	time_t t;
	srand((unsigned) time(&t));
	const char orign[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
	char* rstr = icalloc(size+1, char);
	for (int i = 0; i < size; i++) {
		int idx = rand() % (int) (sizeof(orign) - 1);
		rstr[i] = orign[idx];
	}
	rstr[size] = '\0';
	return rstr;
}


//compute fisher exact test
double gammln(double xx)
{
	static double cof[6] = { 76.18009172947146, -86.50532032941677,
		                 24.01409824083091, -1.231739572450155,
		                 0.1208650973866179e-2, -0.5395239384953e-5 };
	double x, tmp, ser;
	int j;

	x    = xx - 1.0;
	tmp  = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser  = 1.0;
	for (j = 0; j <= 5; j++)
	{
		x   += 1.0;
		ser += cof[j] / x;
	}
	return(-tmp + log(2.50662827465 * ser));
}

double factln(int n)
{
	static double a[10001];

	if (n <= 1)
	{
		return(0.0);
	}
	if (n <= 10000)
	{
		return(a[n] ? a[n] : (a[n] = gammln((double)(n + 1.0))));
	}
	else
	{
		return(gammln((double)(n + 1.0)));
	}
}
/*
   https://towardsdatascience.com/fishers-exact-test-from-scratch-with-python-2b907f29e593

   Fisher’s exact test is used to determine whether there is a significant association between two categorical variables in a contingency table.
   Fisher’s exact test is an alternative to Pearson’s chi-squared test for independence. While actually valid for all sample sizes, Fisher’s exact test is practically applied when sample sizes are small.
   A general recommendation is to use Fisher’s exact test- instead of the chi-squared test - whenever more than 20 % of cells in a contingency table have expected frequencies < 5

   In contrast to other statistical tests (e.g. z-test, chi-squared test), where p-values are calculated using an approximation to the true distribution (e.g. normal distribution, χ2-distribution) Fisher’s exact test yields exact p-values, hence the name. That’s why some authors prefer Fisher’s exact test over Pearson’s chi-squared test whenever (computationally) possible — even with larger sample sizes.

   Given a 2x2 table:
 +---+---+
 | a | b |
 +---+---+
 | c | d |
 +---+---+
 */
double fisher_exact_test(int n11, int n12, int n21, int n22)
{
	double p_value = 0.0;
	int a, b, c, d;
	double x = 0.0;

	for (a = n11; a <= n11 + n12; a++)
	{
		b = n11 + n12 - a;
		c = n11 + n21 - a;
		d = n12 + n22 - b;
		if (a < 0 || b < 0 || c < 0 || d < 0)
		{
			break;
		}
		x = factln(a + b) + factln(c + d) + factln(a + c) + factln(b + d) - factln(a) - factln(b) - factln(c) - factln(d) - factln(a + b + c + d);
		p_value += exp(x);
	}
	return(p_value);
}

double
compute_enrichment_pvalue(int population, int sub_population_with_feature, int sample_population, int sample_population_with_feature)
{
	int k = sample_population; // number of draws
	int x = sample_population_with_feature;  // number of observed success
	int m = sub_population_with_feature; // number of successes in the population
	int n = population - m; // number of unsuccesses in the population
	double pvalue = 0.0;

	for (int i = x; i <= k && i <= m; ++i)
	{
		pvalue += hypergeometric_pmf(i, m, n, k);
	}
	return(pvalue);
}

/*
 * Given a population consisting of `m` items of class M and `n` items of
 * class N, this returns the probability of observing `x` items of class M
 * when sampling `k` times without replacement from the entire population
 * (i.e., {M,N}) p(x) = (choose(m, x) * choose(n, k-x)) / choose(m+n, k)
 */
double hypergeometric_pmf(int x, int m, int n, int k)
{
	double a = log_binomial_coefficient(m, x);
	double b = log_binomial_coefficient(n, k - x);
	double c = log_binomial_coefficient(m + n, k);

	return(exp(a + b - c));
}

double log_binomial_coefficient(int population, int sample)
{
	int s = sample > population - sample ? sample : population - sample;

	assert(s <= population);
	assert(population >= 0);
	if (s == population)
	{
		return(0.0);
	}
	double ret = 0.0;
	for (int i = s + 1; i <= population; ++i)
	{
		ret += log(i) - log(i - s);
	}
	return(ret);
}

// string manuplication
void ic_to_upper(char s[])
{
	int i = 0;

	while (s[i] != '\0')
	{
		if (s[i] >= 'a' && s[i] <= 'z')
		{
			s[i] = s[i] - 32;
		}
		i++;
	}
}

void ic_to_lower(char s[])
{
	int i = 0;

	while (s[i] != '\0')
	{
		if (s[i] >= 'A' && s[i] <= 'Z')
		{
			s[i] = s[i] + 32;
		}
		i++;
	}
}

bool ic_starts_with(const char *str, const char *pattern)
{
	if ((strstr(str, pattern) - str) == 0)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

bool ic_ends_with(const char *str, const char *pattern)
{
	const char *init = pattern;       /* Hold the initial position of *pattern */

	while (*str)
	{
		while (*str == *pattern)
		{
			if (!(*str))
			{
				return(true);
			}
			str++;
			pattern++;
		}
		str++;
		pattern = init;
	}
	return(false);
}

void ic_split_string(const char* istring, char delim, s_vector *array)
{
	kstring_t *s = 0;
	int *fields, n = 0;
	s = icalloc(1, kstring_t);
	ksprintf(s, "%s", istring);
	fields = ksplit(s, delim, &n);
	for (int i = 0; i < n; i++)
	{
		kv_push(char *, *array, strdup(s->s + fields[i]));
	}
	safe_free(fields);
	safe_free(s->s);
	safe_free(s);
}

bool ic_all_space(const char *str)
{
	while (*str)
	{
		if (!isspace(*str++))
		{
			return(false);
		}
	}
	return(true);
}

char* ic_strmove(char *tgt, char *src)
{
	return(memmove(tgt, src, strlen(src) + 1));
}

// remove leading whitespaces from a string
char* ic_strtriml(char *str)
{
	char *p = str;
	while (*p && isspace((int)*p))
	{
		p++;
	}
	if (p - str)
	{
		ic_strmove(str, p);
	}
	return(str);
}

// remove trailing whitespaces from a string
char* ic_strtrimr(char *str)
{
	// The C library function char *strchr(const char *str, int c) searches for
	// the first occurrence of the character c (an unsigned char) in the string
	// pointed to by the argument str.
	char *p = strchr(str, 0);
	while ((p - 1) - str && isspace((int)*(p - 1)))
	{
		p--;
	}
	*p = 0;

	return(str);
}

char* ic_strip(char* str)
{
	str = ic_strtriml(str);
	str = ic_strtrimr(str);
	return(str);
}

void read_refer_cdr3(const char* ifile, c_vector *vector)
{
	FILE *fp = fopen(ifile, "r");

	if (fp == NULL)
	{
		fprintf(stderr, "Could not open file:%s\n", ifile);
		exit(EXIT_FAILURE);
	}

	char *line;
	s_vector fields;

	while ((line = ic_readline(fp)) != NULL)
	{
		if (ic_all_space(line) == 1)
		{
			continue;
		}
		ic_strip(line);
		if (line[0] == '#')
		{
			continue;
		}
		kv_init(fields);
		ic_split_string(line, ',', &fields);
		assert(kv_size(fields) > 2);
		cdr_t *cdr = create_cdr_entry();
		cdr->cdr3 = strdup((char*) kv_A(fields, 0));
		cdr->v = strdup((char*) kv_A(fields, 1));
		cdr->j = strdup((char*) kv_A(fields, 2));
		kv_push(cdr_t*, *vector, cdr);
		destruct_s_vector(&fields);
		safe_free(line);
	}
	fclose(fp);
}


void read_hla(parameter_t *global_options)
{
	global_options->hla_hash = kh_init(svsm);
	t_vector *tvec = global_options->target_cdr3_vector;
	si_hash *sid_dict = kh_init(sim);
	for(int i = 0; i < kv_size(*tvec); i++)
	{
		tcr_t *ent = kv_A(*tvec, i);
		add_to_si_hash(sid_dict, ent->sid);
	}
	s_vector larray, fields, hla_fields;
	kv_init(larray);
	ic_read_into_line_vector(global_options->hla_file, &larray);
	char append_field[] = "01";
	char hla[MAX_HLA_LENGTH];
	int pos = 0;
	for (int lnum= 0; lnum< kv_size(larray); ++lnum)
	{
		char* cline = (char*)kv_A(larray, lnum);
		cline = ic_strtriml(cline);
		cline = ic_strtrimr(cline);
		if (cline[0] == '#')
		{
			continue;
		}
		kv_init(fields);
		if(strchr(cline, '\t') != NULL)
		{
			ic_split_string(cline, '\t', &fields);
		}else{
			ic_split_string(cline, ',', &fields);
		}

		char* sid = (char*)kv_A(fields, 0);
		if(key_exist_si_hash(sid_dict, sid) == false)
		{
			destruct_s_vector(&fields);
			continue;
		}
		for (int i = 1; i < kv_size(fields); ++i)
		{
			char *allele = (char*)kv_A(fields, i);
			ic_to_upper(allele);

			if(strchr(allele, '*') != NULL)
			{
				kv_init(hla_fields);
				ic_split_string(allele, ':', &hla_fields);
				for(int j = kv_size(hla_fields); j < global_options->number_of_hla_field; j++)
				{
					kv_push(char*, hla_fields, strdup(append_field));
				}
				pos = 0;
				pos += snprintf(hla + pos, MAX_HLA_LENGTH - pos, "%s", kv_A(hla_fields, 0));
				for(int j = 1; j < kv_size(hla_fields) && j < global_options->number_of_hla_field; j++)
				{
					pos += snprintf(hla + pos, MAX_HLA_LENGTH - pos, ":%s", kv_A(hla_fields, j));
				}
				destruct_s_vector(&hla_fields);
				append_to_svs_hash(global_options->hla_hash, sid, hla);
			}
		}
		destruct_s_vector(&fields);
	}
	destruct_s_vector(&larray);
	destruct_si_hash(sid_dict);
}

char* purge_cdr3(parameter_t *global_options, const char* chain, const char* cdr3)
{
	bool good = true;
	size_t len = strlen(cdr3);
	char new_cdr3[MAX_CDR3_LENGTH];

	if(strcmp(chain, CHAIN_AG) == 0)
	{
		if(len < global_options->ag_cdr3_length_min_cutoff || len > global_options->ag_cdr3_length_max_cutoff)
		{
			good = false;
		}
	}
	else
	{
		if(len < global_options->bd_cdr3_length_min_cutoff || len > global_options->bd_cdr3_length_max_cutoff)
		{
			good = false;
		}
	}

	for(int i = 0; i < strlen(cdr3); i++)
	{
		char aa = cdr3[i];
		if (strchr(AMINOACIDS, aa) == NULL)
		{
			good = false;
			break;
		}
	}

	if (good == true)
	{
		return strdup(cdr3);
	}
	else
	{
		sprintf(new_cdr3, NA_SYMBOL);
		return strdup(new_cdr3);
	}
}

void format_input_cdr3_file(parameter_t *global_options)
{
	FILE *fp = fopen(global_options->cdr3_file, "r");

	if (fp == NULL)
	{
		fprintf(stderr, "Could not open file:%s\n", global_options->cdr3_file);
		exit(EXIT_FAILURE);
	}

	char *line;
	bool flag = true;
	bool sid_is_found = false;
	s_vector fields, headers;
	tcr_t *entry;
	char na[] = "NA";
	sf_hash *unq = kh_init(sfm);
	char key[H5];

	while ((line = ic_readline(fp)) != NULL)
	{
		if (ic_all_space(line) != 1)
		{
			ic_strip(line);
			if(flag == true)
			{
				//first line
				char *hline = line;
				ic_to_upper(hline);
				kv_init(headers);
				ic_split_string(hline, ',', &headers);
				for (int i = 0; i < kv_size(headers); ++i)
				{
					char* title = (char*)kv_A(headers, i);
					ic_strip(title);
					if(strcmp(title, "SID") == 0)
					{
						sid_is_found = true;
					}
				}
				assert(sid_is_found == true);
				flag = false;
			}else if(line[0] != '#')
			{
				kv_init(fields);
				ic_split_string(line, ',', &fields);
				if(kv_size(headers) > kv_size(fields))
				{
					debug("%s", line);
				}
				assert(kv_size(headers) <= kv_size(fields));
				entry = create_tcr_entry();
				for (int i = 0; i < kv_size(fields); ++i)
				{
					entry->freq = 1;
					char* value = (char*)kv_A(fields, i);
					char* title = (char*)kv_A(headers, i);
					value = ic_strip(value);
					if(strcmp(title, "CDR3B") == 0)
					{
						entry->cdr3b = purge_cdr3(global_options, CHAIN_BD, value);
					}

					if(strcmp(title, "CDR3A") == 0)
					{
						entry->cdr3a = purge_cdr3(global_options, CHAIN_AG, value);
					}

					if(strcmp(title, "VB") == 0)
					{
						entry->vb = strdup(value);
					}

					if(strcmp(title, "JB") == 0)
					{
						entry->jb = strdup(value);
					}

					if(strcmp(title, "VA") == 0)
					{
						entry->va = strdup(value);
					}

					if(strcmp(title, "JA") == 0)
					{
						entry->ja = strdup(value);
					}

					if(strcmp(title, "SID") == 0)
					{
						entry->sid = strdup(value);
					}

					if(strcmp(title, "CONDITION") == 0)
					{
						entry->condition = strdup(value);
					}

					if(strcmp(title, "FREQUENCY") == 0)
					{
						entry->freq = atof(value);
					}
				}
				destruct_s_vector(&fields);
				if(entry->cdr3a == 0)
				{
					entry->cdr3a = strdup(na);
				}

				if(entry->va == 0)
				{
					entry->va = strdup(na);
				}
				if(entry->ja == 0)
				{
					entry->ja = strdup(na);
				}

				if(entry->cdr3b == 0)
				{
					entry->cdr3b = strdup(na);
				}

				if(entry->vb == 0)
				{
					entry->vb = strdup(na);
				}

				if(entry->jb == 0)
				{
					entry->jb = strdup(na);
				}

				if(entry->condition == 0)
				{
					entry->condition = strdup(na);
				}

				if(entry->sid == 0)
				{
					entry->sid = strdup(na);
				}

				if(strcmp(entry->cdr3b, NA_SYMBOL) != 0 || strcmp(entry->cdr3a, NA_SYMBOL) != 0)
				{
					sprintf(key, "%s,%s,%s,%s,%s,%s,%s,%s", entry->cdr3b, entry->vb, entry->jb, entry->cdr3a, entry->va, entry->ja, entry->sid, entry->condition);
					append_to_sf_hash(unq, key, entry->freq);
				}
				destruct_tcr_entry(entry);
			}
		}
		safe_free(line);
	}
	destruct_s_vector(&headers);
	fclose(fp);
	char ofile[MAX_FILENAME_LENGTH];
	sprintf(ofile, "%s.%d", global_options->cdr3_file, getpid());
	fp = fopen(ofile, "w");

	if (fp == NULL)
	{
		fprintf(stderr, "Could not open file:%s\n", ofile);
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "CDR3B,VB,JB,CDR3A,VA,JA,SID,CONDITION,FREQUENCY\n");
	for (khiter_t k = kh_begin(unq); k != kh_end(unq); k++)
	{
		if (kh_exist(unq, k))
		{
			fprintf(fp, "%s,%.1e\n", kh_key(unq, k), kh_val(unq, k));
		}
	}
	fclose(fp);
	destruct_sf_hash(unq);
	sprintf(global_options->cdr3_file, "%s", ofile);
}

void read_target_cdr3(parameter_t *global_options, t_vector *vector)
{
	FILE *fp = fopen(global_options->cdr3_file, "r");

	if (fp == NULL)
	{
		fprintf(stderr, "Could not open file:%s\n", global_options->cdr3_file);
		exit(EXIT_FAILURE);
	}

	char *line;
	bool flag = true;
	bool sid_is_found = false;
	s_vector fields, headers;
	tcr_t *entry;

	while ((line = ic_readline(fp)) != NULL)
	{
		if (ic_all_space(line) != 1)
		{
			ic_strip(line);
			if(flag == true)
			{
				//first line
				char *hline = line;
				ic_to_upper(hline);
				kv_init(headers);
				ic_split_string(hline, ',', &headers);
				for (int i = 0; i < kv_size(headers); ++i)
				{
					char* title = (char*)kv_A(headers, i);
					ic_strip(title);
					if(strcmp(title, "SID") == 0)
					{
						sid_is_found = true;
					}
				}
				assert(sid_is_found == true);
				flag = false;
			}else if(line[0] != '#')
			{
				kv_init(fields);
				ic_split_string(line, ',', &fields);

				assert(kv_size(headers) <= kv_size(fields));
				entry = create_tcr_entry();
				for (int i = 0; i < kv_size(fields); ++i)
				{
					entry->freq = 1;
					char* value = (char*)kv_A(fields, i);
					char* title = (char*)kv_A(headers, i);
					value = ic_strip(value);
					if(strcmp(title, "CDR3B") == 0)
					{
						entry->cdr3b = purge_cdr3(global_options, CHAIN_BD, value);
					}

					if(strcmp(title, "CDR3A") == 0)
					{
						entry->cdr3a = purge_cdr3(global_options, CHAIN_AG, value);
					}

					if(strcmp(title, "VB") == 0)
					{
						entry->vb = strdup(value);
					}

					if(strcmp(title, "JB") == 0)
					{
						entry->jb = strdup(value);
					}

					if(strcmp(title, "VA") == 0)
					{
						entry->va = strdup(value);
					}

					if(strcmp(title, "JA") == 0)
					{
						entry->ja = strdup(value);
					}

					if(strcmp(title, "SID") == 0)
					{
						entry->sid = strdup(value);
					}

					if(strcmp(title, "CONDITION") == 0)
					{
						entry->condition = strdup(value);
					}

					if(strcmp(title, "FREQUENCY") == 0)
					{
						entry->freq = atof(value);
					}
				}
				destruct_s_vector(&fields);
				kv_push(tcr_t*, *vector, entry);
			}
		}
		safe_free(line);
	}
	destruct_s_vector(&headers);
	fclose(fp);
}

void read_cluster_cdr3(parameter_t *global_options, const char* ifile, const char* chain, s_set *set)
{
	FILE *cfp = fopen(ifile, "r");

	if (cfp == NULL)
	{
		fprintf(stderr, "Could not open file:%s for reading\n", ifile);
		exit(EXIT_FAILURE);
	}

	s_vector fields;
	char *line = 0, *cdr3 = 0, *vgene = 0;
	char key[H2];
	while ((line = ic_readline(cfp)) != NULL)
	{

		if (ic_all_space(line) == 1)
		{
			safe_free(line);
		}
		else
		{
			if(line[strlen(line) - 1] == '\r' || line[strlen(line) - 1] =='\n')
			{
				line[strlen(line) - 1] = 0;
			}

			kv_init(fields);
			ic_split_string(line, ',', &fields);
			for(int i =0; i<kv_size(fields); i++)
			{
				ic_strip(kv_A(fields, i));
			}

			if(strcmp(chain, CHAIN_AG) == 0)
			{
				cdr3 = kv_A(fields, 5);
				vgene = kv_A(fields, 6);
			}

			if(strcmp(chain, CHAIN_BD) == 0)
			{
				cdr3 = kv_A(fields, 8);
				vgene = kv_A(fields, 9);
			}

			if(global_options->same_v  == 1)
			{
				sprintf(key, "%s:%s", vgene, cdr3);
			}else{
				sprintf(key, "%s", cdr3);
			}

			if (strcmp(cdr3, NA_SYMBOL) != 0)
			{
				append_to_s_set(set, key);
			}
			destruct_s_vector(&fields);
			safe_free(line);
		}
	}
	fclose(cfp);
}

void output_hla_association(parameter_t *global_options, const char* chain)
{
	if(global_options->hla_availability == false)
	{
		return;
	}
	svs_hash *hla_hash = global_options->hla_hash;

	int population = kh_size(hla_hash);
	char pair[H5];
	si_hash *all_pair_hash = kh_init(sim);
	si_hash *all_allele_hash = kh_init(sim);

	for (khiter_t k = kh_begin(hla_hash); k != kh_end(hla_hash); k++)
	{
		if (kh_exist(hla_hash, k))
		{
			char *sid = (char*)kh_key(hla_hash, k);
			s_vector *allele_vector = kh_val(hla_hash, k);
			for (int a = 0; a < kv_size(*allele_vector); a++)
			{
				char *allele = kv_A(*allele_vector, a);
				sprintf(pair, "%s,%s", sid, allele);
				append_to_si_hash(all_pair_hash, pair, 0);
			}
		}
	}

	for (khiter_t k = kh_begin(all_pair_hash); k != kh_end(all_pair_hash); k++)
	{
		if (kh_exist(all_pair_hash, k))
		{
			char *sid_allele = (char*)kh_key(all_pair_hash, k);
			s_vector vec;
			kv_init(vec);
			ic_split_string(sid_allele, ',', &vec);
			char* allele = kv_A(vec, 1);
			add_to_si_hash(all_allele_hash, allele);
			destruct_s_vector(&vec);
		}
	}
	destruct_si_hash(all_pair_hash);
	char hla_out_file[MAX_FILENAME_LENGTH];
	sprintf(hla_out_file, "%s_hla.csv", global_options->out_prefix);
	FILE *hla_out_fp = NULL;
	if(access(hla_out_file, F_OK) == 0)
	{
		hla_out_fp = fopen(hla_out_file, "a");
		if (hla_out_fp == NULL)
		{
			fprintf(stderr, "Could not open file %s\n", hla_out_file);
			exit(EXIT_FAILURE);
		}
	}else{
		hla_out_fp = fopen(hla_out_file, "w");
		if (hla_out_fp == NULL)
		{
			fprintf(stderr, "Could not open file %s\n", hla_out_file);
			exit(EXIT_FAILURE);
		}
		fprintf(hla_out_fp, "chain,pattern,allele,pvalue,#subjects with HLA in total,#subjects with this allele in total,#subjects in this cluster with HLA,#subjects in this cluster with this allele\n");
	}

	t_vector *target_cdr3_vector = global_options->target_cdr3_vector;

	svi_hash *cdr3toindex = kh_init(svim);
	for(int i = 0; i < kv_size(*target_cdr3_vector); i++)
	{
		tcr_t *entry = kv_A(*target_cdr3_vector, i);
		char* cdr3 = entry->cdr3a;
		if(strcmp(chain, CHAIN_BD) == 0)
		{
			cdr3 = entry->cdr3b;
		}
		append_to_svi_hash(cdr3toindex, cdr3, i);
	}
	char netfile[MAX_FILENAME_LENGTH];
	sprintf(netfile, "%s_%s.txt", global_options->out_prefix, chain);
	FILE *nfp = fopen(netfile, "r");
	if(nfp == NULL)
	{
		fprintf(stderr, "Could not open file %s for reading\n", netfile);
		exit(EXIT_FAILURE);
	}

	char *line = 0;
	char vgene[MAX_CDR3_LENGTH];
	while ((line = ic_readline(nfp)) != NULL)
	{
		if(line[strlen(line) - 1] == '\r' || line[strlen(line) - 1] =='\n')
		{
			line[strlen(line) - 1] = 0;
		}
		pattern_t *p = line_to_pattern(line);
		get_v_gene(p->pattern, vgene);
		si_hash *sid_hash = kh_init(sim);
		si_hash *pair_hash = kh_init(sim);
		si_hash *allele_hash = kh_init(sim);

		for(int i = 0; i < kv_size(*(p->members)); i++)
		{
			char* cdr3 = (char*)kv_A(*(p->members), i);
			i_vector *idx_vec = get_value_svi_hash(cdr3toindex, cdr3);
			for(int j = 0; j<kv_size(*idx_vec); j++)
			{
				int jj = kv_A(*idx_vec, j);
				tcr_t *entry = kv_A(*target_cdr3_vector, jj);
				char *v = entry->va;
				if(strcmp(chain, CHAIN_BD) == 0)
				{
					v = entry->vb;
				}
				if(global_options->same_v == 1 && strcmp(v, vgene) != 0)
				{
					continue;
				}
				khiter_t kk = kh_get(svsm, hla_hash, entry->sid);
				if(kk != kh_end(hla_hash))
				{
					append_to_si_hash(sid_hash, entry->sid, 0);
					s_vector *allele_vector = kh_val(hla_hash, kk);
					for (int a = 0; a < kv_size(*allele_vector); a++)
					{
						char *allele = kv_A(*allele_vector, a);
						sprintf(pair, "%s,%s", entry->sid, allele);
						append_to_si_hash(pair_hash, pair, 0);
					}
				}
			}
		}
		int sample_population = kh_size(sid_hash);
		destruct_si_hash(sid_hash);
		for (khiter_t k = kh_begin(pair_hash); k != kh_end(pair_hash); k++)
		{
			if (kh_exist(pair_hash, k))
			{
				char *sid_allele = (char*)kh_key(pair_hash, k);
				s_vector vec;
				kv_init(vec);
				ic_split_string(sid_allele, ',', &vec);
				char* allele = kv_A(vec, 1);
				// dbg_print("allele %s\n", allele);
				add_to_si_hash(allele_hash, allele);
				destruct_s_vector(&vec);
			}
		}
		destruct_si_hash(pair_hash);
		for (khiter_t k = kh_begin(allele_hash); k != kh_end(allele_hash); k++)
		{
			if (kh_exist(allele_hash, k))
			{
				char *allele = (char*)kh_key(allele_hash, k);
				int sample_population_with_feature = kh_val(allele_hash, k);
				int sub_population_with_feature = get_value_si_hash(all_allele_hash, allele, 0);
				float pvalue = compute_enrichment_pvalue(population, sub_population_with_feature, sample_population, sample_population_with_feature);
				if(pvalue < global_options->hla_association_pvalue_cutoff)
				{
					fprintf(hla_out_fp, "%s,%s,%s,%.1e,%d,%d,%d,%d\n", chain, p->pattern, allele, pvalue, population, sub_population_with_feature, sample_population, sample_population_with_feature);
				}
			}
		}
		destruct_si_hash(allele_hash);
		destruct_pattern(p);
	}

	fprintf(hla_out_fp, "\n");
	destruct_si_hash(all_allele_hash);
	// destruct_svs_hash(hla_hash);
	fclose(hla_out_fp);
	fclose(nfp);
	destruct_svi_hash(cdr3toindex);
}

void output_one_cluster(parameter_t *global_options, svi_hash *cdr3toindex, const char* chain, pattern_t *p, int lnum, FILE *ofp)
{
	int ignored_left = global_options->bd_ignored_v_end;
	int ignored_rght = global_options->bd_ignored_j_end;

	if (strcmp(chain, CHAIN_AG) == 0)
	{
		ignored_left = global_options->ag_ignored_v_end;
		ignored_rght = global_options->ag_ignored_j_end;
	}
	char sub_cdr3[MAX_CDR3_LENGTH];
	char vgene[MAX_CDR3_LENGTH];
	get_v_gene(p->pattern, vgene);
	char sample[H1];
	t_vector *target_cdr3_vector = global_options->target_cdr3_vector;
	s_set* uniq_cdr3 = kh_init(set);
	s_set* uniq_smpl = kh_init(set);

	si_vector vec;
	kv_init(vec);

	for(int i = 0; i < kv_size(*(p->members)); i++)
	{
		char* cdr3 = (char*)kv_A(*(p->members), i);
		i_vector *idx_vec = get_value_svi_hash(cdr3toindex, cdr3);
		for(int j = 0; j<kv_size(*idx_vec); j++)
		{
			int jj = kv_A(*idx_vec, j);
			tcr_t *entry = kv_A(*target_cdr3_vector, jj);
			si_t *obj = icalloc(1, si_t);
			char *rcdr3 = 0;
			if (strcmp(chain, CHAIN_AG) == 0)
			{
				if (global_options->same_v == 1 && strcmp(vgene, entry->va) != 0)
				{
					continue;
				}
				if(strcmp(entry->cdr3a, NA_SYMBOL) != 0)
				{
					append_to_s_set(uniq_cdr3, entry->cdr3a);
				}
				rcdr3 = entry->cdr3b;
			}
			else
			{
				if (global_options->same_v == 1 && strcmp(vgene, entry->vb) != 0)
				{
					continue;
				}
				rcdr3 = entry->cdr3a;
				if(strcmp(entry->cdr3b, NA_SYMBOL) != 0)
				{
					append_to_s_set(uniq_cdr3, entry->cdr3b);
				}
			}
			if(strcmp(rcdr3, NA_SYMBOL) == 0)
			{
				obj->key = strdup(rcdr3);
			}else{
				int pp = 0;
				for(int p = ignored_left; p < strlen(rcdr3) - ignored_rght; p++)
				{
					sub_cdr3[pp++] = rcdr3[p];
				}
				sub_cdr3[pp]= '\0';
				obj->key = strdup(sub_cdr3);
				// kv_push(char*, evec, strdup(sub_cdr3));
			}
			obj->val = jj;
			kv_push(si_t*, vec, obj);
			sprintf(sample, "%s %s", entry->sid, entry->condition);
			append_to_s_set(uniq_smpl, sample);
		}
	}
	sort_cluster(&vec);

	for(int j = 0; j<kv_size(vec); j++)
	{
		int jj = ((si_t*)kv_A(vec, j))->val;
		tcr_t *entry = kv_A(*target_cdr3_vector, jj);
		if(entry->freq + 0.05 > 1)
		{
			fprintf(ofp, "%d,%s,%.0e,%.2f,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d\n", lnum, p->pattern, p->pvalue, p->entfrc, kh_size(uniq_smpl), kh_size(uniq_cdr3), entry->cdr3a, entry->va, entry->ja, entry->cdr3b, entry->vb, entry->jb, entry->sid, entry->condition, (int)(entry->freq+0.05));
		}
		else
		{
			fprintf(ofp, "%d,%s,%.0e,%.2f,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%.1e\n", lnum, p->pattern, p->pvalue, p->entfrc, kh_size(uniq_smpl), kh_size(uniq_cdr3), entry->cdr3a, entry->va, entry->ja, entry->cdr3b, entry->vb, entry->jb, entry->sid, entry->condition, entry->freq);
		}
	}

	destruct_si_vector(&vec);
	destruct_s_set(uniq_cdr3);
	destruct_s_set(uniq_smpl);
}

void output_html(parameter_t *global_options, const char* chain, p_vector *pvec, vi_vector *viv)
{
	char html_file[MAX_FILENAME_LENGTH];
	sprintf(html_file, "%s_%s.html", global_options->out_prefix, chain);
	FILE *html_fp = fopen(html_file, "w");
	if (html_fp == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", html_file);
		exit(EXIT_FAILURE);
	}
	json_t *json = json_array();
	json_t *root=0,*node_arr=0,*link_arr=0;
	s_vector index_array;
	si_hash *unq_dict;
	for(int k0 = 0; k0 < kv_size(*viv); k0++)
	{
		i_vector *vec = kv_A(*viv, k0);
		if (kv_size(*vec) < 4) {
			break;
		}
		kv_init(index_array);
		unq_dict = kh_init(sim);
		root = json_object();
		json_array_append(json, root);
		node_arr = json_array();
		link_arr = json_array();
		json_object_set_new(root, "nodes", node_arr);
		json_object_set_new(root, "links", link_arr);
		for(int k1 = 0; k1 < kv_size(*vec); k1++)
		{
			pattern_t *pi = kv_A(*pvec, kv_A(*vec, k1));
			si_hash *hi = kh_init(sim);
			for(int ii = 0; ii < kv_size(*(pi->members)); ii++)
			{
				char* cdr3 = kv_A(*(pi->members), ii);
				add_to_si_hash(hi, cdr3);
			}
			for(int k2 = k1+1; k2 < kv_size(*vec); k2++)
			{
				pattern_t *pj = kv_A(*pvec, kv_A(*vec, k2));
				int overlap = overlap_cluster(hi, pj);
				if(overlap >= global_options->overlap_cutoff)
				{
					if(key_exist_si_hash(unq_dict, pi->pattern) == false)
					{
						kv_push(char*, index_array, strdup(pi->pattern));
						append_to_si_hash(unq_dict, pi->pattern, kv_size(index_array) - 1);
						json_t* node_obj = json_object();
						json_object_set_new(node_obj, "name", json_string(pi->pattern));
						json_object_set_new(node_obj, "group", json_integer(get_value_si_hash(unq_dict, pi->pattern, 0)));
						json_array_append(node_arr, node_obj);
					}

					if(key_exist_si_hash(unq_dict, pj->pattern) == false)
					{
						kv_push(char*, index_array, strdup(pj->pattern));
						append_to_si_hash(unq_dict, pj->pattern, kv_size(index_array) - 1);
						json_t* node_obj = json_object();
						json_object_set_new(node_obj, "name", json_string(pj->pattern));
						json_object_set_new(node_obj, "group", json_integer(get_value_si_hash(unq_dict, pj->pattern, 0)));
						json_array_append(node_arr, node_obj);
					}

					json_t* link_obj = json_object();
					json_object_set_new(link_obj, "source", json_integer(get_value_si_hash(unq_dict, pi->pattern, 0)));
					json_object_set_new(link_obj, "target", json_integer(get_value_si_hash(unq_dict, pj->pattern, 0)));
					json_object_set_new(link_obj, "weight", json_integer(overlap));
					json_array_append(link_arr, link_obj);
				}
			}
			destruct_si_hash(hi);
		}
		destruct_s_vector(&index_array);
		destruct_si_hash(unq_dict);
	}

	kv_init(index_array);
	unq_dict = kh_init(sim);
	root = json_object();
	json_array_append(json, root);
	node_arr = json_array();
	link_arr = json_array();
	json_object_set_new(root, "nodes", node_arr);
	json_object_set_new(root, "links", link_arr);

	for(int k0 = 0; k0 < kv_size(*viv); k0++)
	{
		i_vector *vec = kv_A(*viv, k0);

		if (kv_size(*vec) >= 4) {
			continue;
		}

		for(int k1 = 0; k1 < kv_size(*vec); k1++)
		{
			pattern_t *pi = kv_A(*pvec, kv_A(*vec, k1));
			si_hash *hi = kh_init(sim);
			for(int ii = 0; ii < kv_size(*(pi->members)); ii++)
			{
				char* cdr3 = kv_A(*(pi->members), ii);
				add_to_si_hash(hi, cdr3);
			}
			for(int k2 = k1+1; k2 < kv_size(*vec); k2++)
			{
				pattern_t *pj = kv_A(*pvec, kv_A(*vec, k2));
				int overlap = overlap_cluster(hi, pj);
				if(overlap >= global_options->overlap_cutoff)
				{
					if(key_exist_si_hash(unq_dict, pi->pattern) == false)
					{
						kv_push(char*, index_array, strdup(pi->pattern));
						append_to_si_hash(unq_dict, pi->pattern, kv_size(index_array) - 1);
						json_t* node_obj = json_object();
						json_object_set_new(node_obj, "name", json_string(pi->pattern));
						json_object_set_new(node_obj, "group", json_integer(get_value_si_hash(unq_dict, pi->pattern, 0)));
						json_object_set_new(node_obj, "weight", json_integer(kv_size(*(pi->members))));
						json_array_append(node_arr, node_obj);
					}

					if(key_exist_si_hash(unq_dict, pj->pattern) == false)
					{
						kv_push(char*, index_array, strdup(pj->pattern));
						append_to_si_hash(unq_dict, pj->pattern, kv_size(index_array) - 1);
						json_t* node_obj = json_object();
						json_object_set_new(node_obj, "name", json_string(pj->pattern));
						json_object_set_new(node_obj, "group", json_integer(get_value_si_hash(unq_dict, pj->pattern, 0)));
						json_object_set_new(node_obj, "weight", json_integer(kv_size(*(pj->members))));
						json_array_append(node_arr, node_obj);
					}

					json_t* link_obj = json_object();
					json_object_set_new(link_obj, "source", json_integer(get_value_si_hash(unq_dict, pi->pattern, 0)));
					json_object_set_new(link_obj, "target", json_integer(get_value_si_hash(unq_dict, pj->pattern, 0)));
					json_object_set_new(link_obj, "weight", json_integer(overlap));
					json_array_append(link_arr, link_obj);
				}
			}
			destruct_si_hash(hi);
		}

	}
	destruct_s_vector(&index_array);
	destruct_si_hash(unq_dict);
	char *s = json_dumps(json, 0);
	fprintf(html_fp, "<!DOCTYPE html>\n\
<!-- adapted from  http://bl.ocks.org/jose187/4733747 -->\n\
<meta charset=\"utf-8\">\n\
<script src=\"https://d3js.org/d3.v3.min.js\" charset=\"utf-8\"></script>\n\
<script src=\"https://code.jquery.com/jquery-3.6.0.slim.min.js\"></script>\n\
<style>\n\
\n\
.link {\n\
stroke: #aaa;\n\
}\n\
\n\
.node text {\n\
stroke:#333;\n\
font-family: Courier New;\n\
font-size: 10px;\n\
cursos:pointer;\n\
}\n\
\n\
.node circle{\n\
stroke:#fff;\n\
stroke-width:3px;\n\
fill:#555;\n\
}\n\
\n\
</style>\n\
<body>\n\
<pre id=\"rdata\" hidden>%s</pre>\n\
<script>\n\
var link;\n\
var node;\n\
function plot(idx){\n\
var width = window.innerWidth;\n\
var height = window.innerHeight;\n\
d3.select(\"svg\").remove();\n\
var force = d3.layout.force()\n\
.size([width, height])\n\
.charge(-400)\n\
.linkDistance(40);\n\
var drag = force.drag()\n\
.on(\"dragstart\", dragstart);\n\
var svg = d3.select(\"body\").append(\"svg\")\n\
.attr(\"width\", width)\n\
.attr(\"height\", height);\n\
link = svg.selectAll(\".link\");\n\
node = svg.selectAll(\".node\");\n\
var data =JSON.parse($(\"#rdata\").html());\n\
var graph = data[idx];\n\
force\n\
.nodes(graph.nodes)\n\
.links(graph.links)\n\
.start();\n\
link0 = link.data(graph.links)\n\
.enter().append(\"g\");\n\
link = link0.append(\"line\")\n\
.attr(\"class\", \"link\")\n\
.style(\"stroke-width\", function(d) { return Math.sqrt(d.weight - 60 + 1); });\n\
node = node.data(graph.nodes)\n\
.enter()\n\
.append(\"g\")\n\
.attr(\"class\", \"node\")\n\
.on(\"dblclick\", dblclick)\n\
.call(drag);\n\
node\n\
.append(\"circle\")\n\
.attr(\"r\",function(d) { return d.weight+4;});\n\
node\n\
.append(\"text\")\n\
.attr(\"dx\", 12)\n\
.attr(\"dy\", \".35em\")\n\
.attr(\"fill\",\"red\")\n\
.attr(\"font-family\",\"Courier New\")\n\
.attr(\"font-size\",\"10px\")\n\
.text(function(d) { return d.name });\n\
var linkText = link0\n\
.append(\"text\")\n\
.attr(\"class\", \"link-label\")\n\
.attr(\"font-family\", \"Arial, Helvetica, sans-serif\")\n\
.attr(\"fill\", \"black\")\n\
.style(\"font\", \"normal 10px Arial\")\n\
.attr(\"dy\", \".0em\")\n\
.attr(\"text-anchor\", \"middle\")\n\
.text(function(d) {\n\
return d.weight+'%%';\n\
});\n\
force.on(\"tick\", function() {\n\
link\n\
.attr(\"x1\", function(d) { return d.source.x; })\n\
.attr(\"y1\", function(d) { return d.source.y; })\n\
.attr(\"x2\", function(d) { return d.target.x; })\n\
.attr(\"y2\", function(d) { return d.target.y; });\n\
node\n\
.attr(\"transform\", function(d) { return \"translate(\" + d.x + \",\" + d.y + \")\"; });\n\
linkText\n\
.attr(\"x\", function(d) {\n\
return ((d.source.x + d.target.x)/2);\n\
})\n\
.attr(\"y\", function(d) {\n\
return ((d.source.y + d.target.y)/2);\n\
});\n\
});\n\
}\n\
function dblclick(d) {\n\
d3.select(this).classed(\"fixed\", d.fixed = false);\n\
}\n\
function dragstart(d) {\n\
d3.select(this).classed(\"fixed\", d.fixed = true);\n\
}\n\
var data =JSON.parse($(\"#rdata\").html());\n\
var links = '<center>';\n\
for(var i = 0; i < data.length; i++)\n\
{\n\
links += '<a href=\"#\" onclick=\"plot('+i+')\">'+(i+1)+'</a>&nbsp;';\n\
if((i+1) %% 50==0)\n\
{\n\
links += '<br/>';\n\
}\n\
}\n\
links += '</center>';\n\
$(\"body\").append(links);\n\
</script>\n\
</body>\n\
</html>\n\
", s);
	safe_free(s);
	json_decref(json);
}

void create_overlap_matrix(parameter_t *global_options, p_vector *pvec, vi_vector *matrix)
{
	ivi_hash *tmp = kh_init(ivim);
	int cdr3_index = 0;
	si_hash *index_dct = kh_init(sim);
	i_vector *vec;
	for(int i = 0; i < kv_size(*pvec); i++)
	{
		vec = icalloc(1, i_vector);
		kv_init(*vec);
		kv_push(i_vector*, *matrix, vec);
		pattern_t *pi = kv_A(*pvec, i);
		for(int ii = 0; ii < kv_size(*(pi->members)); ii++)
		{
			char* cdr3 = kv_A(*(pi->members), ii);
			if(key_exist_si_hash(index_dct, cdr3) == false)
			{
				append_to_si_hash_nd(index_dct, cdr3, cdr3_index++);
			}
			append_to_ivi_hash(tmp, get_value_si_hash(index_dct, cdr3, 0), i);
		}
	}

	si_hash *punq = kh_init(sim);
	char idx_pair[H1];

	for (khiter_t k = kh_begin(tmp); k != kh_end(tmp); k++)
	{
		if (!kh_exist(tmp, k))
		{
			continue;
		}
		i_vector *vector = kh_val(tmp, k);
		for(int i0 = 0; i0 < kv_size(*vector); i0++)
		{
			int i = kv_A(*vector, i0);
			pattern_t *pi = kv_A(*pvec, i);
			si_hash *hi = kh_init(sim);
			for(int ii = 0; ii < kv_size(*(pi->members)); ii++)
			{
				char* cdr3 = kv_A(*(pi->members), ii);
				add_to_si_hash(hi, cdr3);
			}
			for(int j0 = i0+1; j0 < kv_size(*vector); j0++)
			{
				int j = kv_A(*vector, j0);
				sprintf(idx_pair, "%d %d", i, j);
				if(key_exist_si_hash(punq, idx_pair) == true)
				{
					continue;
				}
				add_to_si_hash(punq, idx_pair);

				pattern_t *pj = kv_A(*pvec, j);
				int overlap = overlap_cluster(hi, pj);
				if(overlap >= global_options->overlap_cutoff)
				{
					vec = kv_A(*matrix, i);
					kv_push(int, *vec, j);
					vec = kv_A(*matrix, j);
					kv_push(int, *(vec), i);
				}
			}
			destruct_si_hash(hi);
		}
	}
	destruct_si_hash(punq);
	destruct_ivi_hash(tmp);
	destruct_si_hash_nd(index_dct);

	// for(int i = 0; i < kv_size(*pvec); i++)
	// {
	//      pattern_t *pi = kv_A(*pvec, i);
	//      si_hash *hi = kh_init(sim);
	//      for(int ii = 0; ii < kv_size(*(pi->members)); ii++)
	//      {
	//              char* cdr3 = kv_A(*(pi->members), ii);
	//              add_to_si_hash(hi, cdr3);
	//      }
	//      for(int j = i+1; j < kv_size(*pvec); j++)
	//      {
	//              pattern_t *pj = kv_A(*pvec, j);
	//              int overlap = overlap_cluster(hi, pj);
	//              if(overlap >= global_options->overlap_cutoff)
	//              {
	//                      vec = kv_A(*matrix, i);
	//                      kv_push(int, *vec, j);
	//                      vec = kv_A(*matrix, j);
	//                      kv_push(int, *(vec), i);
	//              }
	//      }
	//      destruct_si_hash(hi);
	// }
}

void output_cluster(parameter_t *global_options, const char* chain)
{
	char cluster_out_file[MAX_FILENAME_LENGTH];
	sprintf(cluster_out_file, "%s_%s_cluster.csv", global_options->out_prefix, chain);
	FILE *cluster_out_fp = fopen(cluster_out_file, "w");
	if (cluster_out_fp == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", cluster_out_file);
		exit(EXIT_FAILURE);
	}
	fprintf(cluster_out_fp, "index,pattern,pvalue,entropy_fraction,#unique_sample,#unique_cdr3,cdr3a,va,ja,cdr3b,vb,jb,sid,condition,frequency\n");

	t_vector *target_cdr3_vector = global_options->target_cdr3_vector;

	svi_hash *cdr3toindex = kh_init(svim);
	for(int i = 0; i < kv_size(*target_cdr3_vector); i++)
	{
		tcr_t *entry = kv_A(*target_cdr3_vector, i);
		char* cdr3 = entry->cdr3a;
		if(strcmp(chain, CHAIN_BD) == 0)
		{
			cdr3 = entry->cdr3b;
		}
		append_to_svi_hash(cdr3toindex, cdr3, i);
	}

	p_vector pvec;
	kv_init(pvec);

	char netfile[MAX_FILENAME_LENGTH];
	sprintf(netfile, "%s_%s.txt", global_options->out_prefix, chain);
	FILE *nfp = fopen(netfile, "r");
	if(nfp == NULL)
	{
		fprintf(stderr, "Could not open file %s for reading\n", netfile);
		exit(EXIT_FAILURE);
	}

	char *line = 0;
	while ((line = ic_readline(nfp)) != NULL)
	{
		if(line[strlen(line) - 1] == '\r' || line[strlen(line) - 1] =='\n')
		{
			line[strlen(line) - 1] = 0;
		}
		pattern_t *p = line_to_pattern(line);
		safe_free(line);
		kv_push(pattern_t*, pvec, p);
	}
	fclose(nfp);
	debug("constructing overlapping matrix");
	vi_vector matrix;
	kv_init(matrix);
	create_overlap_matrix(global_options, &pvec, &matrix);
	// i_vector *vec;
	// for(int i = 0; i < kv_size(pvec); i++)
	// {
	//      vec = icalloc(1, i_vector);
	//      kv_init(*vec);
	//      kv_push(i_vector*, matrix, vec);
	// }
	//
	// for(int i = 0; i < kv_size(pvec); i++)
	// {
	//      pattern_t *pi = kv_A(pvec, i);
	//      si_hash *hi = kh_init(sim);
	//      for(int ii = 0; ii < kv_size(*(pi->members)); ii++)
	//      {
	//              char* cdr3 = kv_A(*(pi->members), ii);
	//              add_to_si_hash(hi, cdr3);
	//      }
	//      for(int j = i+1; j < kv_size(pvec); j++)
	//      {
	//              pattern_t *pj = kv_A(pvec, j);
	//              int overlap = overlap_cluster(hi, pj);
	//              if(overlap >= global_options->overlap_cutoff)
	//              {
	//                      vec = kv_A(matrix, i);
	//                      kv_push(int, *vec, j);
	//                      vec = kv_A(matrix, j);
	//                      kv_push(int, *(vec), i);
	//              }
	//      }
	//      destruct_si_hash(hi);
	// }
	vi_vector viv;
	kv_init(viv);
	debug("searching connected components");
	connected_components(&matrix, &viv);
	destruct_vi_vector(&matrix);

	int lnum = 0;
	ii_hash *used = kh_init(iim);
	debug("output clusters ... ...");
	for(int i = 0; i < kv_size(viv); i++)
	{
		i_vector *vec = kv_A(viv, i);
		sort_i_vector(vec);
		for(int j = 0; j < kv_size(*vec); j++)
		{
			pattern_t *p = kv_A(pvec, kv_A(*vec, j));
			add_to_ii_hash(used, kv_A(*vec, j), 1);
			output_one_cluster(global_options, cdr3toindex, chain, p, lnum++, cluster_out_fp);
			fprintf(cluster_out_fp, "\n");
		}
		fprintf(cluster_out_fp, "\n");
	}

	for(int i = 0; i < kv_size(pvec); i++)
	{
		if(key_exist_ii_hash(used, i) == false)
		{
			pattern_t *p = kv_A(pvec, i);
			output_one_cluster(global_options, cdr3toindex, chain, p, lnum++, cluster_out_fp);
			fprintf(cluster_out_fp, "\n");
		}
	}
	destruct_ii_hash(used);
	fclose(cluster_out_fp);

	destruct_svi_hash(cdr3toindex);
	if(global_options->plot == 1)
	{
		output_html(global_options, chain, &pvec, &viv);
	}
	destruct_p_vector(&pvec);
	destruct_vi_vector(&viv);
}

bool match_pattern(const char* cdr3, const char* pattern, int ignored_left, int ignored_rght)
{
	for(int p = ignored_left; p < strlen(cdr3) - ignored_rght; p++)
	{
		if(pattern[p-ignored_left] != cdr3[p] && pattern[p-ignored_left] != '%')
		{
			return(false);
		}
	}
	return(true);
}

int overlap_cluster(si_hash *h1, pattern_t *p2)
{

	int n0=0, n1 = 0, n2 = 0;
	n1 = kh_size(h1);
	n2 = kv_size(*(p2->members));

	for(int i = 0; i < n2; i++)
	{
		char* cdr3 = (char*)kv_A(*(p2->members), i);
		if(key_exist_si_hash(h1, cdr3)== true)
		{
			n0++;
		}

	}

	if(n1 > n2)
	{
		return (n0 * 100 / n2);
	}else{
		return (n0 * 100 / n1);
	}
}

int count_unique_target_cdr3(parameter_t *global_options, const char* chain)
{
	int number = 0;
	int max_len = global_options->ag_cdr3_length_max_cutoff;
	int min_len = global_options->ag_cdr3_length_min_cutoff;
	t_vector *vector = global_options->target_cdr3_vector;

	if(strcmp(CHAIN_BD, chain) == 0)
	{
		max_len = global_options->bd_cdr3_length_max_cutoff;
		min_len = global_options->bd_cdr3_length_min_cutoff;
	}

	for(int i = 0; i < kv_size(*vector); i++)
	{
		tcr_t *tcr = kv_A(*vector, i);
		char *cdr3 = tcr->cdr3a;

		if(strcmp(CHAIN_BD, chain) == 0)
		{
			cdr3 = tcr->cdr3b;
		}

		if(strlen(cdr3) >= min_len && strlen(cdr3) <= max_len)
		{
			number++;
		}
	}

	return number;
}

int count_unique_refer_cdr3(parameter_t *global_options, const char* chain)
{
	int number = 0;
	int max_len = global_options->ag_cdr3_length_max_cutoff;
	int min_len = global_options->ag_cdr3_length_min_cutoff;
	c_vector *vector = global_options->ag_refer_cdr3_vector;

	if(strcmp(CHAIN_BD, chain) == 0)
	{
		max_len = global_options->bd_cdr3_length_max_cutoff;
		min_len = global_options->bd_cdr3_length_min_cutoff;
		vector = global_options->bd_refer_cdr3_vector;
	}

	for(int i = 0; i < kv_size(*vector); i++)
	{
		cdr_t *cdr = kv_A(*vector, i);
		if(strlen(cdr->cdr3) >= min_len && strlen(cdr->cdr3) <= max_len)
		{
			number++;
		}
	}

	return number;
}


void prepare_data(parameter_t *global_options)
{
	global_options->unique_tra_target_total=0;
	global_options->unique_trb_target_total=0;
	global_options->unique_tra_refer_total=0;
	global_options->unique_trb_refer_total=0;
	global_options->hla_availability = true;


	format_input_cdr3_file(global_options);
	if(strlen(global_options->ag_refer_file) > 0 && access(global_options->ag_refer_file, R_OK) == 0)
	{
		global_options->ag_refer_cdr3_vector = icalloc(1, c_vector);
		kv_init(*global_options->ag_refer_cdr3_vector);
		read_refer_cdr3(global_options->ag_refer_file, global_options->ag_refer_cdr3_vector);
		global_options->unique_tra_refer_total=count_unique_refer_cdr3(global_options, CHAIN_AG);
	}

	if(strlen(global_options->bd_refer_file) > 0 && access(global_options->bd_refer_file, R_OK) == 0)
	{
		global_options->bd_refer_cdr3_vector = icalloc(1, c_vector);
		kv_init(*(global_options->bd_refer_cdr3_vector));
		read_refer_cdr3(global_options->bd_refer_file, global_options->bd_refer_cdr3_vector);
		global_options->unique_trb_refer_total=count_unique_refer_cdr3(global_options, CHAIN_BD);
	}

	global_options->target_cdr3_vector = icalloc(1, t_vector);
	kv_init(*(global_options->target_cdr3_vector));
	read_target_cdr3(global_options, global_options->target_cdr3_vector);

	global_options->unique_tra_target_total=count_unique_target_cdr3(global_options, CHAIN_AG);
	global_options->unique_trb_target_total=count_unique_target_cdr3(global_options, CHAIN_BD);

	if(strlen(global_options->hla_file) == 0 || access(global_options->hla_file, R_OK) != 0)
	{
		global_options->hla_availability = false;
	}
	else
	{
		read_hla(global_options);
	}

}

void sort_cluster_on_entfrc(const char* ifile)
{
	debug("sort_cluster_on_entfrc");
	FILE *ifp = fopen(ifile, "r");
	if (ifp == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", ifile);
		exit(EXIT_FAILURE);
	}
	// debug("sort_cluster_on_entfrc.1");
	s_vector tmpfile_vector;
	kv_init(tmpfile_vector);

	char tmpfile[MAX_FILENAME_LENGTH], cmd[K1];
	int index = 0, pid = getpid();
	FILE *ofp = 0;
	char *line = 0;
	int lnum = 0;
	p_vector pvec;
	kv_init(pvec);
	while ((line = ic_readline(ifp)) != NULL)
	{
		if(line[strlen(line) - 1] == '\r' || line[strlen(line) - 1] =='\n')
		{
			line[strlen(line) - 1] = 0;
		}
		lnum++;
		if(lnum % 10000 == 0)
		{
			debug("processing %d lines", lnum);
		}
		pattern_t *p = line_to_pattern(line);
		kv_push(pattern_t*, pvec, p);
		if (lnum % NUM_OF_LINE_TO_SORT == 0)
		{
			sort_p_vector(&pvec);
			index++;
			sprintf(tmpfile, "/tmp/.ictools_%d_%d", pid, index);
			debug("write to file %s", tmpfile);
			kv_push(char*, tmpfile_vector, strdup(tmpfile));
			ofp = fopen(tmpfile, "w");
			if (ofp == NULL)
			{
				fprintf(stderr, "Could not open file %s for wrting\n", tmpfile);
				exit(EXIT_FAILURE);
			}
			for(int i = 0; i < kv_size(pvec); i++)
			{
				pattern_t *p = (pattern_t*)kv_A(pvec, i);
				kstring_t *kstr = pattern_to_kstring(p);
				fprintf(ofp, "%s", kstr->s);
				safe_free(kstr->s);
				safe_free(kstr);
			}
			destruct_p_vector(&pvec);
			kv_init(pvec);
			fclose(ofp);
		}
	}
	fclose(ifp);
	// debug("sort_cluster_on_entfrc.2");
	if(kv_size(pvec) > 0)
	{
		sort_p_vector(&pvec);
		index++;
		sprintf(tmpfile, "/tmp/.ictools_%d_%d", pid, index);
		kv_push(char*, tmpfile_vector, strdup(tmpfile));
		ofp = fopen(tmpfile, "w");
		if (ofp == NULL)
		{
			fprintf(stderr, "Could not open file %s for wrting\n", tmpfile);
			exit(EXIT_FAILURE);
		}
		for(int i = 0; i < kv_size(pvec); i++)
		{
			pattern_t *p = (pattern_t*)kv_A(pvec, i);
			kstring_t *kstr = pattern_to_kstring(p);
			fprintf(ofp, "%s", kstr->s);
			safe_free(kstr->s);
			safe_free(kstr);
		}
		destruct_p_vector(&pvec);
		fclose(ofp);

	}
	// debug("sort_cluster_on_entfrc.3");
	int num_of_file = kv_size(tmpfile_vector);
	FILE **    fp_array    = (FILE **)malloc(sizeof(FILE *) * num_of_file);
	heap_t *   hp          = init_heap();
	s_vector fields;
	float entfrc = 0.0;
	hnode_t *  node;
	char *     fname;
	for (int i = 0; i < num_of_file; ++i)
	{
		fname        = kv_A(tmpfile_vector, i);
		ifp          = fopen(fname, "r");
		fp_array[i] = ifp;
		line        = ic_readline(ifp);
		if (line != NULL)
		{
			node = icalloc(1, hnode_t);
			kv_init(fields);
			ic_split_string(line, 0, &fields);
			entfrc      = atof(kv_A(fields, 2));
			node->data  = line;
			node->key   = entfrc;
			node->index = i;
			insert_heap_node(hp, node);
			destruct_s_vector(&fields);
		}
	}
	// debug("sort_cluster_on_entfrc.4");
	ofp = fopen(ifile, "w");
	while (hp->size > 0)
	{
		// debug("sort_cluster_on_entfrc.401");
		node = pop_heap_node(hp);
		fprintf(ofp, "%s", node->data);
		index = node->index;
		line = ic_readline(fp_array[index]);
		safe_free(node->data);
		safe_free(node);
		if (line != NULL)
		{
			node = icalloc(1, hnode_t);
			kv_init(fields);
			ic_split_string(line, 0, &fields);
			entfrc      = atof(kv_A(fields, 2));
			node->data  = line;
			node->key   = entfrc;
			node->index = index;
			insert_heap_node(hp, node);
			destruct_s_vector(&fields);

		}
		// debug("sort_cluster_on_entfrc.402");
	}
	// debug("sort_cluster_on_entfrc.4.1");
	fclose(ofp);
	for (int i = 0; i < num_of_file; ++i)
	{
		fclose(fp_array[i]);
		fname = kv_A(tmpfile_vector, i);
		sprintf(cmd, "rm %s", fname);
		system_call(cmd);
	}
	// debug("sort_cluster_on_entfrc.4.2");
	destruct_s_vector(&tmpfile_vector);
	// debug("sort_cluster_on_entfrc.4.3");
	destroy_heap(hp);
	// debug("sort_cluster_on_entfrc.5");
}

void combinations(i_vector *iterable, int r, iv_vector *rslt)
{
	int n = kv_size(*iterable);
	assert(r <= n);
	int i = 0;
	int *indices = icalloc(r, int);
	i_vector *tmp = icalloc(1, i_vector);
	kv_init(*tmp);
	for (i = 0; i < r; i++)
	{
		indices[i] = i;
		kv_push(int, *tmp, kv_A(*iterable, indices[i]));
	}
	kv_push(i_vector*, *rslt, tmp);

	int normal_end = 1;

	while(1)
	{
		normal_end = 1;
		for (i = r-1; i >=0; i--)
		{
			if (indices[i] != i + n - r)
			{
				normal_end = 0;
				break;
			}
		}
		if (normal_end == 1)
		{
			break;
		}
		indices[i] += 1;
		for(int j = i+1; j < r; j++)
		{
			indices[j] = indices[j-1] + 1;
		}
		tmp = icalloc(1, i_vector);
		kv_init(*tmp);
		for (i = 0; i < r; i++)
		{
			kv_push(int, *tmp, kv_A(*iterable, indices[i]));
		}
		kv_push(i_vector*, *rslt, tmp);
	}
}

void get_v_gene(const char* pattern, char* vgene)
{
	s_vector fields;
	kv_init(fields);
	ic_split_string(pattern, ':', &fields);
	sprintf(vgene, "%s", kv_A(fields, 0));
	destruct_s_vector(&fields);
}

bool within(i_vector *small, i_vector *large)
{
	for(int i = 0; i < kv_size(*small); i++)
	{
		int sent = kv_A(*small, i);
		bool found = false;
		for(int j = 0; j < kv_size(*large); j++)
		{
			int lent = kv_A(*large, j);
			if (sent == lent)
			{
				found = true;
				break;
			}
		}
		if (found == false)
		{
			return false;
		}
	}
	return true;
}

// if mstr1 >= mstr2: return 1
// elif mstr1 < mstr2: return -1
// else return 0
int issubset(i_vector *v1, i_vector *v2)
{
	int rlst = 0;

	if (kv_size(*v1) >= kv_size(*v2))
	{
		if(within(v2, v1) == true)
		{
			rlst = 1;
		}
	}else{
		if(within(v1, v2) == true)
		{
			rlst = -1;
		}
	}
	return rlst;
}

int purge_cluster(parameter_t *global_options, s_vector *ivec, ii_hash *dict)
{
	s_vector tmp;
	kv_init(tmp);
	for(int i = 0; i < kv_size(*ivec); i++)
	{
		if(key_exist_ii_hash(dict, i) == true)
		{
			continue;
		}
		kv_push(char*, tmp, kv_A(*ivec, i));
		debug("===%d %s", i, kv_A(*ivec, i));
	}
	float sofar_val = entropy_fraction(&tmp);
	float best_val = sofar_val;
	debug("=== %f", best_val);
	kv_destroy(tmp);

	int best_pos = -1;
	for(int i = 0; i < kv_size(*ivec); i++)
	{
		if(key_exist_ii_hash(dict, i) == true)
		{
			continue;
		}
		kv_init(tmp);
		for(int j = 0; j < kv_size(*ivec); j++)
		{
			if(i == j || key_exist_ii_hash(dict, j) == true) {
				continue;
			}
			kv_push(char*, tmp, kv_A(*ivec, j));
		}
		if(kv_size(tmp) < OTHER_CHAIN_CUTOFF)
		{
			break;
		}
		float current_val = entropy_fraction(&tmp);
		debug("%f %d %s", current_val, i, kv_A(*ivec, i));
		if(current_val < best_val)
		{
			best_val = current_val;
			best_pos = i;
		}
		kv_destroy(tmp);
	}
	// smaller, better
	if(best_val < global_options->purge_fraction * sofar_val)
	{
		add_to_ii_hash(dict, best_pos, 1);
		return best_pos;
	}else
	{
		return -1;
	}
}

float entropy_fraction(s_vector *ivec)
{
	char kmer[ENTROPY_KMER_SIZE+1];
	kmer[ENTROPY_KMER_SIZE] = '\0';
	float tot = 0.0001;
	si_hash *freq = kh_init(sim);
	for(int i = 0; i < kv_size(*ivec); i++)
	{
		char *frag = kv_A(*ivec, i);
		// printf("frag is %s\n", frag);
		if(strlen(frag) < ENTROPY_KMER_SIZE)
		{
			tot += 1;
			add_to_si_hash(freq, frag);
			continue;
		}
		for(int j = 0; j < strlen(frag) - ENTROPY_KMER_SIZE + 1; j++)
		{
			strncpy(kmer, frag+j, ENTROPY_KMER_SIZE);
			tot += 1;
			add_to_si_hash(freq, kmer);
		}
	}

	float ent0 = -log2f(1/tot);
	float ent1 = 0.0;
	for (khiter_t k = kh_begin(freq); k != kh_end(freq); k++)
	{
		if (kh_exist(freq, k))
		{
			int count = kh_val(freq, k)+0.0001;
			char *kmer = kh_key(freq, k);
			ent1 -= count/tot * log2f(count/tot);
			// printf("ent0:%.1f kmer:%s(%d) ent1:%.1f\n", ent0, kmer, count, ent1);
		}
	}
	destruct_si_hash(freq);
	// printf("entropy fraction is %.3f\n", ent1/ent0);
	return ent1/ent0;
}

void connected_components(vi_vector *matrix, vi_vector *rslt)
{
	vi_vector container;
	kv_init(container);
	ii_hash *seen = kh_init(iim);
	queue_t *pq = icalloc(1, queue_t);
	int node;
	init_queue(pq);

	for (int root = 0; root < kv_size(*matrix); root++)
	{
		if(key_exist_ii_hash(seen, root) == false)
		{
			add_to_ii_hash(seen, root, 1);
			i_vector *component = icalloc(1, i_vector);
			kv_init(*component);
			enqueue(pq, root);
			while(!queue_is_empty(pq))
			{
				dequeue(pq, &node);
				kv_push(int, *component, node);
				i_vector *vec = kv_A(*matrix, node);
				for(int i = 0; i < kv_size(*vec); i++)
				{
					int neighbor = kv_A(*vec, i);
					if (key_exist_ii_hash(seen, neighbor) == false)
					{
						add_to_ii_hash(seen, neighbor, 1);
						enqueue(pq, neighbor);
					}
				}
			}
			kv_push(i_vector*, container, component);
		}
	}
	destruct_ii_hash(seen);
	safe_free(pq);
	/*
	      ii_hash *used_dict = kh_init(iim);
	      while(1)
	      {
	              bool done = true;
	              ii_hash *cc = kh_init(iim);
	              for (khiter_t k = kh_begin(matrix); k != kh_end(matrix); k++)
	              {
	                      if (kh_exist(matrix, k))
	                      {
	                              int node = kh_key(matrix, k);
	                              if (key_exist_ii_hash(used_dict, node) == true)
	                              {
	                                      continue;
	                              }
	                              done = false;
	                              add_to_ii_hash(cc, node, 1);
	                              add_to_ii_hash(used_dict, node, 1);
	                              i_vector *vec = kh_val(matrix, k);
	                              for(int i = 0; i < kv_size(*vec); i++)
	                              {
	                                      int _node = kv_A(*vec, i);
	                                      add_to_ii_hash(cc, _node, 1);
	                                      add_to_ii_hash(used_dict, _node, 1);
	                              }
	                              break;
	                      }
	              }
	              while(1)
	              {
	                      int orig_size = kh_size(cc);
	                      ii_hash *tmp = kh_init(iim);
	                      for (khiter_t k = kh_begin(cc); k != kh_end(cc); k++)
	                      {
	                              if (kh_exist(cc, k))
	                              {
	                                      int node = kh_key(cc, k);
	                                      i_vector *vec = get_value_ivi_hash(matrix, node);
	                                      for(int i = 0; i < kv_size(*vec); i++)
	                                      {
	                                              int _node = kv_A(*vec, i);
	                                              add_to_ii_hash(tmp, _node, 1);
	                                              add_to_ii_hash(used_dict, _node, 1);
	                                      }
	                              }
	                      }
	                      for (khiter_t k = kh_begin(tmp); k != kh_end(tmp); k++)
	                      {
	                              if (kh_exist(tmp, k))
	                              {
	                                      int node = kh_key(tmp, k);
	                                      add_to_ii_hash(cc, node, 1);
	                              }
	                      }
	                      destruct_ii_hash(tmp);
	                      if(orig_size == kh_size(cc))
	                      {
	                              break;
	                      }
	              }

	              if(done != true)
	              {
	                      i_vector *vec = icalloc(1, i_vector);
	                      kv_init(*vec);
	                      for (khiter_t k = kh_begin(cc); k != kh_end(cc); k++)
	                      {
	                              if (kh_exist(cc, k))
	                              {
	                                      int node = kh_key(cc, k);
	                                      kv_push(int, *vec, node);
	                              }
	                      }
	                      kv_push(i_vector*, container, vec);
	              }
	              destruct_ii_hash(cc);
	              if (done == true)
	              {
	                      break;
	              }
	      }
	      destruct_ii_hash(used_dict);
	 */
	is_vector isv;
	kv_init(isv);

	for(int i = 0; i < kv_size(container); i++)
	{
		i_vector *vec = kv_A(container, i);
		is_t *ent = icalloc(1, is_t);
		ent->idx = i;
		ent->sze = kv_size(*vec);
		kv_push(is_t*, isv, ent);
	}
	sort_is_vector(&isv);
	for(int i = 0; i < kv_size(isv); i++)
	{
		is_t *ent = kv_A(isv, i);
		kv_push(i_vector*, *rslt, kv_A(container, ent->idx));
	}
	kv_destroy(container);
	destruct_is_vector(&isv);
}
