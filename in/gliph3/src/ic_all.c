#include "ic_all.h"

void count_pattern(parameter_t *global_options, t_vector *tvec, c_vector *rvec, int cdr3_len, const char* chain, unsigned long pattern_mask, av_hash *avh)
{
	char *kmer = icalloc(cdr3_len+1, char);
	kmer[cdr3_len] = '\0';
	char *cdr3 = 0, *vgene = 0, *ocdr3 = 0;
	char sample[H5], vk[H2], vc[H2], _cdr3[H2], _ocdr3[H2], _buf[H2];
	ssi_hash *hash = kh_init(ssim);
	kstring_t *kstr = icalloc(1, kstring_t);
	si_hash *unq = kh_init(sim);
	si_hash *tdict = kh_init(sim);

	int oleft = global_options->ag_ignored_v_end, orght = global_options->ag_ignored_j_end;
	int left = global_options->bd_ignored_v_end, rght = global_options->bd_ignored_j_end;

	if (strcmp(chain, CHAIN_AG) == 0)
	{
		oleft = global_options->bd_ignored_v_end;
		orght = global_options->bd_ignored_j_end;
		left = global_options->ag_ignored_v_end;
		rght = global_options->ag_ignored_j_end;
	}

	bool continuous = true;
	for (int p = left; p < cdr3_len - rght; p++)
	{
		if(!(pattern_mask & (1 << p )))
		{
			continuous = false;
			break;
		}
	}

	for (int i = 0; i < kv_size(*tvec); i++)
	{
		tcr_t* entry = (tcr_t*)kv_A(*tvec, i);
		if (strcmp(chain, CHAIN_BD) == 0)
		{
			cdr3 = entry->cdr3b;
			vgene = entry->vb;
			ocdr3 = entry->cdr3a;
		}else{
			cdr3 = entry->cdr3a;
			vgene = entry->va;
			ocdr3 = entry->cdr3b;
		}

		for (int p = 0; p < cdr3_len; p++)
		{
			if(pattern_mask & (1 << p ))
			{
				kmer[p] = cdr3[p];
			}else{
				kmer[p] = '.';
			}
		}
		if(global_options->same_v == 1)
		{
			sprintf(vk, "%s:%s", vgene, kmer);
			sprintf(vc, "%s:%s", vgene, cdr3);
		}
		else
		{
			sprintf(vk, "%s", kmer);
			sprintf(vc, "%s", cdr3);
		}

		if(continuous == true)
		{
			add_to_si_hash(tdict, vk);
			add_to_si_hash(unq, vc);
		}else{
			if(key_exist_si_hash(unq, vc) == false)
			{
				add_to_si_hash(tdict, vk);
				add_to_si_hash(unq, vc);
			}
		}
		sprintf(sample, "%s@%s@%s@%s", cdr3, ocdr3, entry->sid, entry->condition);
		add_to_ssi_hash(hash, vk, sample);
	}
	int number_unique_target_cdr3 = kh_size(unq);
	int number_unique_refer_cdr3 = kv_size(*rvec);
	destruct_si_hash(unq);

	si_hash *rdict = kh_init(sim);
	for(int i = 0; i < kv_size(*rvec); i++)
	{
		cdr_t* cdr = kv_A(*rvec, i);
		for (int p = 0; p < cdr3_len; p++)
		{
			if(pattern_mask & (1 << p ))
			{
				kmer[p] = cdr->cdr3[p];
			}else{
				kmer[p] = '.';
			}
		}
		if(global_options->same_v == 1)
		{
			sprintf(vk, "%s:%s", cdr->v, kmer);
		}
		else
		{
			sprintf(vk, "%s", kmer);
		}

		if (size_content_for_key_ssi_hash(hash, vk) >= global_options->pattern_unique_sample_cutoff)
		{
			add_to_si_hash(rdict, vk);
		}
	}
	safe_free(kmer);

	if(global_options->purge_cluster == 1)
	{
		for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
		{
			if (kh_exist(hash, k))
			{
				char* pattern = (char*)kh_key(hash, k);
				si_hash *tmp = (si_hash*)kh_val(hash, k);
				if(kh_size(tmp) < global_options->pattern_unique_sample_cutoff)
				{
					continue;
				}

				s_vector ovec;
				kv_init(ovec);
				for (khiter_t kk = kh_begin(tmp); kk != kh_end(tmp); kk++)
				{
					if (kh_exist(tmp, kk))
					{
						char* smp = (char*)kh_key(tmp, kk);
						sscanf(smp, "%[^@]@%[^@]@", _cdr3, _ocdr3);
						if(strcmp(_ocdr3, NA_SYMBOL) != 0)
						{
							int jj = 0;
							for(int ii = oleft; ii < strlen(_ocdr3) - orght; ii++)
							{
								_buf[jj++] = _ocdr3[ii];
							}
							_buf[jj]='\0';
							kv_push(char*, ovec, strdup(_buf));
						}
					}
				}
				if(kv_size(ovec) > OTHER_CHAIN_CUTOFF)
				{
					debug("%s --------------------------", pattern);

					ii_hash *tmp_dct = kh_init(iim);
					int pos = purge_cluster(global_options, &ovec, tmp_dct);

					while(pos != -1)
					{
						pos = purge_cluster(global_options, &ovec, tmp_dct);
					}

					si_hash *s_dct = kh_init(sim);
					for (khiter_t kk = kh_begin(tmp_dct); kk != kh_end(tmp_dct); kk++)
					{
						if (kh_exist(tmp_dct, kk))
						{
							int pos = kh_key(tmp_dct, kk);
							char* core = kv_A(ovec, pos);
							add_to_si_hash(s_dct, core);
							debug("- %s", core);
						}
					}
					debug("--------------------------\n");
					destruct_ii_hash(tmp_dct);
					si_hash *n_tmp = kh_init(sim);
					for (khiter_t kk = kh_begin(tmp); kk != kh_end(tmp); kk++)
					{
						if (kh_exist(tmp, kk))
						{
							char* smp = (char*)kh_key(tmp, kk);
							int count = kh_val(tmp, kk);
							sscanf(smp, "%[^@]@%[^@]@", _cdr3, _ocdr3);
							if(strcmp(_ocdr3, NA_SYMBOL) != 0)
							{
								int jj = 0;
								for(int ii = oleft; ii < strlen(_ocdr3) - orght; ii++)
								{
									_buf[jj++] = _ocdr3[ii];
								}
								_buf[jj]='\0';
								if(key_exist_si_hash(s_dct, _buf) == true)
								{
									continue;
								}
							}
							append_to_si_hash(n_tmp, strdup(smp), count);
						}
					}
					destruct_si_hash(tmp);
					destruct_si_hash(s_dct);
					kh_val(hash, k) = n_tmp;
				}
				destruct_s_vector(&ovec);
			}
		}
	}

	for (khiter_t k = kh_begin(hash); k != kh_end(hash); k++)
	{
		if (kh_exist(hash, k))
		{
			char* pattern = (char*)kh_key(hash, k);
			si_hash *tmp = (si_hash*)kh_val(hash, k);
			if(kh_size(tmp) >= global_options->pattern_unique_sample_cutoff)
			{
				int target_count = get_value_si_hash(tdict, pattern, 0);
				int refer_count = get_value_si_hash(rdict, pattern, 0);
				double p_value = fisher_exact_test(target_count, number_unique_target_cdr3 - target_count, refer_count, number_unique_refer_cdr3 - refer_count);
				int ove = (number_unique_refer_cdr3 * (1+target_count)) / (number_unique_target_cdr3 * (1+refer_count));
				float frac = 1.0;
				bool other_present = false;
				s_vector ovec;
				kv_init(ovec);
				for (khiter_t kk = kh_begin(tmp); kk != kh_end(tmp); kk++)
				{
					if (kh_exist(tmp, kk))
					{
						char* smp = (char*)kh_key(tmp, kk);
						// int count = kh_val(tmp, kk);
						sscanf(smp, "%[^@]@%[^@]@", _cdr3, _ocdr3);
						if(strcmp(_ocdr3, NA_SYMBOL) != 0)
						{
							int jj = 0;
							for(int ii = oleft; ii < strlen(_ocdr3) - orght; ii++)
							{
								_buf[jj++] = _ocdr3[ii];
							}
							_buf[jj]='\0';
							kv_push(char*, ovec, strdup(_buf));
						}
					}
				}
				// updated on 05/06/22
				// change kv_size(ovec) > OTHER_CHAIN_CUTOFF
				// to
				// kv_size(ovec) >= OTHER_CHAIN_CUTOFF
				if(kv_size(ovec) >= OTHER_CHAIN_CUTOFF)
				{
					other_present = true;
					frac = entropy_fraction(&ovec);
					printf("*****ovec size is %d, EF is %.3f\n", kv_size(ovec), frac);
				}

				destruct_s_vector(&ovec);
				if(p_value > global_options->pattern_pvalue_cutoff || ove < global_options->pattern_ove_cutoff || (other_present && frac > global_options->entropy_fraction_cutoff))
				{
					continue;
				}
				all_t *pat = icalloc(1, all_t);
				pat->pattern = strdup(pattern);
				pat->p_value = p_value;
				pat->ove = ove;
				pat->number_dot = 0;
				for(int ii = 0; ii < cdr3_len; ii++)
				{
					if(pattern[ii] == '.')
					{
						pat->number_dot++;
					}
				}
				pat->entropy_fraction = frac;
				pat->target_count = target_count;
				pat->target_total = number_unique_target_cdr3;
				pat->refer_count = refer_count;
				pat->refer_total = number_unique_refer_cdr3;
				s_vector mbuf;
				kv_init(mbuf);

				si_hash *unq = kh_init(sim);

				for (khiter_t kk = kh_begin(tmp); kk != kh_end(tmp); kk++)
				{
					if (kh_exist(tmp, kk))
					{
						char* smp = (char*)kh_key(tmp, kk);
						sscanf(smp, "%[^@]@%[^@]@", _cdr3, _ocdr3);
						if(key_exist_si_hash(unq, _cdr3) == false)
						{
							kv_push(char*, mbuf, strdup(_cdr3));
							add_to_si_hash(unq, _cdr3);
						}
					}
				}

				destruct_si_hash(unq);
				sort_s_vector(&mbuf);
				kstr->l = 0;
				ksprintf(kstr, "%s", kv_A(mbuf, 0));
				for(int ii = 1; ii < kv_size(mbuf); ii++)
				{
					ksprintf(kstr, " %s", kv_A(mbuf, ii));
				}
				append_to_av_hash(avh, kstr->s, pat);
				destruct_s_vector(&mbuf);
			}
		}
	}
	safe_free(kstr->s);
	safe_free(kstr);
	destruct_si_hash(rdict);
	destruct_si_hash(tdict);
	destruct_ssi_hash(hash);
}


bool match_all_pattern(char* pattern, char* cdr3)
{
	assert(strlen(pattern) == strlen(cdr3));

	for(int p = 0; p < strlen(pattern); p++)
	{
		if(cdr3[p] != pattern[p] && pattern[p] != '.')
		{
			return false;
		}
	}
	return true;
}


void discover_enriched_pattern(parameter_t *global_options, const char* chain)
{
	char kmer_file[H5];
	sprintf(kmer_file, "%s_%s.txt", global_options->out_prefix, chain);
	FILE *kmer_fp = fopen(kmer_file, "w");
	if (kmer_fp == NULL)
	{
		fprintf(stderr, "Could not open file %s for writing\n", kmer_file);
		exit(EXIT_FAILURE);
	}

	int ignored_left = global_options->bd_ignored_v_end;
	int ignored_rght = global_options->bd_ignored_j_end;
	int min_length = global_options->bd_cdr3_length_min_cutoff;
	int max_length = global_options->bd_cdr3_length_max_cutoff;

	if(strcmp(chain, CHAIN_AG)==0)
	{
		ignored_left = global_options->ag_ignored_v_end;
		ignored_rght = global_options->ag_ignored_j_end;
		min_length = global_options->ag_cdr3_length_min_cutoff;
		max_length = global_options->ag_cdr3_length_max_cutoff;
	}

	s_vector fields;
	for(int l = min_length; l <= max_length; l++)
	{
		iv_vector rslt;
		kv_init(rslt);
		li_hash *pdict = kh_init(lim);
		i_vector iterable;
		kv_init(iterable);

		for(int p = ignored_left; p < l - ignored_rght; p++)
		{
			kv_push(int, iterable, p);
		}

		for (int m = global_options->min_motif_length; m <= global_options->max_motif_length; m++)
		{
			if(l - ignored_left - ignored_rght - m >= 0 )
			{
				combinations(&iterable, m, &rslt);
			}
		}
		for(int i = 0; i < kv_size(rslt); i++)
		{
			unsigned long pattern_mask = 0;
			i_vector *tmp = kv_A(rslt, i);
			sort_i_vector(tmp);
			int gap = 0;
			for(int j = 1; j < kv_size(*tmp); j++)
			{
				int val0 = kv_A(*tmp, j-1);
				int val1 = kv_A(*tmp, j);
				gap += val1 - val0 - 1;
			}

			if (global_options->continuous_motif == 1 && gap > 0)
			{
				continue;
			}

			for(int j = 0; j < kv_size(*tmp); j++)
			{
				int val = kv_A(*tmp, j);
				pattern_mask |= 1 << val;
			}
			add_to_li_hash(pdict, pattern_mask, 1);
		}

		destruct_iv_vector(&rslt);
		kv_init(rslt);
		for(int d = 0; d <= global_options->max_diff_position; d++)
		{
			int m = l - ignored_left - ignored_rght - d;
			combinations(&iterable, m, &rslt);
		}

		for(int i = 0; i < kv_size(rslt); i++)
		{
			unsigned long pattern_mask = 0;
			i_vector *tmp = kv_A(rslt, i);
			for(int j = 0; j < kv_size(*tmp); j++)
			{
				int val = kv_A(*tmp, j);
				pattern_mask |= 1 << val;
			}
			add_to_li_hash(pdict, pattern_mask, 1);
		}
		destruct_i_vector(&iterable);
		destruct_iv_vector(&rslt);
		t_vector tvec;
		kv_init(tvec);
		t_vector *tmptvector = global_options->target_cdr3_vector;
		for(int i = 0; i < kv_size(*tmptvector); i++)
		{
			tcr_t *tcr = kv_A(*tmptvector, i);
			if(strcmp(chain, CHAIN_AG) == 0 && strlen(tcr->cdr3a) == l)
			{
				kv_push(tcr_t*, tvec, tcr);
			}
			else if (strcmp(chain, CHAIN_BD) == 0 && strlen(tcr->cdr3b) == l)
			{
				kv_push(tcr_t*, tvec, tcr);
			}
		}

		c_vector rvec;
		kv_init(rvec);
		c_vector *tmprvector = global_options->ag_refer_cdr3_vector;
		if(strcmp(chain, CHAIN_BD) == 0)
		{
			tmprvector = global_options->bd_refer_cdr3_vector;
		}

		for(int i = 0; i < kv_size(*tmprvector); i++)
		{
			cdr_t *cdr = kv_A(*tmprvector, i);
			if(strlen(cdr->cdr3) == l)
			{
				kv_push(cdr_t*, rvec, cdr);
			}
		}

		av_hash *avh = kh_init(avm);
		int number_mask = kh_size(pdict);
		int current_pos = 0;
		for (khiter_t k = kh_begin(pdict); k != kh_end(pdict); k++)
		{
			if (kh_exist(pdict, k))
			{
				unsigned long pattern_mask = kh_key(pdict, k);
				current_pos++;
				debug("length: %d pattern mask: %d %d/%d", l, pattern_mask, current_pos, number_mask);
				count_pattern(global_options, &tvec, &rvec, l, chain, pattern_mask, avh);
			}
		}
		s_vector container;
		kv_init(container);

		si_hash *index_dct = kh_init(sim);
		int index = 0;
		s_vector cdr3_vec;
		for (khiter_t k = kh_begin(avh); k != kh_end(avh); k++)
		{
			if (kh_exist(avh, k))
			{
				char *mstr = (char*)kh_key(avh, k);
				kv_init(cdr3_vec);
				ic_split_string(mstr, 0, &cdr3_vec);
				for(int i = 0; i < kv_size(cdr3_vec); i++)
				{
					char* cdr3 = kv_A(cdr3_vec, i);
					if(key_exist_si_hash(index_dct, cdr3) == false)
					{
						append_to_si_hash(index_dct, cdr3, index++);
					}
				}
				destruct_s_vector(&cdr3_vec);
				kv_push(char*, container, mstr);
			}
		}
		debug("container size: %d", kv_size(container));
		// debug("step 2");
		si_hash *tunq = kh_init(sim);
		sort_s_vector_by_length(&container);
		ivi_hash *ivih = kh_init(ivim);

		si_hash *punq = kh_init(sim);
		char idx_pair[H1];

		for(int k1 = 0; k1 < kv_size(container); k1++)
		{
			char *mstr = kv_A(container, k1);
			kv_init(cdr3_vec);
			ic_split_string(mstr, 0, &cdr3_vec);
			for(int i = 0; i < kv_size(cdr3_vec); i++)
			{
				char* cdr3 = kv_A(cdr3_vec, i);
				int idx = get_value_si_hash(index_dct, cdr3, 0);
				append_to_ivi_hash(ivih, idx, k1);
			}
			destruct_s_vector(&cdr3_vec);
		}

		for (khiter_t k = kh_begin(ivih); k != kh_end(ivih); k++)
		{
			if (!kh_exist(ivih, k))
			{
				continue;
			}
			i_vector *vec = kh_val(ivih, k);
			for(int ii = 0; ii < kv_size(*vec); ii++)
			{
				int idx_i = kv_A(*vec, ii);
				char *mstr1 = kv_A(container, idx_i);
				i_vector v1;
				kv_init(v1);
				kv_init(cdr3_vec);
				ic_split_string(mstr1, 0, &cdr3_vec);
				for(int i = 0; i < kv_size(cdr3_vec); i++)
				{
					char* cdr3 = kv_A(cdr3_vec, i);
					kv_push(int, v1, get_value_si_hash(index_dct, cdr3, 0));
				}
				destruct_s_vector(&cdr3_vec);

				for(int jj = ii + 1; jj < kv_size(*vec); jj++)
				{
					int idx_j = kv_A(*vec, jj);
					sprintf(idx_pair, "%d %d", idx_i, idx_j);
					if(key_exist_si_hash(punq, idx_pair) == true)
					{
						continue;
					}
					add_to_si_hash(punq, idx_pair);

					char *mstr2 = kv_A(container, idx_j);
					i_vector v2;
					kv_init(v2);
					kv_init(cdr3_vec);
					ic_split_string(mstr2, 0, &cdr3_vec);
					for(int i = 0; i < kv_size(cdr3_vec); i++)
					{
						char* cdr3 = kv_A(cdr3_vec, i);
						int v2_v = get_value_si_hash(index_dct, cdr3, 0);
						kv_push(int, v2, v2_v);
					}
					destruct_s_vector(&cdr3_vec);

					int flag = issubset(&v1, &v2);
					if (flag == 1)
					{
						add_to_si_hash_nd(tunq, mstr2);
					}
					else if (flag == -1)
					{
						add_to_si_hash_nd(tunq, mstr1);
					}
					kv_destroy(v2);
				}
				kv_destroy(v1);
			}
		}
		destruct_ivi_hash(ivih);
		destruct_si_hash(punq);
/*
                for(int k1 = 0; k1 < kv_size(container); k1++)
                {
                        char *mstr1 = kv_A(container, k1);
                        i_vector v1;
                        kv_init(v1);
                        kv_init(cdr3_vec);
                        ic_split_string(mstr1, 0, &cdr3_vec);
                        for(int i = 0; i < kv_size(cdr3_vec); i++)
                        {
                                char* cdr3 = kv_A(cdr3_vec, i);
                                kv_push(int, v1, get_value_si_hash(index_dct, cdr3, 0));
                        }
                        destruct_s_vector(&cdr3_vec);
                        int v1_0 = kv_A(v1, 0);

                        for(int k2 = k1+1; k2 < kv_size(container); k2++)
                        {
                                char *mstr2 = kv_A(container, k2);
                                i_vector v2;
                                kv_init(v2);
                                kv_init(cdr3_vec);
                                ic_split_string(mstr2, 0, &cdr3_vec);
                                bool found = false;
                                for(int i = 0; i < kv_size(cdr3_vec); i++)
                                {
                                        char* cdr3 = kv_A(cdr3_vec, i);
                                        int v2_v = get_value_si_hash(index_dct, cdr3, 0);
                                        kv_push(int, v2, v2_v);
                                        if (v1_0 == v2_v)
                                        {
                                                found = true;
                                        }
                                }
                                destruct_s_vector(&cdr3_vec);

                                if(found == true)
                                {
                                        int flag = issubset(&v1, &v2);
                                        if (flag == 1)
                                        {
                                                add_to_si_hash_nd(tunq, mstr2);
                                        }
                                        else if (flag == -1)
                                        {
                                                add_to_si_hash_nd(tunq, mstr1);
                                        }
                                }

                                kv_destroy(v2);
                        }
                        kv_destroy(v1);
                }
 */
		destruct_si_hash(index_dct);
		kv_destroy(container);
		for (khiter_t k = kh_begin(avh); k != kh_end(avh); k++)
		{
			if (kh_exist(avh, k))
			{
				char *mstr = (char*)kh_key(avh, k);
				if(key_exist_si_hash(tunq, mstr) == true)
				{
					debug("%s removed", mstr);
					continue;
				}
				all_vector *avec = kh_val(avh, k);
				sort_all_vector(avec);
				all_t *pat = kv_A(*avec, 0);
				fprintf(kmer_fp, "%s %.1e %.2f", pat->pattern, pat->p_value, pat->entropy_fraction);
				kv_init(fields);
				ic_split_string(mstr, 0, &fields);
				for(int j0 = 0; j0 < kv_size(fields); j0++)
				{
					char *ucdr3 = kv_A(fields, j0);
					fprintf(kmer_fp, " %s", ucdr3);
				}
				destruct_s_vector(&fields);
				fprintf(kmer_fp, "\n");
			}
		}
		destruct_si_hash_nd(tunq);
		kv_destroy(tvec);
		kv_destroy(rvec);
		destruct_li_hash(pdict);
		destruct_av_hash(avh);
	}
	fclose(kmer_fp);
	sort_cluster_on_entfrc(kmer_file);
}

void do_all_clustering(parameter_t* global_options)
{
	if (global_options->unique_trb_target_total > 0 && global_options->unique_trb_refer_total > 0)
	{
		discover_enriched_pattern(global_options, CHAIN_BD);
		output_hla_association(global_options, CHAIN_BD);
		output_cluster(global_options, CHAIN_BD);
	}

	if (global_options->unique_tra_target_total > 0 && global_options->unique_tra_refer_total > 0)
	{
		discover_enriched_pattern(global_options, CHAIN_AG);
		output_hla_association(global_options, CHAIN_AG);
		output_cluster(global_options, CHAIN_AG);
	}
}
