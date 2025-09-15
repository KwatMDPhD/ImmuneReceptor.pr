#include "ic_config.h"
#include "ic_util.h"
#include "ic_type.h"
#include "ic_io.h"

void set_default_parameter(parameter_t* parameter)
{
	parameter->out_prefix[0]                         = '\0';
	parameter->hla_file[0]                           = '\0';
	parameter->ag_refer_file[0]                      = '\0';
	parameter->bd_refer_file[0]                      = '\0';
	parameter->cdr3_file[0]                          = '\0';

	parameter->number_of_hla_field                   = 1;
	parameter->hla_association_pvalue_cutoff         = 0.001;
	parameter->ag_ignored_v_end                      = 3;
	parameter->ag_ignored_j_end                      = 3;
	parameter->bd_ignored_v_end                      = 3;
	parameter->bd_ignored_j_end                      = 3;
	parameter->pattern_unique_sample_cutoff          = 2;

	parameter->min_motif_length                      = 3;
	parameter->max_motif_length                      = 5;
	parameter->max_diff_position                     = 2;
	parameter->pattern_ove_cutoff                    = 10;
	parameter->pattern_pvalue_cutoff                 = 0.0001;
	parameter->same_v                                = 1;
	parameter->purge_cluster                         = 1;
	parameter->continuous_motif                      = 1;
	parameter->entropy_fraction_cutoff               = 0.9;
	parameter->purge_fraction                        = 0.95;
	parameter->overlap_cutoff                        = 60;

	parameter->ag_cdr3_length_min_cutoff             = 8;
	parameter->ag_cdr3_length_max_cutoff             = 30;
	parameter->bd_cdr3_length_min_cutoff             = 8;
	parameter->bd_cdr3_length_max_cutoff             = 30;
	parameter->plot                                  = 0;
	parameter->tmpfile_vector                        = icalloc(1, s_vector);
	kv_init(*(parameter->tmpfile_vector));
}

void output_parameter(parameter_t* parameter)
{
	char filename[MAX_FILENAME_LENGTH];
	sprintf(filename, "%s_parameter.txt", parameter->out_prefix);
	FILE *config_fp = fopen(filename, "w");
	if (config_fp == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", filename);
		exit(EXIT_FAILURE);
	}

	fprintf(config_fp, "out_prefix                    = %s\n", parameter->out_prefix                    );
	fprintf(config_fp, "hla_file                      = %s\n", parameter->hla_file                      );
	fprintf(config_fp, "ag_refer_file                 = %s\n", parameter->ag_refer_file                 );
	fprintf(config_fp, "bd_refer_file                 = %s\n", parameter->bd_refer_file                 );
	fprintf(config_fp, "cdr3_file                     = %s\n", parameter->cdr3_file                     );
	fprintf(config_fp, "number_of_hla_field           = %d\n", parameter->number_of_hla_field           );
	fprintf(config_fp, "hla_association_pvalue_cutoff = %f\n", parameter->hla_association_pvalue_cutoff );
	fprintf(config_fp, "ag_ignored_v_end              = %d\n", parameter->ag_ignored_v_end              );
	fprintf(config_fp, "ag_ignored_j_end              = %d\n", parameter->ag_ignored_j_end              );
	fprintf(config_fp, "bd_ignored_v_end              = %d\n", parameter->bd_ignored_v_end              );
	fprintf(config_fp, "bd_ignored_j_end              = %d\n", parameter->bd_ignored_j_end              );
	fprintf(config_fp, "min_motif_length              = %d\n", parameter->min_motif_length              );
	fprintf(config_fp, "max_motif_length              = %d\n", parameter->max_motif_length              );
	fprintf(config_fp, "max_diff_position             = %d\n", parameter->max_diff_position             );
	fprintf(config_fp, "pattern_unique_sample_cutoff  = %d\n", parameter->pattern_unique_sample_cutoff  );
	fprintf(config_fp, "pattern_ove_cutoff            = %d\n", parameter->pattern_ove_cutoff            );
	fprintf(config_fp, "pattern_pvalue_cutoff         = %.1e\n",parameter->pattern_pvalue_cutoff        );
	fprintf(config_fp, "same_v                        = %d\n", parameter->same_v                        );
	fprintf(config_fp, "continuous_motif              = %d\n", parameter->continuous_motif              );
	fprintf(config_fp, "purge_cluster                 = %d\n", parameter->purge_cluster                 );
	fprintf(config_fp, "purge_fraction                = %.2f\n", parameter->purge_fraction              );
	fprintf(config_fp, "entropy_fraction_cutoff       = %.2f\n", parameter->entropy_fraction_cutoff     );
	fprintf(config_fp, "ag_cdr3_length_min_cutoff     = %d\n", parameter->ag_cdr3_length_min_cutoff     );
	fprintf(config_fp, "ag_cdr3_length_max_cutoff     = %d\n", parameter->ag_cdr3_length_max_cutoff     );
	fprintf(config_fp, "bd_cdr3_length_min_cutoff     = %d\n", parameter->bd_cdr3_length_min_cutoff     );
	fprintf(config_fp, "bd_cdr3_length_max_cutoff     = %d\n", parameter->bd_cdr3_length_max_cutoff     );
	fprintf(config_fp, "plot                          = %d\n", parameter->plot                          );
	fprintf(config_fp, "overlap_cutoff                = %d\n", parameter->overlap_cutoff                );

	fclose(config_fp);
}

void check_parameter(parameter_t *parameter)
{
	if (strlen(parameter->out_prefix) == 0)
	{
		fprintf(stderr, "out_prefix is not set!\n");
		exit(EXIT_FAILURE);
	}
	if (strlen(parameter->cdr3_file) == 0 || access(parameter->cdr3_file, R_OK) != 0)
	{
		fprintf(stderr, "cdr3_file is not set or not readable!\n");
		exit(EXIT_FAILURE);
	}
	if (strlen(parameter->ag_refer_file) != 0 && access(parameter->ag_refer_file, R_OK) != 0)
	{
		fprintf(stderr, "ag_refer_file is set but not readable!\n");
		exit(EXIT_FAILURE);
	}

	if (strlen(parameter->bd_refer_file) != 0 && access(parameter->bd_refer_file, R_OK) != 0)
	{
		fprintf(stderr, "bd_refer_file is set but not readable!\n");
		exit(EXIT_FAILURE);
	}

	if (strlen(parameter->hla_file) != 0 && access(parameter->hla_file, R_OK) != 0)
	{
		fprintf(stderr, "ag_refer_file is set but not readable!\n");
		exit(EXIT_FAILURE);
	}

	if(parameter->ag_cdr3_length_min_cutoff <= parameter->ag_ignored_j_end + parameter->ag_ignored_v_end)
	{
		fprintf(stderr, "ag_cdr3_length_min_cutoff(%d) needs to be greater than the sum of ag_ignored_v_end(%d) and ag_ignored_j_end(%d)!\n",
		        parameter->ag_cdr3_length_min_cutoff, parameter->ag_ignored_v_end, parameter->ag_ignored_j_end
		        );
		exit(EXIT_FAILURE);
	}

	if(parameter->bd_cdr3_length_min_cutoff <= parameter->bd_ignored_j_end + parameter->bd_ignored_v_end)
	{
		fprintf(stderr, "bd_cdr3_length_min_cutoff(%d) needs to be greater than the sum of bd_ignored_v_end(%d) and bd_ignored_j_end(%d)!\n",
		        parameter->bd_cdr3_length_min_cutoff, parameter->bd_ignored_v_end, parameter->bd_ignored_j_end
		        );
		exit(EXIT_FAILURE);
	}
}

void read_config_file(const char *config_file, parameter_t *parameter)
{
	s_vector larray, fields;
	kv_init(larray);
	ic_read_into_line_vector(config_file, &larray);

	for (int lnum = 0; lnum < kv_size(larray); lnum++)
	{
		char *cline = kv_A(larray, lnum);
		for (int i = 0; i < strlen(cline); i++)
		{
			if (cline[i] == '#')
			{
				cline[i] = '\0';
				break;
			}
		}

		cline = ic_strip(cline);

		kv_init(fields);
		ic_split_string(cline, '=', &fields);
		if (kv_size(fields) == 2)
		{
			char *ckey = ic_strtrimr((char*)kv_A(fields, 0));
			char *cval = ic_strtriml((char*)kv_A(fields, 1));

			if (strcmp(ckey, "out_prefix") == 0)
			{
				sprintf(parameter->out_prefix, "%s", cval);
			}

			if (strcmp(ckey, "ag_refer_file") == 0)
			{
				sprintf(parameter->ag_refer_file, "%s", cval);
			}

			if (strcmp(ckey, "bd_refer_file") == 0)
			{
				sprintf(parameter->bd_refer_file, "%s", cval);
			}
			// process HLA options
			if (strcmp(ckey, "hla_file") == 0)
			{
				sprintf(parameter->hla_file, "%s", cval);
			}

			if (strcmp(ckey, "number_of_hla_field") == 0)
			{
				parameter->number_of_hla_field = atoi(cval);
			}

			if (strcmp(ckey, "plot") == 0)
			{
				if(cval[0] == 'N' || cval[0] == 'n' || cval[0] == '0')
				{
					parameter->plot = 0;
				}
				else
				{
					parameter->plot = 1;
				}
			}

			if (strcmp(ckey, "hla_association_pvalue_cutoff") == 0)
			{
				parameter->hla_association_pvalue_cutoff = atof(cval);
			}

			if (strcmp(ckey, "max_diff_position") == 0)
			{
				parameter->max_diff_position = atoi(cval);
			}

			if (strcmp(ckey, "min_motif_length") == 0)
			{
				parameter->min_motif_length = atoi(cval);
			}

			if (strcmp(ckey, "max_motif_length") == 0)
			{
				parameter->max_motif_length = atoi(cval);
			}

			if (strcmp(ckey, "pattern_unique_sample_cutoff") == 0)
			{
				parameter->pattern_unique_sample_cutoff = atoi(cval);
			}

			if (strcmp(ckey, "pattern_ove_cutoff") == 0)
			{
				parameter->pattern_ove_cutoff = atoi(cval);
			}

			if (strcmp(ckey, "pattern_pvalue_cutoff") == 0)
			{
				parameter->pattern_pvalue_cutoff = atof(cval);
			}

			if (strcmp(ckey, "overlap_cutoff") == 0)
			{
				parameter->overlap_cutoff = atoi(cval);
			}

			if (strcmp(ckey, "entropy_fraction_cutoff") == 0)
			{
				parameter->entropy_fraction_cutoff = atof(cval);
			}

			if (strcmp(ckey, "purge_fraction") == 0)
			{
				parameter->purge_fraction = atof(cval);
			}

			if (strcmp(ckey, "same_v") == 0)
			{
				if(cval[0] == 'N' || cval[0] == 'n' || cval[0] == '0')
				{
					parameter->same_v = 0;
				}
				else
				{
					parameter->same_v = 1;
				}
			}

			if (strcmp(ckey, "continuous_motif") == 0)
			{
				if(cval[0] == 'N' || cval[0] == 'n' || cval[0] == '0')
				{
					parameter->continuous_motif = 0;
				}
				else
				{
					parameter->continuous_motif = 1;
				}
			}

			if (strcmp(ckey, "purge_cluster") == 0)
			{
				if(cval[0] == 'N' || cval[0] == 'n' || cval[0] == '0')
				{
					parameter->purge_cluster = 0;
				}
				else
				{
					parameter->purge_cluster = 1;
				}
			}

			// process cdr3 parameters
			if (strcmp(ckey, "cdr3_file") == 0)
			{
				sprintf(parameter->cdr3_file, "%s", cval);
			}

			if (strcmp(ckey, "ag_cdr3_length_min_cutoff") == 0)
			{
				parameter->ag_cdr3_length_min_cutoff = atoi(cval);
			}

			if (strcmp(ckey, "ag_cdr3_length_max_cutoff") == 0)
			{
				parameter->ag_cdr3_length_max_cutoff = atoi(cval);
			}

			if (strcmp(ckey, "bd_cdr3_length_min_cutoff") == 0)
			{
				parameter->bd_cdr3_length_min_cutoff = atoi(cval);
			}

			if (strcmp(ckey, "bd_cdr3_length_max_cutoff") == 0)
			{
				parameter->bd_cdr3_length_max_cutoff = atoi(cval);
			}

			if (strcmp(ckey, "ag_ignored_v_end") == 0)
			{
				parameter->ag_ignored_v_end = atoi(cval);
			}

			if (strcmp(ckey, "ag_ignored_j_end") == 0)
			{
				parameter->ag_ignored_j_end = atoi(cval);
			}

			if (strcmp(ckey, "bd_ignored_v_end") == 0)
			{
				parameter->bd_ignored_v_end = atoi(cval);
			}

			if (strcmp(ckey, "bd_ignored_j_end") == 0)
			{
				parameter->bd_ignored_j_end = atoi(cval);
			}

		}
		destruct_s_vector(&fields);
	}
	destruct_s_vector(&larray);
}
