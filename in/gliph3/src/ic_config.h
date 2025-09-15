#ifndef _CONFIG_H_
#define _CONFIG_H_

#include "ic_common.h"
#include "ic_type.h"

typedef struct {
	// input parameters
	char out_prefix[MAX_FILENAME_LENGTH];
	// int overwrite_output;

	// doing HLA association testing, hla file is required
	char hla_file[MAX_FILENAME_LENGTH];
	int number_of_hla_field;
	double hla_association_pvalue_cutoff;

	// number of aas from ends ignored during calculation
	int ag_ignored_v_end;
	int ag_ignored_j_end;
	int bd_ignored_v_end;
	int bd_ignored_j_end;

	// doing cluster based on motifs, reference file is required
	char ag_refer_file[MAX_FILENAME_LENGTH];
	char bd_refer_file[MAX_FILENAME_LENGTH];

	// data cleaning
	char cdr3_file[MAX_FILENAME_LENGTH];
	int ag_cdr3_length_min_cutoff;
	int ag_cdr3_length_max_cutoff;
	int bd_cdr3_length_min_cutoff;
	int bd_cdr3_length_max_cutoff;

	int min_motif_length;
	int max_motif_length;
	int max_diff_position;
	int pattern_unique_sample_cutoff;
	int pattern_ove_cutoff;
	float pattern_pvalue_cutoff;
	int same_v;
	int purge_cluster;
	float entropy_fraction_cutoff;
	float purge_fraction;
	int continuous_motif;

	int plot;
	int overlap_cutoff;

	// dynamic parameters, not provided by user
	bool hla_availability;
	int unique_tra_target_total;
	int unique_trb_target_total;
	int unique_tra_refer_total;
	int unique_trb_refer_total;

	s_vector *tmpfile_vector;
	t_vector *target_cdr3_vector;
	c_vector *ag_refer_cdr3_vector;
	c_vector *bd_refer_cdr3_vector;
	svs_hash *hla_hash;
} parameter_t;

void set_default_parameter(parameter_t* parameter);
void output_parameter(parameter_t* parameter);
void check_parameter(parameter_t *parameter);
void read_config_file(const char* config_file, parameter_t* parameter);


#endif
