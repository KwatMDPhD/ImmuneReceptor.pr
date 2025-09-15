#include "ic_common.h"
#include "ic_config.h"
#include "ic_all.h"
#include "ic_util.h"
#include "ic_type.h"
#include <getopt.h>

parameter_t global_options;

void showspin(void)     // Shows a different charcter every time its called
{
	// static char *spin="-/|\\";
	static char *spin=".........+";
	static int pos=-1;
	if(!spin[++pos]) pos=0; // If we hit the NULL terminator, skip to beginning of string
	// We have to use stderr, so it doesn't buffer, since we're not printing a newline
	// fprintf(stderr,"%c",spin[pos]); // Show character at position 'pos' in string
	fprintf(stderr, "%c", spin[pos]);
}

void dospin()    // Starts showing spinner repeatedly
{
	showspin();
	signal(SIGALRM, dospin); // Register this function to be called on SIGALRM
	alarm(60); // Cause SIGALRM to happen 1 second later
}

void endspin()           // Cancels the spinner
{
	alarm(0);
	fprintf(stderr,"\n"); // Erase the spinner
}

void startspin(const char* indicator)
{
	fprintf(stderr, "%s ",indicator);
	dospin();
}

static void usage(){
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: ictools (Tools for clustering CDR3 sequence data)\n");
	fprintf(stderr, "Version: %s\n", VERSION);
	fprintf(stderr,
	        "\nUsage:   ictools -c parameter_file\n"
	        "Default parameters are as following:\n\n"
	        );
	fprintf(stderr,
	        "\tout_prefix                    = <required>\n"
	        "\thla_file                      = [optional]\n"
	        "\tag_refer_file                 = [optional]\n"
	        "\tbd_refer_file                 = [optional]\n"
	        "\tcdr3_file                     = <required>\n"
	        "\tnumber_of_hla_field           = 1     \n"
	        "\thla_association_pvalue_cutoff = 0.001 \n"
	        "\tag_ignored_v_end              = 3     \n"
	        "\tag_ignored_j_end              = 3     \n"
	        "\tbd_ignored_v_end              = 3     \n"
	        "\tbd_ignored_j_end              = 3     \n"
	        "\tpattern_unique_sample_cutoff  = 2     \n"
	        "\tmin_motif_length              = 3     \n"
	        "\tmax_motif_length              = 5     \n"
	        "\tmax_diff_position             = 2     \n"
	        "\tpattern_ove_cutoff            = 10    \n"
	        "\tpattern_pvalue_cutoff         = 0.0001\n"
	        "\tsame_v                        = 1     \n"
	        "\tpurge_cluster                 = 1     \n"
	        "\tcontinuous_motif              = 1     \n"
	        "\tentropy_fraction_cutoff       = 0.9   \n"
	        "\tpurge_fraction                = 0.95  \n"
	        "\toverlap_cutoff                = 60    \n"
	        "\tag_cdr3_length_min_cutoff     = 8     \n"
	        "\tag_cdr3_length_max_cutoff     = 30    \n"
	        "\tbd_cdr3_length_min_cutoff     = 8     \n"
	        "\tbd_cdr3_length_max_cutoff     = 30    \n"
	        );
	fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{
	char *config_file = 0;
	int n = 0;
	while ((n = getopt(argc, argv, "c:h")) >= 0) {
		switch (n) {
		case 'h': break;
		case 'c': config_file = strdup(optarg); break;
		}
	}
	if (!config_file)
	{
		usage();
		return EXIT_FAILURE;
	}

	set_default_parameter(&global_options);
	read_config_file(config_file, &global_options);
	check_parameter(&global_options);

	startspin("Preparing data\t\t\t");
	prepare_data(&global_options);
	endspin();

	startspin("Working on pattern\t\t");
	do_all_clustering(&global_options);
	endspin();

	output_parameter(&global_options);

	char cmd[H5];
	sprintf(cmd, "rm %s", global_options.cdr3_file);
	system_call(cmd);
	return(0);
}
