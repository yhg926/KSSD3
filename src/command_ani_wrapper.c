#include "command_ani.h"
// #include "command_sketch.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <argp.h>
#include <argz.h>
#include <err.h>
#include <errno.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>

/*** argp wrapper ***/
struct arg_ani
{
	struct arg_global *global;

	char *name;
};

static struct argp_option opt_ani[] =
	{
		{"ref", 'r', "<DIR>", 0, "Path of reference sketches, do not set if need trianlge\v", 1},
		{"query", 'q', "<DIR>", 0, "Path of query sketches \v", 1},
		{"index", 'i', "<FILE>", 0, "Inverted indexing of combined sketch.\v", 2},
		//		{"model", 'M', "<FILE>", 0, "specify the trained ani model.\v", 2},
		{"naive", 'v', 0, 0, "Use naive distane / ani caculation.\v", 2},
		{"outfmt", 'm', "<0/1/2>", 0, "print results detail/matrix/triangle (0/1/2),[0]\v", 3},
		{"slmetrics", 's', "<+-1..5>", 0, "Select metrics: Standard(1) / MashD(2) / AafD(3) / MashD_if_far(4) / AafD_if_far(5), ANI(+) and distance(-) [1]\v", 3},
		{"afcut", 'f', "<FLOAT>", 0, "When report, Skip alignment fraction < [0.1] \v", 3},
		{"anicut", 'n', "<FLOAT>", 0, "When report, Skip ani < [0.5] \v", 3},
		{"control", 'c', "<FLOAT>", 0, "Skip duplicated samples (distance < c) [0] \v", 3},
		{"glist", 'g', "<FILE>", 0, "Sample outfile for kssd set grouping \v", 4},
		{"outfile", 'o', "<FILE>", 0, "outfile path [STDOUT]\v", 5},
		{"threads", 'p', "<INT>", 0, "Threads number to use \v", 6},
		{"diagonal", 'd', 0, 0, "set diagonal\v", 7},
		{"exception", 'e', "<INT>", 0, "set distance value when skipped [0]\v", 8},
		{"pair", 777, 0, 0, "Compute ANI directly from genomes (.fna/.fq) \v", 8},
		{0}};

static char doc_ani[] =
	""
	""
	"\v"
	"";

ani_opt_t ani_opt = {
	.fmt = 0, // 0:detail, 1:matrix 2: triangle
	.c = 0.0, // control duplicated sample by skip distance < c;
	.p = 1,
	.d = 0,	   // diagonal
	.v = 0,	   // naive model
	.s = 1,	   // select metrics: 1:standard, 2:MashD, 3:AafD, 4:MashD_if_far, 5:AafD_if_far
	.pair = 0, // pairwise compute
	.afcut = 0.1,
	.anicut = 0.3,
	.e = 1, // abort
	.index[0] = '\0',
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.outf[0] = '\0',
	.gl[0] = '\0',
	//	.model[0] = '\0',
	.num_remaining_args = 0, // int num_remaining_args; no option arguments num.
	.remaining_args = NULL	 // char **remaining_args; no option arguments array.
};

static error_t parse_ani(int key, char *arg, struct argp_state *state)
{
	struct arg_ani *ani = state->input;

	switch (key)
	{
	case 'm':
	{
		ani_opt.fmt = atoi(arg);
		break;
	}
	case 'c':
	{
		if (!isdigit(*arg))
		{
			printf("-c requires a positive float numer\n");
			exit(1);
		}
		ani_opt.c = atof(arg);
		break;
	}
	case 'e':
	{
		int v = atof(arg);
		if (v < 1)
		{
			printf("-e requires a positive numer >= 1 \n");
			exit(1);
		}
		ani_opt.e = v;
		break;
	}
	case 's':
	{
		ani_opt.s = atoi(arg);
		break;
	}
	case 'p':
	{
		ani_opt.p = atoi(arg);
		break;
	}
	case 'f':
	{
		ani_opt.afcut = atof(arg);
		break;
	}
	case 'n':
	{
		ani_opt.anicut = atof(arg);
		break;
	}
	case 'i':
	{
		strcpy(ani_opt.index, arg);
		break;
	}
	case 'M':
	{
		strcpy(ani_opt.model, arg);
		break;
	}
	case 'q':
	{
		strcpy(ani_opt.qrydir, arg);
		break;
	}
	case 'r':
	{
		strcpy(ani_opt.refdir, arg);
		break;
	}
	case 'o':
	{
		strcpy(ani_opt.outf, arg);
		break;
	}
	case 'g':
	{
		strcpy(ani_opt.gl, arg);
		break;
	}
	case 'd':
	{
		ani_opt.d = 1;
		break;
	}
	case 'v':
	{
		ani_opt.v = 1;
		break;
	}
	case 777:
	{
		ani_opt.pair = 1;
		break;
	}
	case ARGP_KEY_ARGS:
	{
		ani_opt.num_remaining_args = state->argc - state->next;
		ani_opt.remaining_args = state->argv + state->next;
		break;
	}
	case ARGP_KEY_END:
	{

		if (!ani_opt.pair && ani_opt.qrydir[0] == '\0')
		{
			printf("\nError: Mandatory options: '-q' are missing unless use mode '--pair'.\n\n");
			argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
			argp_usage(state);
		}
		if (ani_opt.s < -5 || ani_opt.s > 5 || ani_opt.s == 0)
			printf("\nError: -s option should be within range 1..5 or -5..-1\n\n");

		break;
		/*
					if(state->argc<2)
					{
				printf("\v");
						argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
						printf("\v");
				argp_state_help(state,stdout,ARGP_HELP_LONG);
				printf("\v");
				return EINVAL;
					}
		*/
	}
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp_ani =
	{
		opt_ani,
		parse_ani,
		0, //  "[arguments ...]",
		doc_ani};
extern const char sorted_comb_ctxgid64obj32[];

get_generic_dist_from_features_fn get_generic_dist_from_features = NULL;

extern simple_sketch_t *simple_genomes2mem2sortedctxobj64_mem(infile_tab_t *infile_stat, int drfold);

int cmd_ani(struct argp_state *state)
{

	struct arg_ani ani = {
		0,
	};
	int argc = state->argc - state->next + 1;
	char **argv = &state->argv[state->next - 1];
	ani.global = state->input;
	argp_parse(&argp_ani, argc, argv, ARGP_IN_ORDER, &argc, &ani);
	state->next += argc - 1;
	size_t file_size;

		// instantilize distance fn according selection of naive model or not
	get_generic_dist_from_features = ani_opt.v ? get_naive_dist : lm3ways_dist_from_features;

	if (ani_opt.qrydir[0] != '\0')
	{
		dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt.qrydir, sketch_stat), &file_size);
		if (qry_dim_sketch_stat->conflict) //imply raw reads? which not trained in lm3ways_dist_from_features.  
			get_generic_dist_from_features = get_naive_dist; 

		const_comask_init(qry_dim_sketch_stat);
		if (ani_opt.refdir[0] == '\0')
			return 1; // compute_triangle(&ani_opt);
		else
		{
			if (file_exists_in_folder(ani_opt.refdir, sorted_comb_ctxgid64obj32))
			{
				return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(&ani_opt); // compute_ani(&ani_opt);
			}
			else
			{
				comb_sortedsketch64Xcomb_sortedsketch64_sorted_per_q(&ani_opt);
				//comb_sortedsketch64Xcomb_sortedsketch64(&ani_opt);
				return 1;
			}
			// else err(EXIT_FAILURE, "%s(): Failed to detect index file '%s/%s'\nrun kssd3 sketch -i <sketch_folder> to create one", __func__, ani_opt.refdir, sorted_comb_ctxgid64obj32);
		}
	}
	else if (ani_opt.num_remaining_args >= 2)
	{
		// initializing
		dim_sketch_stat_t lco_stat_val;
		lco_stat_val.klen = 3 * NUM_CODENS + 1; // klen is 3*NUM_CODENS+1
		lco_stat_val.coden_len = NUM_CODENS;
		const_comask_init(&lco_stat_val);

		infile_tab_t *genomes_infiletab = organize_infile_frm_arg(ani_opt.num_remaining_args, (&ani_opt)->remaining_args, 1);

		simple_sketch_t *simple_sketch = simple_genomes2mem2sortedctxobj64_mem(genomes_infiletab, 8);
		simple_sortedsketch64Xcomb_sortedsketch64(simple_sketch, genomes_infiletab, (&ani_opt));
	}
}
