#include "command_ani.h"
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
		{"ref", 'r', "<DIR>", 0, "Path of reference sketches, do not set if need trianlge\v", 9},
		{"query", 'q', "<DIR>", 0, "Path of query sketches \v", 1},
		{"index", 'i', "<FILE>", 0, "Inverted indexing of combined sketch.\v", 2},
		{"model", 'M', "<FILE>", 0, "specify the trained ani model.\v", 2},
		{"outfmt", 'm', "<0/1/2>", 0, "print results detail/matrix/triangle (0/1/2) [0]\v", 2},
		{"afcut", 'f', "<FLOAT>", 0, "When report, Skip alignment fraction < [0.1] \v", 3},
		{"anicut", 'n', "<FLOAT>", 0, "When report, Skip ani < [0.5] \v", 3},
		{"control", 'c', "<FLOAT>", 0, "Skip duplicated samples (distance < c) [0] \v", 3},
		{"glist", 'g', "<FILE>", 0, "Sample outfile for kssd set grouping \v", 4},
		{"outfile", 'o', "<FILE>", 0, "outfile path [STDOUT]\v", 5},
		{"threads", 'p', "<INT>", 0, "Threads number to use \v", 6},
		{"diagonal", 'd', 0, 0, "set diagonal\v", 7},
		{"exception", 'e', "<INT>", 0, "set distance value when XnY == 0 \v", 8},
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
	.d = 0, // diagonal
	.afcut = 0.1,
	.anicut = 0.3,
	.e = 0, // abort
	.index[0] = '\0',
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.outf[0] = '\0',
	.gl[0] = '\0',
	.model[0] = '\0',
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
	case ARGP_KEY_ARGS:
	{
		ani_opt.num_remaining_args = state->argc - state->next;
		ani_opt.remaining_args = state->argv + state->next;
		break;
	}
	case ARGP_KEY_END:
	{
		if (ani_opt.qrydir[0] == '\0')
		{
			printf("\nError: Mandatory options: '-q' are missing.\n\n");
			argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
			argp_usage(state);
		}
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

	if (ani_opt.qrydir[0] != '\0')
	{
		dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt.qrydir, sketch_stat), &file_size);
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
				err(EXIT_FAILURE, "%s(): Failed to detect index file '%s/%s'\nrun kssd3 sketch -i <sketch_folder> to create one", __func__, ani_opt.refdir, sorted_comb_ctxgid64obj32);
		}
	}
}
