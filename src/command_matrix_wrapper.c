#include "command_matrix.h"
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
struct arg_matrix
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_matrix[] =
{
	{"ref",'r',"<DIR>", 0, "Path of reference sketches, do not set if need trianlge\v",9},
	{"query",'q',"<DIR>", 0, "Path of query sketches \v",1},
	{"metric",'m',"<0/1>", 0, "Using mashD or aafD (0/1) [0]\v",2},
	{"control",'c',"<FLOAT>",0,"Skip duplicated samples (distance < c) [0] \v",3},
	{"glist",'g',"<FILE>",0,"Sample outfile for kssd set grouping \v",4},
	{"outfile",'o',"<FILE>",0,"Matrix outfile path [STDOUT]\v",5},
	{"threads",'p',"<INT>", 0, "Threads number to use \v",6},
	{"diagonal",'d',0, 0, "set diagonal\v",7},
	{"exception",'e',"<INT>", 0, "set distance value when XnY == 0 \v",8},
  { 0 }
};

static char doc_matrix[] =
  ""
  ""
  "\v"
  ""
  ;

matrix_opt_t matrix_opt ={
	.metric = 0, // 0:mashD, 1:aafD
	.c = 0.0, //control duplicated sample by skip distance < c;
	.p = 1,
	.d = 0, //diagonal
	.e = -1, //abort
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.outf[0] = '\0',
	.gl[0]= '\0',
	.num_remaining_args = 0, //int num_remaining_args; no option arguments num.
	.remaining_args = NULL //char **remaining_args; no option arguments array.
};

static error_t parse_matrix(int key, char* arg, struct argp_state* state) {
  struct arg_matrix* matrix = state->input;
	
  switch(key)
  {		
		case 'm':
		{
				matrix_opt.metric = atoi(arg);
				break;
		}
		case 'c':
		{
			if(!isdigit(*arg)) {
				printf("-c requires a positive float numer\n");			 
				exit(1);
			}
			matrix_opt.c = atof(arg);
			break;
		}
		case 'e':
		{
			int v = atof(arg);
			if( v < 1 ) {
        printf("-e requires a positive numer >= 1 \n");
        exit(1);
      }
      matrix_opt.e = v;
			break;
		}
		case 'p':
		{
			matrix_opt.p = atoi(arg);
			break;
		}
		case 'q':
		{
			strcpy(matrix_opt.qrydir, arg);
			break;
		}
		case 'r':
    {
      strcpy(matrix_opt.refdir, arg);
      break;
    }
		case 'o':
		{	
			strcpy(matrix_opt.outf, arg);
			break;
		}
		case 'g':
		{
			strcpy(matrix_opt.gl, arg);
			break;
		}
		case 'd':
		{
			matrix_opt.d = 1 ;
			break;
		}
		case ARGP_KEY_ARGS:
		{
			matrix_opt.num_remaining_args = state->argc - state->next;
			matrix_opt.remaining_args  = state->argv + state->next;
			break;
		}
    case ARGP_KEY_END:
    {	
			if(matrix_opt.qrydir[0] == '\0'){
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

static struct argp argp_matrix =
{
  opt_matrix,
  parse_matrix,
	0,//  "[arguments ...]",
  doc_matrix
};

int cmd_matrix(struct argp_state* state)
{
  struct arg_matrix matrix = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  matrix.global = state->input;
  argp_parse(&argp_matrix, argc, argv, ARGP_IN_ORDER, &argc, &matrix);

  state->next += argc - 1;
  if (matrix_opt.qrydir[0] != '\0') {
  	if (matrix_opt.refdir[0] == '\0')
	  	return compute_triangle(&matrix_opt);
		else
			return compute_matrix(&matrix_opt);
	}
  return 1;
}
