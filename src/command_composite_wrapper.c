#include "command_composite.h"
#include "global_basic.h"
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

/*** argp wrapper ***/
struct arg_composite
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_composite[] =
{
	{"ref",'r',"<DIR>", 0, "Path of species specific pan uniq kmer database (reference) \v",1},
	{"query",'q',"<DIR>", 0, "Path of query sketches with abundances \v",2},
	{"outfile",'o',"<DIR>",0,"Output path \v",3},
	{"threads",'p',"<INT>", 0, "Threads number to use \v",4},
	{"binVec",'b',0,0,"Output species abundances in Binary Vector format (.abv) \v",5},
	{"idxbv",'i',0,0,"build index of abundance Binary Vector \v",6},
	{"search",'s',"<0-2>",0,"search for similar abundance Binary Vectors using L1norm(1)/L2norm(2)/cosine(0)\v",7},
	{"readabv",'d',0,0,"Read .abv file \v",8},
  { 0 }
};

static char doc_composite[] =
  "\n"
  "The composite doc prefix."
  "\v"
  "The composite doc suffix."
  ;

composite_opt_t composite_opt ={
	.b = 0, //write out abundance binary vector (1) or not (0)
	.i = 0, // index abundance binary vectors (1) or not (0)
	.s = -1, 	
	.d = 0, // read .abv file
	.p = 1,
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.outdir = "./",
	.num_remaining_args = 0, //int num_remaining_args; no option arguments num.
	.remaining_args = NULL //char **remaining_args; no option arguments array.
};

static error_t parse_composite(int key, char* arg, struct argp_state* state) {
  struct arg_composite* composite = state->input;
  assert( composite );
  assert( composite->global );
	
  switch(key)
  {
		case 'b':
		{
			composite_opt.b = 1;
			break;
		}
		case 'i':
		{
			 composite_opt.i = 1;
				break;
		}
		case 's':
		{
				composite_opt.s = atoi(arg);
				break;
		}
		case 'd':
		{
			composite_opt.d = 1;
			break;
		}
		case 'p':
		{
			composite_opt.p = atoi(arg);
			break;
		}
		case 'r':
		{
			strcpy(composite_opt.refdir, arg);
			break;
		}
		case 'q':
		{
			strcpy(composite_opt.qrydir, arg);
			break;
		}
		case 'o':
		{
			strcpy(composite_opt.outdir, arg);
			break;
		}
		case ARGP_KEY_ARGS:
		{
			composite_opt.num_remaining_args = state->argc - state->next;
			composite_opt.remaining_args  = state->argv + state->next;
			break;
		}
    case ARGP_KEY_NO_ARGS:
    {	
			if(state->argc<2)
			{
      	printf("\v");
				argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
				printf("\v");
      	argp_state_help(state,stdout,ARGP_HELP_LONG);
      	printf("\v");
      	return EINVAL;
			}
    }
		break;
    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp_composite =
{
  opt_composite,
  parse_composite,
	0,//  "[arguments ...]",
  doc_composite
};

int cmd_composite(struct argp_state* state)
{
  struct arg_composite composite = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  composite.global = state->input;
	
  argp_parse(&argp_composite, argc, argv, ARGP_IN_ORDER, &argc, &composite);
	
  state->next += argc - 1;
	if( composite_opt.refdir[0] != '\0' ){
		if (composite_opt.qrydir[0] != '\0') //if(argc >1)
			return get_species_abundance (&composite_opt);
		else if (composite_opt.i) 
			return index_abv (&composite_opt);
		else if (composite_opt.s != -1){
			if( composite_opt.s >=0 && composite_opt.s <3 && composite_opt.num_remaining_args >0 )	return abv_search(&composite_opt);
			else printf("\vUsage: %s composite -r <ref> -s <0|1|2> <query.abv>\n\v", state->name);
		}
		else printf("\vUsage: %s composite -r <ref> < mode: -q | -i | -s >\n\v",state->name);
	}
	else if(composite_opt.d){
		if(composite_opt.num_remaining_args < 1)	printf("\vUsage: %s composite -d <query.abv>\n\v", state->name);
		else return read_abv (&composite_opt);	
	}
	else
		printf("\vUsage: %s composite -r <ref> < mode: -q | -i | -s >\n\v", state->name);	

	return 1;
}//cmd_composite();


