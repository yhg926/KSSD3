#include "command_set_wrapper.h"
#include "command_set.h"
#include "command_operate.h"
/*** argp wraper ***/
struct arg_set
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_set[] =
{
  {"union",'u', 0,  0, "get union set of the sketches.\v",1 },
	{"subtract",'s',"<pan>", 0,"subtract the pan-sketch from all input sketches.\v",2 },
	{"intsect",'i',"<pan>", 0, "intersect with the pan-sketch for all input sketches.\v",2},
	{"uniq_union",'q',0,  0, "get uniq union set of the sketches.\v",3 },
    {"markerdb",333,0,  0, "generate markerdb instead of uniq union set, -q must be set.\v",4 },
//	{"combin_pan",'c',0,  0, "combine pan files to combco file.\v",4 },
	{"threads",'p',"<INT>",  0, "number of threads.\v",5 },
	{"print",'P',0,  0, "print genome names.\v",5 },
	{"psketch",777,0,  0, "print sketch content.\v",5 },
	{"grouping",'g',"<file.tsv>",0,"grouping genomes by input category file.\v",5},
	{"outdir",'o',"<path>",0,"specify the output directory.\v",6},
  { 0 }
};

static char doc_set[] =
  "\n"
  "The set doc prefix."
  "\v"
  "The set doc suffix."
  ;


set_opt_t set_opt = {
.operation = -1,//0:subtract,1:intersect,2 union, 3 uniq_union, 4 combin_pan
.q2markerdb = 0, // when -q set, generate markerdb instead of uniq union set, only for lco sketch 
.p = 1,
.P = 0,
.show = 0,
.num_remaining_args = 0,
.remaining_args = NULL,			
.insketchpath[0] = '\0',
.pansketchpath[0]='\0',
.subsetf[0] = '\0',
.outdir = "./"
};

static error_t parse_set(int key, char* arg, struct argp_state* state) {

  struct arg_set* set = state->input;
  assert( set );
  assert( set->global );
	
  switch(key)
  {

    case 'u':
		{ 
			if (set_opt.operation != -1 ) printf("set operation is already set, -u is ignored.\n");
			else set_opt.operation = 2 ;
			break;
		}
		case 'q':
    {
			
      if (set_opt.operation != -1 ) printf("set operation is already set, -q is ignored.\n");
      else set_opt.operation = 3 ;
      break;
    }
		case 's':
		{
			if (set_opt.operation != -1) printf("set operation is already set, -s is ignored.\n");
			else {
				set_opt.operation = 0;
				strcpy(set_opt.pansketchpath, arg);
			}
			break;
		}
		case 'i':
		{
		
			if (set_opt.operation != -1 ) printf("set operation is already set, -i is ignored.\n");
			else {
				set_opt.operation = 1 ;            
				strcpy(set_opt.pansketchpath, arg);
			}
			break;
		}
		case 'c':
		{
			if (set_opt.operation != -1 ) printf("set operation is already set, -c is ignored.\n");
			else set_opt.operation = 4 ;
			break;
		}
		case 'o':
		{
		
			strcpy(set_opt.outdir, arg);
			 
			break;
		}
		case 'p':
		{
			set_opt.p = atoi(arg); 
			break;
		}
		case 'P':
		{
			set_opt.P = 1;	
			break;
		}
		case 'g':
		{
			strcpy(set_opt.subsetf, arg);
			break;
		}
		case 333:
		{
			set_opt.q2markerdb = 1; 
			break;
		}
		case 777:
		{
			set_opt.show = 1;
			break;	
		}
		case ARGP_KEY_ARGS:
			strcpy(set_opt.insketchpath, state->argv[state->next]);
			set_opt.num_remaining_args = state->argc - state->next;
			set_opt.remaining_args  = state->argv + state->next;
			break;
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

static struct argp argp_set =
{
  opt_set,
  parse_set,
	"<combined sketch>", //0  "[arguments ...]",
  doc_set
};

int cmd_set(struct argp_state* state)
{
  struct arg_set set = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char*  argv0 =  argv[0];

  set.global = state->input;
	argv[0] = malloc(strlen(state->name) + strlen(" set") + 1);

  if(!argv[0])
    argp_failure(state, 1, ENOMEM, 0);
	sprintf(argv[0], "%s set", state->name);
  argp_parse(&argp_set, argc, argv, ARGP_IN_ORDER, &argc, &set);

	free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;
	// operation and arg control
	if(argc >1){	
		if(set_opt.operation == 2){
			if(file_exists_in_folder(set_opt.insketchpath,co_dstat) )
				return sketch_union(&set_opt); 
			else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat))
				return lsketch_union(&set_opt);			
		}
		else if(set_opt.operation == 3){
	  	if(file_exists_in_folder(set_opt.insketchpath,co_dstat))
				return uniq_sketch_union(&set_opt) ;
      else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat))
        return lsketch_union(&set_opt);
		}
		else if(set_opt.operation == 4){
			return combin_pans(&set_opt);
		}
		else if(set_opt.operation == 0 || set_opt.operation == 1 ){
			if(file_exists_in_folder(set_opt.pansketchpath,co_dstat))
				return sketch_operate(&set_opt) ;
			else if (file_exists_in_folder(set_opt.pansketchpath,sketch_stat))
				 return lsketch_operate(&set_opt) ;
		}
		else {
			if(set_opt.P) {
				if(file_exists_in_folder(set_opt.insketchpath,co_dstat)) print_gnames(&set_opt);
				else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat)) print_lco_gnames(&set_opt);
				else printf("%s is not a valid sketch\n",set_opt.insketchpath );

			}
			else if(set_opt.show > 0){
					if(file_exists_in_folder(set_opt.insketchpath,sketch_stat)) show_content(&set_opt);
			}
			else if (set_opt.subsetf[0]!='\0') {
				if(file_exists_in_folder(set_opt.insketchpath,co_dstat))
					return grouping_genomes(&set_opt); // combin_subset_pans(set_opt.subsetf);
				else if(file_exists_in_folder(set_opt.insketchpath,sketch_stat))
					return lgrouping_genomes(&set_opt);

			}
			else printf("set operation use : -u, -q, -i or -s\n");
			return -1 ;
		}
	}
	else
		return -1;
}
