#include "command_triangle.h"
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
#include <omp.h>
#include <stdatomic.h>

/*** argp wrapper ***/
struct arg_triangle
{
  struct arg_global* global;

  char* name;
};

static struct argp_option opt_triangle[] =
{
	{"query",'q',"<DIR>", 0, "Path of query sketches \v",1},
	{"metric",'m',"<0/1>", 0, "Using mashD or aafD (0/1) [0]\v",2},
	{"control",'c',"<FLOAT>",0,"Skip duplicated samples (distance < c) [0] \v",3},
	{"glist",'g',"<FILE>",0,"Output path\v",4},
	{"outfile",'o',"<FILE>",0,"Output path [STDOUT]\v",5},
	{"threads",'p',"<INT>", 0, "Threads number to use \v",6},
	{"diagonal",'d',0, 0, "set diagonal\v",7},
  { 0 }
};

static char doc_triangle[] =
  ""
  ""
  "\v"
  ""
  ;

triangle_opt_t triangle_opt ={
	.metric = 0, // 0:mashD, 1:aafD
	.c = 0.0, //control duplicated sample by skip distance < c;
	.p = 1,
	.d = 0, //diagonal
	.qrydir[0] = '\0',
	.outf[0] = '\0',
	.gl[0]= '\0',
	.num_remaining_args = 0, //int num_remaining_args; no option arguments num.
	.remaining_args = NULL //char **remaining_args; no option arguments array.
};

static error_t parse_triangle(int key, char* arg, struct argp_state* state) {
  struct arg_triangle* triangle = state->input;
	
  switch(key)
  {		
		case 'm':
		{
				triangle_opt.metric = atoi(arg);
				break;
		}
		case 'c':
		{
				triangle_opt.c = atof(arg);
				break;
		}
		case 'p':
		{
			triangle_opt.p = atoi(arg);
			break;
		}
		case 'q':
		{
			strcpy(triangle_opt.qrydir, arg);
			break;
		}
		case 'o':
		{
			strcpy(triangle_opt.outf, arg);
			break;
		}
		case 'g':
		{
			strcpy(triangle_opt.gl, arg);
			break;
		}
		case 'd':
		{
			triangle_opt.d = 1 ;
			break;
		}
		case ARGP_KEY_ARGS:
		{
			triangle_opt.num_remaining_args = state->argc - state->next;
			triangle_opt.remaining_args  = state->argv + state->next;
			break;
		}
    case ARGP_KEY_END:
    {	
			if(triangle_opt.qrydir[0] == '\0'){
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

static struct argp argp_triangle =
{
  opt_triangle,
  parse_triangle,
	0,//  "[arguments ...]",
  doc_triangle
};

int cmd_triangle(struct argp_state* state)
{
  struct arg_triangle triangle = { 0, };
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  triangle.global = state->input;
  argp_parse(&argp_triangle, argc, argv, ARGP_IN_ORDER, &argc, &triangle);

  state->next += argc - 1;
  if (triangle_opt.qrydir[0] != '\0')
    return compute_triangle(&triangle_opt);
  return 1;
}

// core functions

double get_mashD (int K, int X, int Y, int XnY){
//	printf("K=%d,X=%d,Y=%d,XnY=%d,mashd=%lf\n",K,X,Y,XnY, -log( 2*JCD(X,Y,XnY) / (1 + JCD(X,Y,XnY) )) / (K) );
  return (-log( 2*JCD(X,Y,XnY) / (1 + JCD(X,Y,XnY) )) / (K) );
}

double get_aafD (int K, int X, int Y, int XnY){
  return  (-log(CTM(X,Y,XnY)) / (K));
}

typedef double (*Dist) (int,int,int,int);

int compute_triangle(triangle_opt_t *triangle_opt){
	const char *qry_dstat_fpath = test_get_fullpath(triangle_opt->qrydir,co_dstat);
	if (qry_dstat_fpath == NULL) err(errno,"cannot find %s under %s ",co_dstat, qry_dstat_fpath);
	FILE *qry_dstat_fp;
  if( (qry_dstat_fp = fopen(qry_dstat_fpath,"rb")) == NULL ) err(errno, "compute_triangle():%s",qry_dstat_fpath);
	co_dstat_t qry_dstat ;
	fread( &qry_dstat, sizeof(co_dstat_t),1,qry_dstat_fp);

	if(qry_dstat.comp_num > 1)
		err(errno,"this version compute_triangle() was for comp_num == 1, query sketch has comp_num:%d",qry_dstat.comp_num);

	FILE *output; 

  if (triangle_opt->outf[0]!='\0') {
        output = fopen( triangle_opt->outf, "w");
        if (output == NULL) {
            printf("Failed to open file: %s", triangle_opt->outf);
            return 1;
        }
  } else {
        output = stdout;  // Default to stdout
  }


	int kmerlen = qry_dstat.kmerlen;

	ctx_obj_ct_t* tmp_ct_list = malloc(sizeof(ctx_obj_ct_t) * qry_dstat.infile_num);
	fread(tmp_ct_list,sizeof(ctx_obj_ct_t),qry_dstat.infile_num,qry_dstat_fp);

	char (*qryname)[PATHLEN] = malloc(PATHLEN * qry_dstat.infile_num);
	fread(qryname,PATHLEN,qry_dstat.infile_num, qry_dstat_fp);
	fclose(qry_dstat_fp);

  char tmpfname[PATHLEN];
  struct stat tmpstat;
  FILE *tmpfp;
//set distance function
	Dist get_distance;
	if (triangle_opt->metric == 0) 
		get_distance = get_mashD;
	else if (triangle_opt->metric == 1)
	  get_distance = get_aafD;
	else
		err(errno,"compute_triangle(): triangle_opt->metric should be 0 or 1");

	for(int c = 0;  c < 1; c++){ //c < qry_dstat.comp_inum for each kmer subspace component c
		//read qry sketch
		sprintf(tmpfname,"%s/%s.%d",triangle_opt->qrydir,skch_prefix,c);
  	if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"compute_triangle():%s",tmpfname);
    stat(tmpfname, &tmpstat);
    unsigned int *qry_combco = malloc(tmpstat.st_size);
    fread(qry_combco, tmpstat.st_size, 1, tmpfp);
    fclose(tmpfp);
    //read qry index
    sprintf(tmpfname,"%s/%s.%d",triangle_opt->qrydir,idx_prefix,c);
    if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"compute_triangle():%s",tmpfname);
    stat(tmpfname, &tmpstat);
    size_t *qry_index_combco = malloc(tmpstat.st_size);
    fread(qry_index_combco,tmpstat.st_size, 1, tmpfp);
    fclose(tmpfp);

		int *ht_size = malloc(qry_dstat.infile_num * sizeof(int));
		unsigned int **kmer_ht_array = malloc( sizeof(unsigned int *) * qry_dstat.infile_num); 
		int *enrolled_qry = calloc( qry_dstat.infile_num, sizeof(int));
		int enrolled_num = 0; 
		int *overlap_kmer_cnt = malloc(qry_dstat.infile_num * sizeof(int));
		double *tmp_distance = malloc(qry_dstat.infile_num * sizeof(double));

		for(int qn = 0; qn < qry_dstat.infile_num; qn++) {						
			//compare to enrolled_num qry
		//	double dist = 1;
			int Y_size = tmp_ct_list[qn];
			int stop_flag = 0; 
#pragma omp parallel for shared(overlap_kmer_cnt, tmp_distance) num_threads(triangle_opt->p)
			for( int er_qn = 0;  er_qn < enrolled_num; er_qn++) { // enrolled qry num
				if (stop_flag) continue; // Skip if stop condition is met
				unsigned int *tmp_ht = kmer_ht_array[er_qn] ;
				int hash_sz = ht_size[er_qn]; 
				overlap_kmer_cnt[er_qn] = 0;

				for (size_t qi = qry_index_combco[qn]; qi< qry_index_combco[qn+1]; qi++){
					//hash table lookup
					for(int i = 0 ; i< hash_sz; i++) {
						unsigned int key = HASH(qry_combco[qi],i, hash_sz);
						if (tmp_ht[key] == 0) break;
						else if (tmp_ht[key] == qry_combco[qi]){
							overlap_kmer_cnt[er_qn]++;
							break;
						}
					}
				}
				int X_size = tmp_ct_list[enrolled_qry[er_qn]];
//				dist = tmp_distance[er_qn] = get_distance(kmerlen,X_size,Y_size,overlap_kmer_cnt[er_qn]);
				tmp_distance[er_qn] = get_distance(kmerlen,X_size,Y_size,overlap_kmer_cnt[er_qn]);
	    // Check if the calculated distance meets the condition
    		if (tmp_distance[er_qn] < triangle_opt->c) {
#pragma omp atomic write
        	stop_flag = 1; // Signal other threads to stop processing
    		}
			}
			//skip if duplicated detected: dist to any enrolled sample < triangle_opt->c 
			//if(dist < triangle_opt->c) continue;
			if (stop_flag) continue;
			fprintf(output, "%s",qryname[qn]); 			
			for(int er_qn = 0;  er_qn < enrolled_num; er_qn++) fprintf(output,"\t%lf",tmp_distance[er_qn]);
			if(triangle_opt->d) fprintf(output,"\t%lf",0.0);
			fprintf(output,"\n");	

			//hash			
      int hash_sz = nextPrime( (int)((double)(qry_index_combco[qn+1] - qry_index_combco[qn]) / LD_FCTR) );		
      unsigned int *tmp_ht = (unsigned int *)calloc( hash_sz, sizeof(unsigned int) ) ;

			for (size_t qi = qry_index_combco[qn]; qi< qry_index_combco[qn+1]; qi++){
          //hash table insert
          for(int i = 0 ; i< hash_sz; i++) {
            unsigned int key = HASH(qry_combco[qi],i, hash_sz);
            if (tmp_ht[key] == 0) {
							tmp_ht[key] = qry_combco[qi];
							break;
						}
            else if (tmp_ht[key] == qry_combco[qi])
              break;            
          }	
			}
			ht_size[enrolled_num] = hash_sz ;
			kmer_ht_array[enrolled_num] = tmp_ht ; 				
			enrolled_qry[enrolled_num++] = qn ;

		}// loop qn end 

		if(triangle_opt->gl[0] != '\0') { //output genome selection code
			FILE *glout;
			if( (glout = fopen(triangle_opt->gl,"w")) == NULL ){
				printf("Failed to open file:%s\n",triangle_opt->gl);
				return 1;
			}
			int *glist = calloc( qry_dstat.infile_num, sizeof(int));
			for( int er_qn = 0;  er_qn < enrolled_num; er_qn++) 
				glist[enrolled_qry[er_qn]] = er_qn+1 ; //make sure er_qn start from 0
			
			for(int i = 0; i < qry_dstat.infile_num ; i++ ){
				if(glist[i] == 0)	
					fprintf(glout, "%d\tNULL\t%s\n", glist[i],qryname[i]);
				else
					fprintf(glout, "%d\t%s\n", glist[i],qryname[i]);				
			}
			fclose(glout);
			free(glist);

		}

		for( int er_qn = 0;  er_qn < enrolled_num; er_qn++) free(kmer_ht_array[er_qn] );
		free(ht_size);
		free(kmer_ht_array);
		free(enrolled_qry);
		free(overlap_kmer_cnt);
		free(tmp_distance);

	}// loop component c end

	free(qryname);
	free(tmp_ct_list);

	if (output != stdout) 
  	fclose(output);
 

	return 1;

}


