#include "command_matrix.h"
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

// core functions

double get_mashD (int K, int X, int Y, int XnY){
//	printf("K=%d,X=%d,Y=%d,XnY=%d,mashd=%lf\n",K,X,Y,XnY, -log( 2*JCD(X,Y,XnY) / (1 + JCD(X,Y,XnY) )) / (K) );
  return (-log( 2*JCD(X,Y,XnY) / (1 + JCD(X,Y,XnY) )) / (K) );
}

double get_aafD (int K, int X, int Y, int XnY){
  return  (-log(CTM(X,Y,XnY)) / (K));
}

typedef double (*Dist) (int,int,int,int);

int compute_triangle(matrix_opt_t *matrix_opt){
	const char *qry_dstat_fpath = test_get_fullpath(matrix_opt->qrydir,co_dstat);
	if (qry_dstat_fpath == NULL) err(errno,"cannot find %s under %s ",co_dstat, qry_dstat_fpath);
	FILE *qry_dstat_fp;
  if( (qry_dstat_fp = fopen(qry_dstat_fpath,"rb")) == NULL ) err(errno, "compute_triangle():%s",qry_dstat_fpath);
	co_dstat_t qry_dstat ;
	fread( &qry_dstat, sizeof(co_dstat_t),1,qry_dstat_fp);

	if(qry_dstat.comp_num > 1)
		err(errno,"this version compute_triangle() was for comp_num == 1, query sketch has comp_num:%d",qry_dstat.comp_num);

	FILE *output; 

  if (matrix_opt->outf[0]!='\0') {
        output = fopen( matrix_opt->outf, "w");
        if (output == NULL) {
            printf("Failed to open file: %s", matrix_opt->outf);
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
	if (matrix_opt->metric == 0) 
		get_distance = get_mashD;
	else if (matrix_opt->metric == 1)
	  get_distance = get_aafD;
	else
		err(errno,"compute_triangle(): matrix_opt->metric should be 0 or 1");

	for(int c = 0;  c < 1; c++){ //c < qry_dstat.comp_inum for each kmer subspace component c
		//read qry sketch
		sprintf(tmpfname,"%s/%s.%d",matrix_opt->qrydir,skch_prefix,c);
  	if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"compute_triangle():%s",tmpfname);
    stat(tmpfname, &tmpstat);
    unsigned int *qry_combco = malloc(tmpstat.st_size);
    fread(qry_combco, tmpstat.st_size, 1, tmpfp);
    fclose(tmpfp);
    //read qry index
    sprintf(tmpfname,"%s/%s.%d",matrix_opt->qrydir,idx_prefix,c);
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
#pragma omp parallel for shared(overlap_kmer_cnt, tmp_distance) num_threads(matrix_opt->p)
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
				if(overlap_kmer_cnt[er_qn] == 0){
					if( matrix_opt->e == -1 ) {
						fprintf(stderr, "XnY == 0 abort: %s\t%s\tK=%d\tX=%d\tY=%d\tXnY=%d\n",qryname[enrolled_qry[er_qn]],qryname[qn], kmerlen,X_size,Y_size,overlap_kmer_cnt[er_qn]);
						exit(1);
					}
					else 
						tmp_distance[er_qn] = matrix_opt->e;							
				}
				else
					tmp_distance[er_qn] = get_distance(kmerlen,X_size,Y_size,overlap_kmer_cnt[er_qn]);
	
	
	    // Check if the calculated distance meets the condition
    		if (tmp_distance[er_qn] < matrix_opt->c) {
#pragma omp atomic write
        	stop_flag = 1; // Signal other threads to stop processing
    		}
			}
			//skip if duplicated detected: dist to any enrolled sample < matrix_opt->c 
			//if(dist < matrix_opt->c) continue;
			if (stop_flag) continue;
			fprintf(output, "%s",qryname[qn]); 			
			for(int er_qn = 0;  er_qn < enrolled_num; er_qn++) fprintf(output,"\t%lf",tmp_distance[er_qn]);
			if(matrix_opt->d) fprintf(output,"\t%lf",0.0);
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

		if(matrix_opt->gl[0] != '\0') { //output genome selection code
			FILE *glout;
			if( (glout = fopen(matrix_opt->gl,"w")) == NULL ){
				printf("Failed to open file:%s\n",matrix_opt->gl);
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


int compute_matrix(matrix_opt_t *matrix_opt){ //ref is the sketch(es) to be hashed.
	const char *ref_dstat_fpath = test_get_fullpath(matrix_opt->refdir,co_dstat);
	if (ref_dstat_fpath == NULL) err(errno,"cannot find %s under %s ",co_dstat, ref_dstat_fpath);
  const char *qry_dstat_fpath = test_get_fullpath(matrix_opt->qrydir,co_dstat);
  if (qry_dstat_fpath == NULL) err(errno,"cannot find %s under %s ",co_dstat, qry_dstat_fpath);

  FILE *ref_dstat_fp, *qry_dstat_fp;

  if( (ref_dstat_fp = fopen(ref_dstat_fpath,"rb")) == NULL ) err(errno, "compute_matrix():%s",ref_dstat_fpath);
  co_dstat_t ref_dstat,qry_dstat ;
  fread( &ref_dstat, sizeof(co_dstat_t),1,ref_dstat_fp);	
	
  if( (qry_dstat_fp = fopen(qry_dstat_fpath,"rb")) == NULL ) err(errno, "compute_matrix():%s",qry_dstat_fpath);
  fread( &qry_dstat, sizeof(co_dstat_t),1,qry_dstat_fp);

	if(ref_dstat.shuf_id != qry_dstat.shuf_id)
		err(errno,"ref shuf_id %d != qry shuf_id %d",ref_dstat.shuf_id,qry_dstat.shuf_id);	
  if(qry_dstat.comp_num > 1 || ref_dstat.comp_num > 1)
    err(errno,"this version compute_matrix() was for comp_num == 1, ref/query sketch has comp_num:%d/%d",ref_dstat.comp_num,qry_dstat.comp_num);

  FILE *output;

  if (matrix_opt->outf[0]!='\0') {
        output = fopen( matrix_opt->outf, "w");
        if (output == NULL) {
            printf("Failed to open file: %s", matrix_opt->outf);
            return 1;
        }
  } else {
        output = stdout;  // Default to stdout
  }
  int kmerlen = qry_dstat.kmerlen;

  ctx_obj_ct_t* ref_ct_list = malloc(sizeof(ctx_obj_ct_t) * ref_dstat.infile_num);
  fread(ref_ct_list,sizeof(ctx_obj_ct_t),ref_dstat.infile_num,ref_dstat_fp);

  ctx_obj_ct_t* qry_ct_list = malloc(sizeof(ctx_obj_ct_t) * qry_dstat.infile_num);
  fread(qry_ct_list,sizeof(ctx_obj_ct_t),qry_dstat.infile_num,qry_dstat_fp);

	char (*refname)[PATHLEN] = malloc(PATHLEN * ref_dstat.infile_num);
  fread(refname,PATHLEN,ref_dstat.infile_num, ref_dstat_fp);
  fclose(ref_dstat_fp);
  char (*qryname)[PATHLEN] = malloc(PATHLEN * qry_dstat.infile_num);
  fread(qryname,PATHLEN,qry_dstat.infile_num, qry_dstat_fp);
  fclose(qry_dstat_fp);

  char tmpfname[PATHLEN];
  struct stat tmpstat;
  FILE *tmpfp;
//set distance function
  Dist get_distance;
  if (matrix_opt->metric == 0)
    get_distance = get_mashD;
  else if (matrix_opt->metric == 1)
    get_distance = get_aafD;
  else
    err(errno,"compute_matrix(): matrix_opt->metric should be 0 or 1");
		//read ref genome
    sprintf(tmpfname,"%s/%s.0",matrix_opt->refdir,skch_prefix);
    if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"compute_matrix():%s",tmpfname);
    stat(tmpfname, &tmpstat);
    unsigned int *ref_combco = malloc(tmpstat.st_size);
    fread(ref_combco, tmpstat.st_size, 1, tmpfp);
    fclose(tmpfp);
    //read ref index
    sprintf(tmpfname,"%s/%s.0",matrix_opt->refdir,idx_prefix);
    if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"compute_matrix():%s",tmpfname);
    stat(tmpfname, &tmpstat);
    size_t *ref_index_combco = malloc(tmpstat.st_size);
    fread(ref_index_combco,tmpstat.st_size, 1, tmpfp);
    fclose(tmpfp);

    //read qry genome
    sprintf(tmpfname,"%s/%s.0",matrix_opt->qrydir,skch_prefix);
    if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"compute_matrix():%s",tmpfname);
    stat(tmpfname, &tmpstat);
    unsigned int *qry_combco = malloc(tmpstat.st_size);
    fread(qry_combco, tmpstat.st_size, 1, tmpfp);
    fclose(tmpfp);
    //read qry index
    sprintf(tmpfname,"%s/%s.0",matrix_opt->qrydir,idx_prefix);
    if( (tmpfp = fopen(tmpfname,"rb"))==NULL) err(errno,"compute_matrix():%s",tmpfname);
    stat(tmpfname, &tmpstat);
    size_t *qry_index_combco = malloc(tmpstat.st_size);
    fread(qry_index_combco,tmpstat.st_size, 1, tmpfp);
    fclose(tmpfp);
	// print header
	for( int qn = 0;  qn < qry_dstat.infile_num; qn++)	
		fprintf(output,"\t%s",qryname[qn]);
	fprintf(output,"\n");

	int *overlap_kmer_cnt = malloc(qry_dstat.infile_num * sizeof(int));
	for(int rn = 0; rn < ref_dstat.infile_num; rn++) {
  	int X_size = ref_ct_list[rn];
		int hash_sz = nextPrime( (int)((double)(ref_index_combco[rn+1] - ref_index_combco[rn]) / LD_FCTR) );
		unsigned int *tmp_ht = (unsigned int *)calloc( hash_sz, sizeof(unsigned int) ) ;
		// hash rn-th ref genome 
		for (size_t ri = ref_index_combco[rn]; ri< ref_index_combco[rn+1]; ri++){
          //hash table insert
          for(int i = 0 ; i< hash_sz; i++) {
            unsigned int key = HASH(ref_combco[ri],i, hash_sz);
            if (tmp_ht[key] == 0) {
              tmp_ht[key] = ref_combco[ri];
              break;
            }
          }
      }
    
		memset(overlap_kmer_cnt,0,qry_dstat.infile_num * sizeof(int));
#pragma omp parallel for num_threads(matrix_opt->p)
    for( int qn = 0;  qn < qry_dstat.infile_num; qn++) {         
     	
    	for(size_t qi = qry_index_combco[qn]; qi< qry_index_combco[qn+1]; qi++){
          //hash table lookup
        for(int i = 0 ; i< hash_sz; i++) {
        	unsigned int key = HASH(qry_combco[qi],i, hash_sz);
          if (tmp_ht[key] == 0) break;
          else if (tmp_ht[key] == qry_combco[qi]){
            overlap_kmer_cnt[qn]++;
            break;
          }
        }
      } // for qi
	  } // for qn
		free(tmp_ht);
		
		fprintf(output,"%s",refname[rn]);
		for( int qn = 0;  qn < qry_dstat.infile_num; qn++) {
			double dist = overlap_kmer_cnt[qn] == 0 ? matrix_opt->e
			: get_distance(kmerlen,X_size,qry_ct_list[qn],overlap_kmer_cnt[qn]) ;
			fprintf(output,"\t%lf",dist);			
		}		
		fprintf(output,"\n");

  }// loop rn end

	free(refname);
  free(qryname);
  free(ref_ct_list);
	free(qry_ct_list);
	free(ref_combco);
	free(ref_index_combco);
  free(qry_combco);
  free(qry_index_combco);
  fclose(output);
  return 1;
}

