#include "command_matrix.h"
#include "global_basic.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <err.h>
#include <errno.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>
#include "klib/khash.h"

KHASH_SET_INIT_INT64(kmer_set)
// core functions
double get_mashD (uint32_t K, uint32_t X, uint32_t Y, uint32_t XnY){
  return (-log( 2*JCD(X,Y,XnY) / (1 + JCD(X,Y,XnY) )) / (K) );
}
double get_aafD (uint32_t K, uint32_t X, uint32_t Y, uint32_t XnY){
  return  (-log(CTM(X,Y,XnY)) / (K));
}
typedef double (*Dist) (uint32_t,uint32_t,uint32_t,uint32_t);
static size_t file_size; 

int compute_triangle(matrix_opt_t *matrix_opt){

	unify_sketch_t* result = generic_sketch_parse(matrix_opt->qrydir) ;		
	uint32_t* tmp_ct_list = malloc(sizeof(uint32_t)*result->infile_num);
    for(int i = 0 ; i < result->infile_num ; i++) tmp_ct_list[i] = result->sketch_index[i+1] - result->sketch_index[i];

    FILE *output = matrix_opt->outf[0]=='\0' ? stdout: fopen( matrix_opt->outf, "w"); 
    if (output == NULL) err(errno,"%s(): %s",__func__, matrix_opt->outf);      
//set distance function
	Dist get_distance;
	if (matrix_opt->metric == 0) get_distance = get_mashD;
	else if (matrix_opt->metric == 1) get_distance = get_aafD;
	else err(errno,"%s(): matrix_opt->metric should be 0 or 1",__func__);
	//read qry sketch and index	
	khash_t(kmer_set) **h_arr = (khash_t(kmer_set) **)malloc(result->infile_num * sizeof(khash_t(kmer_set) *));
	uint32_t *enrolled_qry = calloc(result->infile_num, sizeof(uint32_t)); uint32_t enrolled_num = 0; 
	uint32_t *overlap_kmer_cnt = malloc(result->infile_num * sizeof(uint32_t));
	double *tmp_distance = malloc(result->infile_num * sizeof(double));
    //compare to enrolled_num qry
	for(int qn = 0; qn < result->infile_num; qn++) {							
		uint32_t Y_size = tmp_ct_list[qn];
		int stop_flag = 0; 
#pragma omp parallel for shared(overlap_kmer_cnt, tmp_distance) num_threads(matrix_opt->p)
		for( int er_qn = 0;  er_qn < enrolled_num; er_qn++) { // enrolled qry num
			if (stop_flag) continue; // Skip if stop condition is met
			khash_t(kmer_set) *tmp_ht = h_arr[er_qn];
			overlap_kmer_cnt[er_qn] = 0;

			for (uint32_t qi = result->sketch_index[qn]; qi< result->sketch_index[qn+1]; qi++){
				if( kh_get(kmer_set, tmp_ht, result->comb_sketch[qi]) != kh_end(tmp_ht)) overlap_kmer_cnt[er_qn]++;
			}			
			uint32_t X_size = tmp_ct_list[enrolled_qry[er_qn]];
			if(overlap_kmer_cnt[er_qn] == 0){
				if( matrix_opt->e == -1 ) err(errno, "XnY == 0 abort: %s\t%s\tK=%u\tX=%u\tY=%u\tXnY=%u\n",result->gname[enrolled_qry[er_qn]],result->gname[qn], result->kmerlen,X_size,Y_size,overlap_kmer_cnt[er_qn]);
				else tmp_distance[er_qn] = matrix_opt->e;							
			}
			else tmp_distance[er_qn] = get_distance(result->kmerlen,X_size,Y_size,overlap_kmer_cnt[er_qn]);
	    // Check if the calculated distance meets the condition
    		if (tmp_distance[er_qn] < matrix_opt->c) {
#pragma omp atomic write
        		stop_flag = 1; // Signal other threads to stop processing
    		}
		}//er_qn loop
		//skip if duplicated detected: dist to any enrolled sample < matrix_opt->c //if(dist < matrix_opt->c) continue;
		if (stop_flag) continue;
		fprintf(output, "%s",result->gname[qn]); 			
		for(uint32_t er_qn = 0;  er_qn < enrolled_num; er_qn++) fprintf(output,"\t%lf",tmp_distance[er_qn]);
		if(matrix_opt->d) fprintf(output,"\t%lf",0.0);	fprintf(output,"\n");	
		//hash	
		khash_t(kmer_set) *h = kh_init(kmer_set);	int ret;	
		for (uint32_t qi = result->sketch_index[qn]; qi< result->sketch_index[qn+1]; qi++) kh_put(kmer_set, h, result->comb_sketch[qi], &ret);			
		h_arr[enrolled_num] = h ; 				
		enrolled_qry[enrolled_num++] = qn ;
	}// loop qn end 

	if(matrix_opt->gl[0] != '\0') { //output genome selection code
		FILE *glout;
		if( (glout = fopen(matrix_opt->gl,"w")) == NULL ) err(errno,"%s(): Failed to open file:%s\n",__func__,matrix_opt->gl);		
		int *glist = calloc( result->infile_num, sizeof(int));
		for( int er_qn = 0;  er_qn < enrolled_num; er_qn++) glist[enrolled_qry[er_qn]] = er_qn+1 ; //make sure er_qn start from 0
			
		for(int i = 0; i < result->infile_num ; i++ ){
			if(glist[i] == 0)	fprintf(glout, "%d\tNULL\t%s\n", glist[i],result->gname[i]);
			else fprintf(glout, "%d\t%s\n", glist[i],result->gname[i]);				
		}
		fclose(glout);free(glist);
	}
	//clean
	for( int er_qn = 0;  er_qn < enrolled_num; er_qn++) kh_destroy(kmer_set, h_arr[er_qn]);
	if (output != stdout) fclose(output);
	free_unify_sketch(result);
	return 1;
}


int compute_matrix(matrix_opt_t *matrix_opt){ //ref is the sketch(es) to be hashed.

	unify_sketch_t* ref_result = generic_sketch_parse(matrix_opt->refdir), *qry_result = generic_sketch_parse(matrix_opt->qrydir);
	if(ref_result->stat_type != qry_result->stat_type) err(EXIT_FAILURE,"%s(): ref sketch type %u != qry %u",__func__,ref_result->stat_type,qry_result->stat_type);
	else if(ref_result->hash_id != qry_result->hash_id) err(EXIT_FAILURE,"%s(): ref hash_id %u != qry %u",__func__,ref_result->hash_id,qry_result->hash_id);	
	int kmerlen = ref_result->kmerlen;
//set distance function
	Dist get_distance;
  	if (matrix_opt->metric == 0) get_distance = get_mashD;
  	else if (matrix_opt->metric == 1) get_distance = get_aafD;
  	else err(errno,"compute_matrix(): matrix_opt->metric should be 0 or 1");
// print header
	FILE *output = matrix_opt->outf[0]=='\0' ? stdout: fopen( matrix_opt->outf, "w");
    if (output == NULL) err(errno,"%s(): %s",__func__, matrix_opt->outf);
	for( int qn = 0;  qn < qry_result->infile_num; qn++)	fprintf(output,"\t%s",qry_result->gname[qn]);
	fprintf(output,"\n");

	uint32_t *overlap_kmer_cnt = malloc(qry_result->infile_num * sizeof(uint32_t));
	for(int rn = 0; rn < ref_result->infile_num; rn++) {
        int X_size = ref_result->sketch_index[rn+1] - ref_result->sketch_index[rn];
		khash_t(kmer_set) *h = kh_init(kmer_set);   int ret;
		// hash rn-th ref genome 
		for (uint64_t ri = ref_result->sketch_index[rn]; ri< ref_result->sketch_index[rn+1]; ri++) 
    		kh_put(kmer_set, h, ref_result->comb_sketch[ri], &ret);
		memset(overlap_kmer_cnt,0,qry_result->infile_num * sizeof(uint32_t));
#pragma omp parallel for num_threads(matrix_opt->p)
    	for( int qn = 0;  qn < qry_result->infile_num; qn++) {         
    		for(uint64_t qi = qry_result->sketch_index[qn]; qi< qry_result->sketch_index[qn+1]; qi++){
				if( kh_get(kmer_set, h, qry_result->comb_sketch[qi]) != kh_end(h)) overlap_kmer_cnt[qn]++;
      		} // for qi
		} // for qn
		fprintf(output,"%s",ref_result->gname[rn]);
		for( int qn = 0;  qn < qry_result->infile_num; qn++) {
			int Y_size =  qry_result->sketch_index[qn+1] - qry_result->sketch_index[qn];
			double dist = overlap_kmer_cnt[qn] == 0 ? matrix_opt->e
						: fabs(get_distance(kmerlen,X_size,Y_size,overlap_kmer_cnt[qn])) ;
			fprintf(output,"\t%lf",dist);			
		}		
		fprintf(output,"\n");
		kh_destroy(kmer_set,h);	
  	}// loop rn end
	free_unify_sketch(ref_result); free_unify_sketch(qry_result);
 	fclose(output);
  return 1;
}

