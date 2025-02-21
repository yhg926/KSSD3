#include "command_ani.h"
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
#include "../klib/khash.h"

const char gid_obj_prefix[] = "gidobj";
const char ctx_idx_prefix[] = "ctx.index";

inline obj_t extract_outobj_frm_kmer64 (uint64_t kmer, int klen, int holen){ //holen: one side outter obj len
	 	return ((kmer >> (2*( klen - holen )) ) << (2*holen)) | ( (kmer <<(64-2*holen) ) >> (64-2*holen));	
}

int compute_ani(ani_opt_t *ani_opt){ //ref is the sketch(es) to be hashed.
	gen_inverted_index4comblco(ani_opt->refdir);	

}






KHASH_MAP_INIT_INT64(vec_ptr, Vector *)
void gen_inverted_index4comblco(const char *refdir) {
	khash_t(vec_ptr) *h = kh_init(vec_ptr);
	unify_sketch_t* ref_result = generic_sketch_parse(refdir);
	int klen = ref_result->stats.lco_stat_val.klen;
	int holen = ref_result->stats.lco_stat_val.holen;

	for(int rn = 0; rn < ref_result->infile_num; rn++) {
		for (uint64_t ri = ref_result->sketch_index[rn]; ri< ref_result->sketch_index[rn+1]; ri++){
			int ret; uint64_t key = ref_result->comb_sketch[ri];
			khiter_t k = kh_put(vec_ptr, h, key, &ret);
			if (ret == 1){
				Vector vec; vector_init(&vec, sizeof(id_obj_t)); 
				kh_value(h, k) = &vec;
			}
			else if(ret != 0) {				
				kh_destroy(vec_ptr, h);
				err(EXIT_FAILURE,"%s(): Failed to insert key into hashtable",__func__);				
			}				
			id_obj_t id_obj = { rn, extract_outobj_frm_kmer64(key, klen,holen) };
			vector_push(kh_value(h, k), &id_obj); 			
		}
	}  
	//free
	for (khint_t k = kh_begin(h); k != kh_end(h); ++k){
		if (kh_exist(h, k)) {
			printf("%lx",kh_key(h, k));
			for(int i = 0; i< kh_value(h, k)->size ;i++ )  {
			  id_obj_t *id_obj = vector_get(kh_value(h, k),i) ; 
				printf("\t%d|%u", id_obj->gid,id_obj->obj);
			}
			printf("\n");
			vector_free(kh_value(h, k));
		}
	}
	kh_destroy(vec_ptr, h);
}

/*
KHASH_SET_INIT_INT64(kmer_set)
int ani_matrix(ani_opt_t *ani_opt){ //ref is the sketch(es) to be hashed.

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


#define CONFLICT_OBJ (0)
KHASH_MAP_INIT_INT64(u64, uint64_t)
int compute_ani_matrix(matrix_opt_t *matrix_opt){ //ref is the sketch(es) to be hashed.

  unify_sketch_t* ref_result = generic_sketch_parse(matrix_opt->refdir), *qry_result = generic_sketch_parse(matrix_opt->qrydir);
  if(ref_result->stat_type != qry_result->stat_type) err(EXIT_FAILURE,"%s(): ref sketch type %u != qry %u",__func__,ref_result->stat_type,qry_result->stat_type);
  else if(ref_result->hash_id != qry_result->hash_id) err(EXIT_FAILURE,"%s(): ref hash_id %u != qry %u",__func__,ref_result->hash_id,qry_result->hash_id);
	dim_sketch_stat_t *lco_stat_readin = (dim_sketch_stat_t *) ref_result->mem_stat;
	int obj_len = lco_stat_readin->klen - 2* lco_stat_readin->hclen;
	if(obj_len == 0 ) err(EXIT_FAILURE,"%s() abort!: sketching mode has 0bp object",__func__);

	uint64_t tmp_var =  UINT64_MAX >> (64 - 2*lco_stat_readin->hclen) ;
	uint64_t ctxmask = (tmp_var << (2*(lco_stat_readin->klen - lco_stat_readin->hclen - lco_stat_readin->holen)) ) | (tmp_var << (2*(lco_stat_readin->holen) ));
//  int kmerlen = ref_result->kmerlen;
// print header
  FILE *output = matrix_opt->outf[0]=='\0' ? stdout: fopen( matrix_opt->outf, "w");
    if (output == NULL) err(errno,"%s(): %s",__func__, matrix_opt->outf);
  for( int qn = 0;  qn < qry_result->infile_num; qn++)  fprintf(output,"\t%s",qry_result->gname[qn]);
  fprintf(output,"\n");

  co_distance_t *ctx_diff_obj_cnt = malloc(qry_result->infile_num * sizeof(co_distance_t));
  for(int rn = 0; rn < ref_result->infile_num; rn++) {
   //     int X_size = ref_result->sketch_index[rn+1] - ref_result->sketch_index[rn];
    khash_t(u64) *h = kh_init(u64);   int ret;
    // hash rn-th ref genome
    for (uint64_t ri = ref_result->sketch_index[rn]; ri< ref_result->sketch_index[rn+1]; ri++) {
       khiter_t k = kh_put(u64, h, ref_result->comb_sketch[ri] & ctxmask , &ret);
			 if( ret == 1)  kh_value(h, k) = ref_result->comb_sketch[ri];
			 else if (ret == 0) kh_value(h, k) = CONFLICT_OBJ ; //ret == 0: mark confict obj using 0; assuming no duplicted k-mers in a sketch   
		}
    memset(ctx_diff_obj_cnt,0,qry_result->infile_num * sizeof(co_distance_t));
#pragma omp parallel for num_threads ((matrix_opt->p))
      for( int qn = 0;  qn < qry_result->infile_num; qn++) {
        for(uint64_t qi = qry_result->sketch_index[qn]; qi< qry_result->sketch_index[qn+1]; qi++){
					khiter_t it = kh_get(u64, h, qry_result->comb_sketch[qi] & ctxmask );
        	if( it != kh_end(h) && kh_value(h, it) != CONFLICT_OBJ ) {
						 ctx_diff_obj_cnt[qn].ctx_ct++;
						if (kh_value(h, it) != qry_result->comb_sketch[qi]) 	ctx_diff_obj_cnt[qn].diff_obj++;
					}
         } // for qi
    } // for qn
    fprintf(output,"%s",ref_result->gname[rn]);
    for( int qn = 0;  qn < qry_result->infile_num; qn++) {
      //int Y_size =  qry_result->sketch_index[qn+1] - qry_result->sketch_index[qn];
			double ani = 0;
			if(ctx_diff_obj_cnt[qn].ctx_ct != 0) {										
				double dist = (double) ctx_diff_obj_cnt[qn].diff_obj / ctx_diff_obj_cnt[qn].ctx_ct ;
				ani = pow((1 - dist),(1.0/obj_len));
			}
		//	double dist = ctx_diff_obj_cnt[qn].ctx_ct == 0 ? matrix_opt->e : (double) ctx_diff_obj_cnt[qn].diff_obj / ctx_diff_obj_cnt[qn].ctx_ct ;		
      fprintf(output,"\t%lf", ani);
    }
    fprintf(output,"\n");
    kh_destroy(u64,h);
    }// loop rn end
  free_unify_sketch(ref_result); free_unify_sketch(qry_result);
  fclose(output);
  return 1;
}
*/
