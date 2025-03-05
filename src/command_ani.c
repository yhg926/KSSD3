#include "command_ani.h"
#include "command_matrix.h"
#include "global_basic.h"
#include "kssdlib_sort.h"
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
//#include "../klib/khashl.h"
#define GID_NBITS 20 //2^20, 1M
#define CONFLICT_OBJ UINT32_MAX
//pulic vars
const char gid_obj_prefix[] = "gidobj", ctx_idx_prefix[] = "ctx.index";
uint64_t ho_mask_len, hc_mask_len, io_mask_len, ho_mask_left,hc_mask_left, io_mask, hc_mask_right, ho_mask_right,ctxmask;
int iolen, klen,hclen,holen ; uint32_t gid_mask; bitslen_t Bitslen; // context, gid, and obj bits len

KHASH_MAP_INIT_INT64(u64, uint32_t)
int compute_ani(ani_opt_t *ani_opt){ //ref is the sketch(es) to be hashed.

	gen_inverted_index4comblco(ani_opt->refdir);	

}

void check_comb_sortedsketch64 ( unify_sketch_t* result);
void gen_inverted_index4comblco(const char *refdir) {

	time_t seconds = time(NULL) ;
	unify_sketch_t* ref_result = generic_sketch_parse(refdir);
	const_comask_init(&ref_result->stats.lco_stat_val);
	uint64_t sketch_size = ref_result->sketch_index[ref_result->infile_num] ; 
    if (sketch_size > (float) UINT32_MAX * LD_FCTR ) err(EXIT_FAILURE,"%s():sketch_index maximun %lu exceed UINT32_MAX*LF;%f",__func__, sketch_size,(float) UINT32_MAX * LD_FCTR);
	if( ref_result->infile_num >= ( 1<<GID_NBITS )) err(EXIT_FAILURE,"%s(): genome numer %d exceed maximum:%u",__func__,ref_result->infile_num,1<<GID_NBITS);
	if(GID_NBITS + 4*hclen >64) err(EXIT_FAILURE,"%s(): context_bits_len(%d)+gid_bits_len(%d) exceed 64",__func__,4*hclen,GID_NBITS);
	printf("time1=%ld\n", time(NULL) - seconds );//	check_comb_sortedsketch64 ( ref_result);
	comb_sortedsketch64Xcomb_sortedsketch64 ( ref_result, ref_result );
  printf("time2=%ld\n", time(NULL) - seconds );
	ctxgidobj_t *ctxgidobj_arr = ctxobj64_2ctxgidobj(ref_result->sketch_index, ref_result->comb_sketch, ref_result->infile_num, sketch_size, klen,hclen,holen);
  printf("time3=%ld\n", time(NULL) - seconds );
	ctxgidobj_sort_array(ctxgidobj_arr, sketch_size) ; 
  printf("time4=%ld\n", time(NULL) - seconds );
	sort_sketch_summary_t *sort_sketch_summary = summarize_ctxgidobj_arr(ctxgidobj_arr, ref_result->sketch_index, sketch_size, ref_result->infile_num);
  printf("time5=%ld\n", time(NULL) - seconds );
  sorted_ctxgidobj_arrXcomb_sortedsketch64 ( ref_result, ctxgidobj_arr, sort_sketch_summary );
printf("time6=%ld\n", time(NULL) - seconds );
	free_sort_sketch_summary(sort_sketch_summary);
  free_unify_sketch(ref_result);	
	free(ctxgidobj_arr);
}

void const_comask_init(dim_sketch_stat_t *lco_stat_val ){
//	init all public vars ;
	klen =  lco_stat_val->klen;
	holen = lco_stat_val->holen;
	hclen = lco_stat_val->hclen; 
	iolen = klen - 2*(hclen + holen);
 
  Bitslen.ctx = 4*hclen;
	Bitslen.gid = GID_NBITS;
	Bitslen.obj = 2*klen-4*hclen ;

	ho_mask_len = holen == 0? 0: UINT64_MAX >> (64 - 2*holen);
  hc_mask_len = hclen == 0? 0: UINT64_MAX >> (64 - 2*hclen) ;
  io_mask_len = iolen == 0? 0: UINT64_MAX >> (64 - 2*iolen) ; //UINT64_MAX >> 64 is undefined  not 0

  ho_mask_left = ho_mask_len << (2*(klen-holen));//(2*(holen + hclen + iolen + hclen));
  hc_mask_left = hc_mask_len << (2*(holen + hclen + iolen));//2*(holen + hclen + iolen));
  io_mask = io_mask_len << (2*(holen + hclen));
  hc_mask_right = hc_mask_len << (2*holen);
  ho_mask_right = ho_mask_len;
	
  ctxmask = hc_mask_left | hc_mask_right ;
  gid_mask = (1U<< GID_NBITS) - 1 ;
}



void sketch64_2ctxobj64 (uint64_t *sketch_index, uint64_t *sketch64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen ){
#pragma omp parallel for num_threads(32) schedule(guided)
        for (uint32_t ri = 0; ri< arrlen; ri++) {
            sketch64[ri] = ((sketch64[ri] & hc_mask_left) << 2*holen) |
                    	 ((sketch64[ri] & hc_mask_right) << 2*(holen + iolen)  )         |
                    	  ((sketch64[ri] & ho_mask_left) >> 4*hclen ) | 
						((sketch64[ri] & ho_mask_right) << (2*iolen) )|
						((sketch64[ri] & io_mask) >> (2*(holen + hclen ))) ;

        }
   
}

uint96_t * sketch64_2uint96co(uint64_t *sketch_index, uint64_t *sketch64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen ){
//	int iolen = klen - 2*(hclen + holen);
	uint96_t *uint96co =  malloc( arrlen * sizeof(uint96_t)) ;
#pragma omp parallel for num_threads(32) schedule(guided)
	for( uint32_t rn = 0; rn < infile_num; rn++){
        for (uint32_t ri = sketch_index[rn]; ri< sketch_index[rn+1]; ri++) {
			uint64_t kmer = sketch64[ri];	tmp_ctxgidobj_t tmp_ctxgidobj;	
			tmp_ctxgidobj.obj = ((kmer & ho_mask_left) >> 4*hclen ) | ((kmer & ho_mask_right) << (2*iolen) )|((kmer & io_mask) >> (2*(holen + hclen )))  ; //obj , lowest 
			tmp_ctxgidobj.gid = rn;
			tmp_ctxgidobj.ctx = ((kmer & hc_mask_left) >> 2*(holen  + iolen) ) | ((kmer & hc_mask_right) >> (2*holen));
			uint96co[ri] = concat_lower_bits(&tmp_ctxgidobj, Bitslen);
        }
    }
	return uint96co;
}


ctxgidobj_t* ctxobj64_2ctxgidobj(uint64_t *sketch_index, uint64_t *ctxobj64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen ){
	ctxgidobj_t* ctxgidobj = malloc( arrlen * sizeof(ctxgidobj_t));
	uint64_t obj_len_mask = (1LU<< 2*( 2*holen + iolen)) - 1; 
#pragma omp parallel for num_threads(32) schedule(guided)
	for( uint32_t rn = 0; rn < infile_num; rn++){
		for (uint32_t ri = sketch_index[rn]; ri< sketch_index[rn+1]; ri++) {
			ctxgidobj[ri].ctxgid = (ctxobj64[ri]  >> 2*( 2*holen + iolen) << GID_NBITS) | (rn & gid_mask );
			ctxgidobj[ri].obj = ctxobj64[ri] & obj_len_mask ;
			
		}
	}	
	return ctxgidobj;
}


sort_sketch_summary_t * summarize_ctxgidobj_arr(ctxgidobj_t* ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num){

	uint64_t gidmask = UINT64_MAX >> (64-GID_NBITS) ;
	num_ctx_cfltobj_t *num_ctx_cfltobj_arr = malloc(infile_num*sizeof(num_ctx_cfltobj_t));
	for(int i = 0 ; i < infile_num; i++) {
		num_ctx_cfltobj_arr[i].num_ctx = sketch_index[i+1] - sketch_index[i];	
		num_ctx_cfltobj_arr[i].num_conflictobj = 0;
	}
	
	uint32_t num_ctx = 1;
	for(uint32_t i = 1 ; i < arrlen  ; i++ ){
		if(ctxgidobj_arr[i].ctxgid >> GID_NBITS !=  ctxgidobj_arr[i-1].ctxgid >> GID_NBITS) num_ctx++;
	    if (ctxgidobj_arr[i].ctxgid  == ctxgidobj_arr[i-1].ctxgid) {
			uint32_t gid = (uint32_t)ctxgidobj_arr[i].ctxgid & gidmask;		
			if(num_ctx_cfltobj_arr[gid].num_ctx == (sketch_index[gid+1] - sketch_index[gid])) 
				num_ctx_cfltobj_arr[gid].num_conflictobj++ ; 
				
			num_ctx_cfltobj_arr[gid].num_ctx-- ;
	    } 
	}
//	for(int i = 0 ; i < infile_num; i++) num_conflictobj +=  num_ctx_cfltobj_arr[i].num_conflictobj;    	
	sort_sketch_summary_t *sort_sketch_summary = malloc(sizeof(sort_sketch_summary_t));
	sort_sketch_summary->arrlen = arrlen; sort_sketch_summary->numgids_perctx = (double)arrlen/num_ctx ; 
	sort_sketch_summary->num_ctx = num_ctx; sort_sketch_summary->infile_num = infile_num; 
	sort_sketch_summary->num_ctx_cfltobj_arr = num_ctx_cfltobj_arr;
	return sort_sketch_summary;
}



void free_sort_sketch_summary(sort_sketch_summary_t * sort_sketch_summary){
	free(sort_sketch_summary->num_ctx_cfltobj_arr);
	free(sort_sketch_summary);
}

#include <assert.h>
#define CTX(X, Y)  (ctx[(size_t)(((X)*((X)+1))/2 + (Y))])
#define OBJ(X, Y)  (obj[(size_t)(((X)*((X)+1))/2 + (Y))])

/* use sketch variants to caculate distance */
//1.global sorted comb_sketch64 (i.e. sorted_ctxgidobj_arr or inverted index): 
//slow when dist matrix is large, low memory cache efficient, but may be very fast when matrix is small?
void sorted_ctxgidobj_arr2triangle (ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary){
	uint64_t gidmask = UINT64_MAX >> (64-GID_NBITS);
	uint32_t arrlen = sort_sketch_summary->arrlen;
	int infile_num = sort_sketch_summary->infile_num;
	uint16_t *ctx = calloc( (size_t)infile_num*(infile_num+1)/2 , sizeof(uint16_t));
	uint16_t *obj = calloc( (size_t)infile_num*(infile_num+1)/2 , sizeof(uint16_t));

	for(uint32_t i= 0, j ;i < arrlen - 1 ; i = j ){
		j  = i + 1 ;
		// find range i..j;
		for( ; j< arrlen && (ctxgidobj_arr[i].ctxgid >>GID_NBITS  == ctxgidobj_arr[j].ctxgid >>GID_NBITS ); j++) ;

#pragma omp parallel for num_threads(32) schedule(guided)
		for(uint32_t a = i+1 ;  a < j ; a++ ){
			if (ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a-1].ctxgid || ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a+1].ctxgid ) continue ; // with confclit object
            uint32_t x = ctxgidobj_arr[a].ctxgid & gidmask ;
#pragma omp parallel for num_threads(32) schedule(guided)
			for(uint32_t b = i ;  b < a ; b++ ){
			
				 if ( (b > 0 && ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b-1].ctxgid) || ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b+1].ctxgid ) continue ;			
			     uint32_t y = ctxgidobj_arr[b].ctxgid & gidmask ;
				 CTX(x, y) ++ ;
				if(ctxgidobj_arr[a].obj != ctxgidobj_arr[b].obj) 
					OBJ(x, y)++; 
			}			
		}
printf("\ri=%d",i);				
	}
printf("flg6:OK\n");
	for(int x = 1 ; x <  infile_num; x++){
		 
		for(int y = 0; y< x; y++ ) {
			if( CTX(x, y) > 0) {
				printf("%d\t%d\t%d\t%d\t%f\n",x,y,CTX(x, y), OBJ(x,y), (float) OBJ(x,y)/CTX(x, y));
			}
		}

	}
printf("flg7:OK\n");
		
}
//2. inverted index(i.e. global sorted_ctxgidobj_arr) X common index(i.e. genome-wise sorted comb_sortedsketch64 ) :
// ** the fatest method when dist matrix is sparse
void sorted_ctxgidobj_arrXcomb_sortedsketch64 ( unify_sketch_t* qry_result, ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary ){
	
	uint64_t* qry_sketch_index = qry_result->sketch_index;
	uint64_t* qry_comb_sketch =  qry_result->comb_sketch;
 	int qry_infile_num = qry_result->infile_num;

	uint64_t gidmask = (1UL << Bitslen.gid) - 1 ;
	uint64_t objmask = (1UL << Bitslen.obj) - 1;
    uint32_t ref_arrlen = sort_sketch_summary->arrlen;
    int ref_infile_num = sort_sketch_summary->infile_num;
// for self comparision only
	assert(qry_infile_num == ref_infile_num);

    uint16_t *ctx = calloc( (size_t)ref_infile_num*(ref_infile_num+1)/2 , sizeof(uint16_t));
    uint16_t *obj = calloc( (size_t)ref_infile_num*(ref_infile_num+1)/2 , sizeof(uint16_t));
#pragma omp parallel for num_threads(32) schedule(guided)
	for( uint32_t rn = 0; rn < qry_infile_num; rn++){
//		if(rn > 0) break;
		const uint64_t *a = qry_comb_sketch + qry_sketch_index[rn]; 
		size_t a_size = qry_sketch_index[rn+1] - qry_sketch_index[rn];		
		size_t* idx = find_first_occurrences_AT_ctxgidobj_arr (a, a_size, ctxgidobj_arr , ref_arrlen);
      for(int i = 0; i< a_size; i++ ){
			  if(idx[i]==SIZE_MAX) continue;
			for(int d = idx[i];  ; d++){			
				uint32_t gid = ctxgidobj_arr[d].ctxgid & gidmask;
				if (gid >= rn || (ctxgidobj_arr[d].ctxgid >> Bitslen.gid) != (a[i] >> Bitslen.obj)) break;
				CTX(rn, gid)++;
				if(a[i] & objmask != ctxgidobj_arr[d].obj) OBJ(rn,gid)++;
			}			
		}				
	   free(idx);
   }	

/*
   for( uint32_t rn = 0; rn < qry_infile_num; rn++){
        for(int i = 0 ; i < rn;i++){
            if( CTX(rn,i) > 0 ) printf("\t%d|%d|%d|%f",i, CTX(rn,i),OBJ(rn,i),(float) OBJ(rn,i)/CTX(rn, i));
        }
        printf("\n");
	}
*/
}//end

//3. common index X common index (i.e. genome-wise sorted comb_sortedsketch64 ):
//    code hits: use paris wise small sortted arrays overlapping
/*  Advatages:
		** immediate output (no need precompute dist matrix )	
		** most convience and memory efficience, no invert indxeing needed  
    ** when genome are highly similar (dense dist matrix) speed is even faster than inverted index(i.e. global sorted_ctxgidobj_arr) X common index
		** small sortted arrays overlapping is ~ 2X-3X faster than hashtable lookup based overlapping.
//...
*/

void comb_sortedsketch64Xcomb_sortedsketch64 ( unify_sketch_t* ref_result, unify_sketch_t* qry_result ){
//#pragma omp parallel for num_threads(32) schedule(guided)	
	for(uint32_t rn = 0; rn < ref_result->infile_num; rn++){

		uint64_t *arr_ref = ref_result->comb_sketch + ref_result->sketch_index[rn];
		size_t len_ref = ref_result->sketch_index[rn + 1] - ref_result->sketch_index[rn]; 
		printf("%s",ref_result->gname[rn]);

		for(uint32_t qn = 0; qn < qry_result->infile_num; qn++){
		    uint64_t *arr_qry = qry_result->comb_sketch + qry_result->sketch_index[qn];
    		size_t len_qry = qry_result->sketch_index[qn + 1] - qry_result->sketch_index[qn];			
			  size_t XnY = count_overlaps(arr_ref,len_ref,arr_qry,len_qry);
			  double dist = get_mashD(klen,len_ref,len_qry,XnY);
				printf("\t%lf", 1-dist);
		}
		printf("\n");
	}

}
void check_comb_sortedsketch64 ( unify_sketch_t* result){
	for(uint32_t rn = 0; rn < result->infile_num; rn++){
		uint64_t *arr = result->comb_sketch + result->sketch_index[rn];			
    size_t len = result->sketch_index[rn + 1] - result->sketch_index[rn];
		for(uint32_t i = 1 ; i < len; i++){
			if (arr[i] < arr[i-1] )	err(EXIT_FAILURE,"%s(): %dth genome %dth kmer < %dth kmer (%lx<%lx)", __func__,rn,i,i-1,arr[i],arr[i-1] )	;
		}
	} 
}













//

/**
 * Finds the first occurrence index in large sorted array `b` for each element in a small sorted array `a`.
 * 
 * @param a       Sorted array of uint64_t elements (ascending order)
 * @param a_size  Number of elements in array `a`
 * @param b       Sorted array of uint64_t elements (ascending order)
 * @param b_size  Number of elements in array `b`
 * @return        Array of indices (size_t*) where indices[i] = first occurrence of a[i] in b,
 *                SIZE_MAX if not found. Caller must free() the returned array.
 */
size_t* find_first_occurrences_AT_ctxgidobj_arr (const uint64_t *a, size_t a_size, 
                               const ctxgidobj_t *b, size_t b_size) {
    size_t *indices = malloc(a_size * sizeof(size_t));
    if (!indices) return NULL;

    size_t low = 0;  // Track lower bound for binary search
    
    for (size_t i = 0; i < a_size; ++i) {
		const uint64_t target = a[i]  >> (2*(2*holen + iolen ));

        size_t high = b_size;
        size_t pos = SIZE_MAX;  // Default: not found

        // Leftmost binary search within [low, high)
        while (low < high) {
            size_t mid = low + (high - low) / 2;
            if ((b[mid].ctxgid >> GID_NBITS) < target) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }

        // Check if target was found: modified from: if (low < b_size && b[low] == target) { 
        if (low < b_size && (b[low].ctxgid >> GID_NBITS == target ) ) {
            pos = low;
        }

        indices[i] = pos;
    }

    return indices;
}




