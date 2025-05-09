#include "command_ani.h"
#include "command_matrix.h"
#include "global_basic.h"
#include "kssdlib_sort.h"
#include "sketch_rearrange.h"
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
extern const char sorted_comb_ctxgid64obj32[];
extern double C9O7_98[6], C9O7_96[6];
size_t file_size;
/*
int compute_ani(ani_opt_t *ani_opt){
	//initialize 
	dim_sketch_stat_t *dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat),&file_size);
	const_comask_init(dim_sketch_stat);
	// read index
	uint64_t *sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir,idx_sketch_suffix),&file_size);	
	assert(file_size == (dim_sketch_stat->infile_num + 1)*sizeof(sketch_index[0]));
	size_t ref_sketch_size = sketch_index[dim_sketch_stat->infile_num];
	ctxgidobj_t *sortedcomb_ctxgid64obj32 = read_from_file(test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32),&file_size); 
	assert(file_size == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]) );

	unify_sketch_t* result = generic_sketch_parse(ani_opt->refdir);
	const_comask_init(&result->stats.lco_stat_val);

	ctxgidobj_t * sortedcomb_ctxgid64obj32 =  comb_sortedsketch64_2sortedcomb_ctxgid64obj32(result);
printf("FLG3:%lu\n",time(NULL) - second);
	uint64_t sketch_size = result->sketch_index[result->infile_num];
	sort_sketch_summary_t *sort_sketch_summary = summarize_ctxgidobj_arr(sortedcomb_ctxgid64obj32, result->sketch_index, sketch_size, result->infile_num);
printf("FLG4:%lu\n",time(NULL) - second);
    sorted_ctxgidobj_arrXcomb_sortedsketch64 ( result, sortedcomb_ctxgid64obj32, sort_sketch_summary );
 printf("FLG5:%lu\n",time(NULL) - second);
}
*/

//ani linear learning coeffs
double C9O7_98[6] = {-88.96, 3.701e-5, 1.537, 1.324, 8.118, 1.856};
double C9O7_96[6] = {-34.01, -1.225e-4, 5.524, 6.456, 39.88, 4.289};

inline double get_learned_ani (int XnY_ctx, float af_qry, float af_ref, float dist, float ani){
  double learned_ani = 0; double coeffs[6] = {0};
  if(hclen == 9 && holen == 7){
    if (ani >= 98) memcpy(coeffs,C9O7_98, sizeof(C9O7_98));
    else if(ani >= 96) memcpy(coeffs,C9O7_96, sizeof(C9O7_96));
    else return learned_ani;
    learned_ani = coeffs[0] + coeffs[1]*XnY_ctx + coeffs[2]*af_qry + coeffs[3]*af_ref + coeffs[4]*dist + coeffs[5]*ani;
  }
  if(learned_ani > 100 ) learned_ani = 100;
  return learned_ani;
}
//-----------

int compare_idani_desc(const void *a, const void *b) {
    const idani_t *itemA = (const idani_t *)a;
    const idani_t *itemB = (const idani_t *)b;
    return (itemA->ani < itemB->ani) - (itemA->ani > itemB->ani);
}
#define MCTX(L, X, Y)  (ctx[(size_t)((L)*(X) + (Y))])
#define MOBJ(L, X, Y)  (obj[(size_t)((L)*(X) + (Y))])
inline void count_ctx_obj_frm_comb_sketch_section(uint32_t *ctx,uint32_t *obj, ctxgidobj_t *ctxgidobj_arr, size_t ref_sksize, int ref_gnum, int section_gnum, uint64_t *section_sk, uint64_t *section_skidx, uint32_t *num_passid_block, idani_t **sort_idani_block, ani_opt_t *ani_opt){
  uint64_t gidmask = UINT64_MAX >> (64-GID_NBITS), objmask = (1UL << Bitslen.obj) - 1;
#pragma omp parallel for  num_threads(ani_opt->p) schedule(guided)
    for (int i = 0 ; i < section_gnum; i++){
        uint64_t *a = section_sk + (section_skidx[i] - section_skidx[0]);
        size_t a_size = section_skidx[i+1] - section_skidx[i];
				assert(a_size > 0);
        size_t* idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, ctxgidobj_arr, ref_sksize);

        for(int j = 0; j< a_size; j++ ){
            if(idx[j]==SIZE_MAX) continue;
                //skip conlict object;
                if ((j >0) && ( (a[j] >> Bitslen.obj) == (a[j-1] >> Bitslen.obj) ) ) continue;
                if ( (j < a_size -1 ) && ( (a[j] >> Bitslen.obj) == (a[j+1] >> Bitslen.obj) ) ) continue;

                for(int d = idx[j];  ; d++){
                    if ((ctxgidobj_arr[d].ctxgid >> Bitslen.gid) != (a[j] >> Bitslen.obj)) break;
                    uint32_t gid = ctxgidobj_arr[d].ctxgid & gidmask; MCTX(ref_gnum, i, gid)++;
                    if((a[j] & objmask) != ctxgidobj_arr[d].obj) MOBJ(ref_gnum, i, gid)++;
                }
        }
        free(idx);
			  if(ani_opt->fmt == 0) { // print details  >0:matrix no need to sort
      		//sorting ref gid by ani descendingly
        	num_passid_block[i] = 0;
        	for(int j = 0; j< ref_gnum;j++){
          	if( (float)MCTX(ref_gnum,i,j)/a_size < ani_opt->afcut) continue;
          	float dist = (float)MOBJ(ref_gnum,i,j)/MCTX(ref_gnum,i,j) ;
          	float ani = pow((1-dist), (float)1/(2*holen));
          	if(ani <  ani_opt->anicut) continue;
          	sort_idani_block[i][num_passid_block[i]].id = j;
          	sort_idani_block[i][num_passid_block[i]].ani = ani;
          	num_passid_block[i]++;
      	}
      	qsort(sort_idani_block[i], num_passid_block[i], sizeof(idani_t), compare_idani_desc);
			}
    }
}

#define BLOCK_SIZE (4096) // #of qry genomes per batch, for mem_eff handling
int mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt){
	dim_sketch_stat_t *ref_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat),&file_size);
 	int ref_infile_num = ref_dim_sketch_stat->infile_num;
  // read index
  size_t ctxgidobj_arr_fsize;
  uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir,idx_sketch_suffix),&file_size);
  assert(file_size == (ref_infile_num + 1)*sizeof(ref_sketch_index[0]));
  size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
  ctxgidobj_t * sortedcomb_ctxgid64obj32 = read_from_file(test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32),&ctxgidobj_arr_fsize);
  assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]) );

  dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat),&file_size);
  int qry_infile_num = qry_dim_sketch_stat->infile_num;
  assert(qry_dim_sketch_stat->hash_id == ref_dim_sketch_stat->hash_id);
  uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir,idx_sketch_suffix),&file_size);
  size_t qry_sketch_size = qry_sketch_index[qry_infile_num];

  int block_size = BLOCK_SIZE; int offset_gid = 0 ;FILE *fp;
  assert( (fp = fopen(test_get_fullpath(ani_opt->qrydir,combined_sketch_suffix),"rb") ) != NULL) ;
  uint64_t *tmp_ctxobj = malloc(qry_sketch_size * sizeof(uint64_t));
  uint32_t *ctx = malloc( ref_infile_num * block_size * sizeof(uint32_t))  ;
  uint32_t *obj = malloc( ref_infile_num * block_size * sizeof(uint32_t))  ;
//for order id by descending ani
  uint32_t *num_passid_block = malloc(block_size * sizeof(uint32_t));
  idani_t **sort_idani_block = malloc(block_size * sizeof(idani_t*));
  for(int i = 0; i< block_size; i++) sort_idani_block[i] = malloc( ref_infile_num * sizeof(idani_t));

  char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_dim_sketch_stat + 1);
  char (*qryname)[PATHLEN] = (char (*)[PATHLEN])(qry_dim_sketch_stat + 1);

  FILE *outfp = ani_opt->outf[0]=='\0' ? stdout: fopen( ani_opt->outf, "w");
//printf header
	if(ani_opt->fmt) { // matrix format
		for(int i = 0;i < ref_infile_num;i++) fprintf(outfp, "\t%s",refname[i]);
		fprintf(outfp,"\n");
	}
	else fprintf(outfp,"Qry\tRef\tXnY_ctx\tQry_align_fraction\tRef_align_fraction\tco-distance\tANI\tlearned_ANI(if > 0)\n");

  for( int b = 0; b <= qry_infile_num/block_size ; b++){

    int this_block_size = (b == qry_infile_num/block_size) ? (qry_infile_num % block_size):block_size;
    int this_sketch_size = qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid];
    int read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t),this_sketch_size,fp);
    uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
    assert(this_sketch_size == read_sketch_size);

    memset(ctx,0,ref_infile_num * block_size * sizeof(uint32_t));
    memset(obj,0,ref_infile_num * block_size * sizeof(uint32_t));
    count_ctx_obj_frm_comb_sketch_section(ctx,obj,sortedcomb_ctxgid64obj32, ref_sketch_size, ref_infile_num, this_block_size, tmp_ctxobj, this_sketch_index, num_passid_block,sort_idani_block,ani_opt);
    if(ani_opt->fmt)  ani_block_print_matrix (ref_infile_num,  b*block_size,  this_block_size, qry_sketch_index, ctx, obj,qryname,outfp, ani_opt);
		else ani_block_print(ref_infile_num, b*block_size,this_block_size,ref_sketch_index,qry_sketch_index,ctx,obj, refname,qryname,num_passid_block,sort_idani_block,outfp);

    offset_gid+=this_block_size;
  }

	for(int i = 0 ; i <block_size;i++) free(sort_idani_block[i]);
  free_all(ref_dim_sketch_stat,ref_sketch_index,qry_dim_sketch_stat,qry_sketch_index,tmp_ctxobj,ctx,obj,num_passid_block,sort_idani_block,NULL);
  free_read_from_file(sortedcomb_ctxgid64obj32,ctxgidobj_arr_fsize);
  fclose(fp); fclose(outfp);
  return ctxgidobj_arr_fsize;
}

void ani_block_print(int ref_infile_num, int qry_gid_offset, int this_block_size, uint64_t *ref_sketch_index, uint64_t *qry_sketch_index, uint32_t *ctx, uint32_t *obj,char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],uint32_t *num_passid_block, idani_t **sort_idani_block, FILE *outfp ){

  for(int i=0; i< this_block_size; i++){
    int qry_gid = qry_gid_offset + i;
    int qry_sketch_size = qry_sketch_index[qry_gid+1] - qry_sketch_index[qry_gid];

    for(int n = 0; n < num_passid_block[i] ;n++){
			int j = sort_idani_block[i][n].id;			
      int XnY_ctx = MCTX(ref_infile_num,i,j);
      float af_qry = (float)XnY_ctx / qry_sketch_size ;
      int ref_sketch_size = ref_sketch_index[j+1] - ref_sketch_index[j];
      float af_ref = (float)XnY_ctx / ref_sketch_size ;
      double dist = (double)MOBJ(ref_infile_num,i,j)/XnY_ctx;
      double ani = sort_idani_block[i][n].ani * 100;
			double learned_ani = get_learned_ani(XnY_ctx, af_qry, af_ref, dist,ani);
      fprintf(outfp,"%s\t%s\t%d\t%f\t%f\t%lf\t%lf\t%lf\n",qryfname[qry_gid],refname[j], XnY_ctx,af_qry,af_ref,dist,ani,learned_ani);
    }
  }
}

void ani_block_print_matrix (int ref_infile_num, int qry_gid_offset, int this_block_size, uint64_t *qry_sketch_index, uint32_t *ctx, uint32_t *obj,char (*qryfname)[PATHLEN],FILE *outfp,ani_opt_t *ani_opt){
	for(int i=0; i< this_block_size; i++){
		int qry_gid = qry_gid_offset + i;
    int qry_sketch_size = qry_sketch_index[qry_gid+1] - qry_sketch_index[qry_gid];		
		fprintf(outfp,"%s", qryfname[qry_gid]);

		for(int j = 0; j < ref_infile_num ; j++){
			double ani; int XnY_ctx = MCTX(ref_infile_num,i,j);
			float af_qry = (float)XnY_ctx / qry_sketch_size;
			if ( af_qry < ani_opt->afcut) ani = ani_opt->e;
			else{
				double dist = (double)MOBJ(ref_infile_num,i,j)/XnY_ctx;
				ani = pow((1-dist), (float)1/(2*holen)) * 100;	
				if ( af_qry < ani_opt->ani) ani = ani_opt->e;
			}
			fprintf(outfp,"\t%lf", ani);				
		}
		fprintf(outfp,"\n");
	}
}


//
#define DIFF_OBJ_BITS 1
size_t dedup_with_ctxobj_counts(uint32_t *arr, size_t n, co_distance_t **ctxobj_cnt) {
    *ctxobj_cnt = NULL;
    if (n == 0) return 0;

    co_distance_t *tmp_ctxobj_cnt = malloc(n * sizeof(co_distance_t));
    if (!tmp_ctxobj_cnt) err(EXIT_FAILURE,"%s(): tmp_ctxobj_cnt malloc failure",__func__); 
    
    size_t j = 0;
    tmp_ctxobj_cnt[0].ctx_ct = 1;
    tmp_ctxobj_cnt[0].diff_obj = arr[0] % 2;
    arr[0] >>= DIFF_OBJ_BITS;

    for (size_t i = 1; i < n; i++) {
        if ((arr[i] >> DIFF_OBJ_BITS) == arr[j]) {
            tmp_ctxobj_cnt[j].ctx_ct++;
            tmp_ctxobj_cnt[j].diff_obj += (arr[i] % 2);
        } else {
            j++;
            tmp_ctxobj_cnt[j].ctx_ct = 1;
            tmp_ctxobj_cnt[j].diff_obj = arr[i] % 2;
            arr[j] = arr[i] >> DIFF_OBJ_BITS;
        }
    }
    // Trim arrays
    co_distance_t *ctxobj_tmp = realloc(tmp_ctxobj_cnt, (j + 1) * sizeof(co_distance_t));
    *ctxobj_cnt = ctxobj_tmp ? ctxobj_tmp : tmp_ctxobj_cnt;
    
    return j + 1;
}

ctxgidobj_t * comb_sortedsketch64_2sortedcomb_ctxgid64obj32(unify_sketch_t* ref_result){
	//const_comask_init(&ref_result->stats.lco_stat_val);
	uint64_t sketch_size = ref_result->sketch_index[ref_result->infile_num] ;
    if (sketch_size > (float) UINT32_MAX * LD_FCTR ) err(EXIT_FAILURE,"%s():sketch_index maximun %lu exceed UINT32_MAX*LF;%f",__func__, sketch_size,(float) UINT32_MAX * LD_FCTR);
    if( ref_result->infile_num >= ( 1<<GID_NBITS )) err(EXIT_FAILURE,"%s(): genome numer %d exceed maximum:%u",__func__,ref_result->infile_num,1<<GID_NBITS);
    if(GID_NBITS + 4*hclen >64) err(EXIT_FAILURE,"%s(): context_bits_len(%d)+gid_bits_len(%d) exceed 64",__func__,4*hclen,GID_NBITS);
	ctxgidobj_t *ctxgidobj_arr = ctxobj64_2ctxgidobj(ref_result->sketch_index, ref_result->comb_sketch, ref_result->infile_num, sketch_size, klen,hclen,holen);
	ctxgidobj_sort_array(ctxgidobj_arr, sketch_size) ;
	return ctxgidobj_arr;		
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


// sparse_mem_eff.. seems slower than mem_eff,
int sparse_mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt){
    //initialize
    dim_sketch_stat_t *ref_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat),&file_size);
    const_comask_init(ref_dim_sketch_stat);

	uint64_t gidmask = UINT64_MAX >> (64-GID_NBITS);
	uint64_t objmask = (1UL << Bitslen.obj) - 1;
 
	int ref_infile_num = ref_dim_sketch_stat->infile_num;
    // read index
    size_t ctxgidobj_arr_fsize;
    uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir,idx_sketch_suffix),&file_size);
    assert(file_size == (ref_infile_num + 1)*sizeof(ref_sketch_index[0]));
    size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
    ctxgidobj_t * sortedcomb_ctxgid64obj32 = read_from_file(test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32),&ctxgidobj_arr_fsize);
    assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]) );

    dim_sketch_stat_t *qry_dim_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat),&file_size);
    int qry_infile_num = qry_dim_sketch_stat->infile_num;
    assert(qry_dim_sketch_stat->hash_id == ref_dim_sketch_stat->hash_id);
    uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir,idx_sketch_suffix),&file_size);
    size_t qry_sketch_size = qry_sketch_index[qry_infile_num];

    int block_size = BLOCK_SIZE; int offset_gid = 0 ;FILE *fp;
    assert( (fp = fopen(test_get_fullpath(ani_opt->qrydir,combined_sketch_suffix),"rb") ) != NULL) ;
    uint64_t *tmp_ctxobj = malloc(qry_sketch_size * sizeof(uint64_t));
	// arry of dynamic gids array per query, note gid <<= DIFF_OBJ_BITS + (ref_obj==qry_obj?0:1);  >>sparse only
//    Vector *ref_gids_perqry_arr  = malloc( block_size * sizeof(Vector));
	//ref_gids_perqry_arr[i][0] is capacity, ref_gids_perqry_arr[i][1] is length of the inner array 
	uint32_t **ref_gids_perqry_arr = malloc( block_size * sizeof(uint32_t *)); 
	co_distance_t ** ctxobj_cnt_perqry_arr = malloc( block_size * sizeof(co_distance_t *)); 
	size_t *lens  = malloc( block_size * sizeof(size_t)) ;     
	for(int i =0 ; i< block_size; i++) {
		ref_gids_perqry_arr[i] = malloc( 6 * sizeof(uint32_t)) ; // 4 for innitialized capacity
		ref_gids_perqry_arr[i][0] = 4;
		ref_gids_perqry_arr[i][1] = 0;
	}
	
    for( int b = 0; b <= qry_infile_num/block_size ; b++){

        int this_block_size = (b == qry_infile_num/block_size) ? (qry_infile_num % block_size):block_size;
        int this_sketch_size = qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid];
        int read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t),this_sketch_size,fp);
        uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
        assert(this_sketch_size == read_sketch_size);

#pragma omp parallel for  num_threads(ani_opt->p) schedule(guided)
    	for (int i = 0 ; i < this_block_size; i++){
        	uint64_t *a = tmp_ctxobj + (this_sketch_index[i] - this_sketch_index[0]);
        	size_t a_size = this_sketch_index[i+1] - this_sketch_index[i];
        	size_t* idx = find_first_occurrences_AT_ctxgidobj_arr (a, a_size, sortedcomb_ctxgid64obj32, ref_sketch_size);

        	for(int j = 0; j< a_size; j++ ){
            	if(idx[j]==SIZE_MAX) continue;
                //skip conlict object;
                	if ((j >0) && ( (a[j] >> Bitslen.obj) == (a[j-1] >> Bitslen.obj) ) ) continue;
                	if ( (j < a_size -1 ) && ( (a[j] >> Bitslen.obj) == (a[j+1] >> Bitslen.obj) ) ) continue;

	                for(int d = idx[j];  ; d++){
                    	if ((sortedcomb_ctxgid64obj32[d].ctxgid >> Bitslen.gid) != (a[j] >> Bitslen.obj)) break;
                    	uint32_t gid01 = sortedcomb_ctxgid64obj32[d].ctxgid & gidmask <<  DIFF_OBJ_BITS;
						if((a[j] & objmask) != sortedcomb_ctxgid64obj32[d].obj) gid01 |= 1;
						if(ref_gids_perqry_arr[i][1] + 2 == ref_gids_perqry_arr[i][0] ){
							ref_gids_perqry_arr[i][0] += 100 ;
							ref_gids_perqry_arr[i] = realloc(ref_gids_perqry_arr[i],ref_gids_perqry_arr[i][0] * sizeof(ref_gids_perqry_arr[i][0]));
						}
						ref_gids_perqry_arr[i][2 + ref_gids_perqry_arr[i][1] ] = gid01;
						ref_gids_perqry_arr[i][1]++;
						//vector_push(&ref_gids_perqry_arr[i],&gid01);						
                	}
        	}
	        free(idx);
			qsort(ref_gids_perqry_arr[i]+2, ref_gids_perqry_arr[i][1], sizeof(uint32_t), qsort_comparator_uint32);
			lens[i] = dedup_with_ctxobj_counts( ref_gids_perqry_arr[i]+2, ref_gids_perqry_arr[i][1], &ctxobj_cnt_perqry_arr[i]);       	 
    	}

        for ( int i = 0 ; i< this_block_size; i++ ){
            int qry_gid = b*block_size+i;
            int qry_sketch_size = ref_sketch_index[qry_gid+1] - ref_sketch_index[qry_gid];
//sparse only code >>
            for (int l=0;l< lens[i]; l++){
				double dist = ctxobj_cnt_perqry_arr[i][l].diff_obj/ctxobj_cnt_perqry_arr[i][l].ctx_ct;
				double ani = pow((1-dist), (double)1/(2*holen));
#define VEC_GET_AS(type, vec, idx) (((type*)(vec).data)[idx])
	//			if(ctxobj_cnt_perqry_arr[i][l].ctx_ct > 20)
      //          	printf ("%d|%d|%d|%d|%d|%lf|%lf\t",qry_gid,VEC_GET_AS(uint32_t, ref_gids_perqry_arr[i], l),qry_sketch_size,ctxobj_cnt_perqry_arr[i][l].ctx_ct,ctxobj_cnt_perqry_arr[i][l].diff_obj,dist,ani);
            }
    //        printf("\n");
//<<  

        }
	//	for(int i =0 ; i< this_block_size; i++) vector_free(&ref_gids_perqry_arr[i]);//ref_gids_perqry_arr[i].size = 0;
#pragma omp parallel for  num_threads(ani_opt->p) schedule(guided)
  		for(int i =0 ; i< this_block_size; i++) {
        	ref_gids_perqry_arr[i][1] = 0 ;
//			ref_gids_perqry_arr[i][0] = 4;
        	free(ctxobj_cnt_perqry_arr[i]);
    	}
	    offset_gid+=this_block_size;
    }// all blocks loop end

	for(int i = 0; i< block_size; i++) free(ref_gids_perqry_arr[i]);
	free_all(ref_gids_perqry_arr,ctxobj_cnt_perqry_arr,NULL);////sparse only code >>
    free_all(ref_dim_sketch_stat,ref_sketch_index,qry_dim_sketch_stat,qry_sketch_index,tmp_ctxobj,NULL);
    free_read_from_file(sortedcomb_ctxgid64obj32,ctxgidobj_arr_fsize);
    fclose(fp);

    return ctxgidobj_arr_fsize;
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
		uint64_t *a = qry_comb_sketch + qry_sketch_index[rn]; 
		size_t a_size = qry_sketch_index[rn+1] - qry_sketch_index[rn];		
		size_t* idx = find_first_occurrences_AT_ctxgidobj_arr (a, a_size, ctxgidobj_arr , ref_arrlen);
      for(int i = 0; i< a_size; i++ ){
			if(idx[i]==SIZE_MAX) continue;
			//skip conlict object;
			if ((i >0) && ( (a[i] >> Bitslen.obj) == (a[i-1] >> Bitslen.obj) ) ) continue;
            if ( (i < a_size -1 ) && ( (a[i] >> Bitslen.obj) == (a[i+1] >> Bitslen.obj) ) ) continue;

			for(int d = idx[i];  ; d++){			
				uint32_t gid = ctxgidobj_arr[d].ctxgid & gidmask;
				if (gid > rn || (ctxgidobj_arr[d].ctxgid >> Bitslen.gid) != (a[i] >> Bitslen.obj)) break;
				CTX(rn, gid)++;
				if((a[i] & objmask) != ctxgidobj_arr[d].obj) //wrong: if(a[i] & objmask != ctxgidobj_arr[d].obj) ...
					OBJ(rn,gid)++;
				
			}			
		}				
	   free(idx);
/*
#pragma omp critical 
		{
	     	for(int i = 0 ; i <= rn;i++){
            	if( CTX(rn,i) > 0 ) printf("\t%d|%d|%d|%d|%f",i,rn, CTX(rn,i),OBJ(rn,i),(float) OBJ(rn,i)/CTX(rn, i));
        	}
			printf("\n");
		}
*/
   }	


   for( uint32_t rn = 0; rn < qry_infile_num; rn++){
		printf("%u",rn);
        for(int i = 0 ; i < rn;i++){
            if( CTX(rn,i) > 0 ) printf("\t%d|%d|%d|%f",i, CTX(rn,i),OBJ(rn,i),(float) OBJ(rn,i)/CTX(rn, i));
        }
        printf("\n");
	}

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
		const uint64_t target = a[i]  >> Bitslen.obj; // (2*(2*holen + iolen ));

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




