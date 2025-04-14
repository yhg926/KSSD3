#ifndef COMMAND_ANI
#define COMMAND_ANI
#include <argp.h>
#include <math.h>
#include <tgmath.h>
#include "global_basic.h"
#include "kssdlib_sort.h"
#include "sketch_rearrange.h"
typedef struct ani_opt
{
	int metric; //0 mashD; 1, aafD; 2. 1-ani
	double c; //minimal distance to enrolled sketches
	int p; //threads
	bool d; //diagnal
	bool ani;
	int e;
	float afcut;
	float anicut;
	char index[PATHLEN];
	char qrydir[PATHLEN];
	char refdir[PATHLEN];
  	char outf[PATHLEN];
	char gl[PATHLEN]; // genome list with selection code 	
  	int num_remaining_args;
  	char **remaining_args;
} ani_opt_t;

typedef struct  {
     uint32_t num_ctx; // including with confilict obj
     uint32_t num_conflictobj; } num_ctx_cfltobj_t;

typedef struct  {
    uint64_t arrlen;
    uint32_t num_ctx;
    uint32_t num_conflictobj;
    double numgids_perctx;
    int infile_num;
    num_ctx_cfltobj_t *num_ctx_cfltobj_arr;
} sort_sketch_summary_t;

/* mv to sketch_rearrange.h
typedef uint32_t obj_t ; // optional uint64_t
typedef struct id_obj{uint32_t gid;obj_t obj;} id_obj_t;
typedef struct { int num; uint32_t *idx ;} inverted_t;
typedef struct __attribute__((packed)) {  uint32_t part[3];} uint96_t;
typedef struct {uint64_t ctxgid; uint32_t obj;} ctxgidobj_t;
typedef struct {uint8_t ctx; uint8_t gid; uint8_t obj;} bitslen_t;
typedef struct {uint64_t ctx; uint32_t gid; uint32_t obj;} tmp_ctxgidobj_t;
*/
int cmd_ani(struct argp_state* state);
//void const_comask_init(dim_sketch_stat_t *lco_stat_val );
int compute_ani(ani_opt_t *ani_opt);
int mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt);
int sparse_mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt);
int compare_u32(const void *a, const void *b);
void sort_arrays(uint64_t *a, uint32_t *b, uint32_t n);
ctxgidobj_t * comb_sortedsketch64_2sortedcomb_ctxgid64obj32(unify_sketch_t* result);
//uint64_t * kmer_arr2regco(uint64_t *array, uint32_t arrlen, int klen, int hclen, int holen);
//void gen_inverted_index4comblco(const char *refdir);
/* mv to sketch_rearrange.h
void sketch64_2ctxobj64 (uint64_t *sketch_index, uint64_t *sketch64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen );
uint96_t* sketch64_2uint96co(uint64_t *sketch_index, uint64_t *sketch64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen );
ctxgidobj_t* ctxobj64_2ctxgidobj(uint64_t *sketch_index, uint64_t *ctxobj64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen );
*/
sort_sketch_summary_t *summarize_ctxgidobj_arr(ctxgidobj_t* ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num);
void free_sort_sketch_summary(sort_sketch_summary_t * sort_sketch_summary);
//void sorted_ctxgidobj_arr2triangle (ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary);
void sorted_ctxgidobj_arrXcomb_sortedsketch64 ( unify_sketch_t* qry_result, ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary );
void comb_sortedsketch64Xcomb_sortedsketch64 ( unify_sketch_t* ref_result, unify_sketch_t* qry_result );

size_t* find_first_occurrences_AT_ctxgidobj_arr (const uint64_t *a, size_t a_size, const ctxgidobj_t *b, size_t b_size);


//print functions
void ani_block_print(int ref_infile_num, int qry_gid_offset, int this_block_size, uint64_t *ref_sketch_index, uint64_t *qry_sketch_index, uint32_t *ctx, uint32_t *obj,char (*refname)[PATHLEN], char (*qryfname)[PATHLEN], FILE *outfp, float af_cutoff, float ani_cutoff );
#endif
