#ifndef COMMAND_ANI
#define COMMAND_ANI
#include <argp.h>
#include <math.h>
#include <tgmath.h>
#include "global_basic.h"

typedef struct ani_opt
{
	int metric; //0 mashD; 1, aafD; 2. 1-ani
	double c; //minimal distance to enrolled sketches
	int p; //threads
	bool d; //diagnal
	bool ani;
	int e;
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

typedef uint32_t obj_t ; // optional uint64_t
typedef struct id_obj{uint32_t gid;obj_t obj;} id_obj_t;
typedef struct { int num; uint32_t *idx ;} inverted_t;
typedef struct __attribute__((packed)) {  uint32_t part[3];} uint96_t;

typedef struct {uint64_t ctxgid; uint32_t obj;} ctxgidobj_t;
typedef struct {uint8_t ctx; uint8_t gid; uint8_t obj;} bitslen_t;
typedef struct {uint64_t ctx; uint32_t gid; uint32_t obj;} tmp_ctxgidobj_t;

int cmd_ani(struct argp_state* state);
void const_comask_init(dim_sketch_stat_t *lco_stat_val );
int compute_ani(ani_opt_t *ani_opt);
int compare_u32(const void *a, const void *b);
void sort_arrays(uint64_t *a, uint32_t *b, uint32_t n);
uint64_t * kmer_arr2regco(uint64_t *array, uint32_t arrlen, int klen, int hclen, int holen);
void gen_inverted_index4comblco(const char *refdir);

void sketch64_2ctxobj64 (uint64_t *sketch_index, uint64_t *sketch64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen );
uint96_t* sketch64_2uint96co(uint64_t *sketch_index, uint64_t *sketch64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen );

ctxgidobj_t* ctxobj64_2ctxgidobj(uint64_t *sketch_index, uint64_t *ctxobj64, int infile_num,  uint32_t arrlen, int klen, int hclen, int holen );
sort_sketch_summary_t *summarize_ctxgidobj_arr(ctxgidobj_t* ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num);
void free_sort_sketch_summary(sort_sketch_summary_t * sort_sketch_summary);
//void sorted_ctxgidobj_arr2triangle (ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary);
void sorted_ctxgidobj_arrXcomb_sortedsketch64 ( unify_sketch_t* qry_result, ctxgidobj_t* ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary );
void comb_sortedsketch64Xcomb_sortedsketch64 ( unify_sketch_t* ref_result, unify_sketch_t* qry_result );

size_t* find_first_occurrences_AT_ctxgidobj_arr (const uint64_t *a, size_t a_size, const ctxgidobj_t *b, size_t b_size);

//inline functions
inline obj_t extract_outobj_frm_kmer64 (uint64_t kmer, int klen, int holen){ //holen: one side outter obj len
        return ((kmer >> (2*( klen - holen )) ) << (2*holen)) | ( (kmer <<(64-2*holen) ) >> (64-2*holen));
}

inline uint96_t concat_lower_bits( tmp_ctxgidobj_t *tmp_ctxgidobj, bitslen_t bitslen) {

    // 2. Calculate concatenated value (using 128-bit arithmetic)
    __uint128_t concat = 0;
    concat |= (__uint128_t)tmp_ctxgidobj->ctx << (bitslen.gid + bitslen.obj);  // MSB part
    concat |= (__uint128_t)tmp_ctxgidobj->gid << bitslen.obj;             // Middle part
    concat |= tmp_ctxgidobj->obj;                                    // LSB part

    // 3. Extract lower 96 bits
    uint96_t result;
    result.part[0] = (concat >> 64) & 0xFFFFFFFF;  // Highest 32 bits of 96
    result.part[1] = (concat >> 32) & 0xFFFFFFFF;  // Middle 32 bits
    result.part[2] = concat & 0xFFFFFFFF;          // Lowest 32 bits
    return result;
}

inline uint64_t extract_bits(uint96_t num, uint32_t n, uint32_t m) {
    // Combine parts into 128-bit temporary for bit manipulation
    __uint128_t value = ((__uint128_t)num.part[0] << 64) | 
                        ((__uint128_t)num.part[1] << 32) | 
                        num.part[2];

    // Calculate bit range and validate inputs
    const uint32_t total_bits = m - n + 1;
    const uint32_t shift_amount = n;
    
    // Shift and mask to extract desired bits
    value >>= shift_amount;
    const uint64_t mask = (total_bits >= 64) ? 
                          UINT64_MAX : 
                          ((1ULL << total_bits) - 1);

    return (uint64_t)(value & mask);
}


//int ani_triangle(ani_opt_t *);
//int ani_matrix(ani_opt_t *);
#endif
