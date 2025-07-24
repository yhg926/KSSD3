#ifndef SKETCH_REARRANGE_H
#define SKETCH_REARRANGE_H

#include "global_basic.h"
#include <math.h>
#include <tgmath.h>
#define GID_NBITS 20
//glovbal public vars
//uint32_t FILTER,hash_id;
extern uint64_t ctxmask,tupmask,ho_mask_len, hc_mask_len, io_mask_len, ho_mask_left,hc_mask_left, io_mask, hc_mask_right, ho_mask_right;
extern uint8_t iolen, klen,hclen,holen ;
extern uint32_t gid_mask; 

typedef uint32_t obj_t ; // optional uint64_t
typedef struct id_obj{uint32_t gid;obj_t obj;} id_obj_t;
typedef struct { int num; uint32_t *idx ;} inverted_t;
typedef struct __attribute__((packed)) {  uint32_t part[3];} uint96_t;
//typedef struct  {uint64_t ctxgid; uint32_t obj;} ctxgidobj_t;
typedef struct  __attribute__((packed)) {uint64_t ctxgid; uint32_t obj;} ctxgidobj_t;
typedef struct {uint8_t ctx; uint8_t gid; uint8_t obj;} bitslen_t;
extern bitslen_t Bitslen;
typedef struct {uint64_t ctx; uint32_t gid; uint32_t obj;} tmp_ctxgidobj_t;

void const_comask_init(dim_sketch_stat_t *lco_stat_val );

ctxgidobj_t *ctxobj64_2ctxgidobj(uint64_t *sketch_index, uint64_t *ctxobj64, int infile_num, uint32_t arrlen);

//other versions of the functions are ignored
//void sketch64_2ctxobj64(uint64_t *sketch64, uint32_t arrlen);
//uint96_t *sketch64_2uint96co(uint64_t *sketch_index, uint64_t *sketch64, int infile_num, uint32_t arrlen);
//uint64_t * kmer_arr2regco(uint64_t *array, uint32_t arrlen, int klen, int hclen, int holen);
//inline declarations:
//inline obj_t extract_outobj_frm_kmer64 (uint64_t kmer, int klen, int holen);
//inline uint96_t concat_lower_bits( tmp_ctxgidobj_t *tmp_ctxgidobj, bitslen_t bitslen);
//inline uint64_t extract_bits(uint96_t num, uint32_t n, uint32_t m);

// core kmer rearrange to ctxobj64
inline obj_t extract_outobj_frm_kmer64 (uint64_t kmer, int klen, int holen){ //holen: one side outter obj len
        return ((kmer >> (2*( klen - holen )) ) << (2*holen)) | ( (kmer <<(64-2*holen) ) >> (64-2*holen));
}

static inline uint64_t uint64_kmer2ctxobj (uint64_t unituple){

 return (      ((unituple & hc_mask_left) << 2*holen)            |
               ((unituple & hc_mask_right) << 2*(holen + iolen)) |
               ((unituple & ho_mask_left) >> 4*hclen )           |
               ((unituple & ho_mask_right) << (2*iolen) )        |
               ((unituple & io_mask) >> (2*(holen + hclen )))    );
}


static inline uint96_t concat_lower_bits( tmp_ctxgidobj_t *tmp_ctxgidobj, bitslen_t bitslen) {

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

static inline uint64_t extract_bits(uint96_t num, uint32_t n, uint32_t m) {
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

static inline uint96_t uint64_kmer2uint96(uint64_t unituple,uint32_t gid){
	tmp_ctxgidobj_t tmp_ctxgidobj;
	tmp_ctxgidobj.obj = ((unituple & ho_mask_left) >> 4*hclen ) | ((unituple & ho_mask_right) << (2*iolen) )|((unituple & io_mask) >> (2*(holen + hclen )))  ; //obj , lowest
    tmp_ctxgidobj.gid = gid;
    tmp_ctxgidobj.ctx = ((unituple & hc_mask_left) >> 2*(holen  + iolen) ) | ((unituple & hc_mask_right) >> (2*holen));
    return concat_lower_bits(&tmp_ctxgidobj, Bitslen);
}

static inline ctxgidobj_t uint64_ctxobj2ctxgidobj96( uint64_t ctxobj64, uint32_t gid,uint64_t obj_len_mask){
	ctxgidobj_t ctxgidobj;
	ctxgidobj.ctxgid = (ctxobj64 >> Bitslen.obj << GID_NBITS) | (gid & gid_mask );
	ctxgidobj.obj = ctxobj64 & obj_len_mask ;
	return ctxgidobj;
}



//coden aware context object pattern 
#ifndef NUM_CODENS
#define NUM_CODENS 9
#endif   
static inline uint64_t generate_coden_pattern64 (){
    uint64_t pattern = 0;
    for (int i = 0; i < NUM_CODENS; ++i) {
        pattern <<= 6;
        pattern |= 0b111100;    
    }
    return pattern;
}
// Reorders a right-aligned pattern: 00_111100 x NUM_CODENS
// Each block = [high4][low2] -> 6 bits
// Pattern occupies lowest (6 * NUM_BLOCKS + 2) bits
static inline uint64_t reorder_unituple_by_coden_pattern64 (uint64_t unituple) {
// need reorder NUM_CODENS of blocks high (to ctx) and NUM_CODENS + 1 of blocks low (to obj)
    uint64_t high = 0, low = unituple & 0x3; 
    #pragma unroll
    for (int i = 0; i < NUM_CODENS; ++i) {
     
        unituple >>= 2;
        high = (high << 4) | unituple & 0xF;
        unituple >>= 4;
        low = (low << 2) | (unituple & 0x3);
    }
    return (high << 2*(NUM_CODENS+1)) | low; // Combine high and low parts
}


typedef uint64_t (*uint64kmer2generic_ctxobj_fn)(uint64_t); 
extern uint64kmer2generic_ctxobj_fn uint64kmer2generic_ctxobj;
void set_uint64kmer2generic_ctxobj(bool is_coden_ctxobj_pattern);









#endif
