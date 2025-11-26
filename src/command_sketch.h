#ifndef SKETCH_H
#define SKETCH_H
#include "global_basic.h"
#include "sketch_rearrange.h"
#include "command_sketch_wrapper.h"
#include "../klib/khash.h"
#include "../klib/kseq.h"
#include "MurmurHash3.h"
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h> // open function
#include <unistd.h>
#include <stdint.h>
#include <stdio.h>

#ifndef XXHASH_H
#define XXHASH_H
#define XXH_STATIC_LINKING_ONLY /* access advanced declarations */
#define XXH_IMPLEMENTATION      /* access definitions */
#include "xxhash.h"
#include "xxh_x86dispatch.h"

#endif

// hash functions for sketching
// hash 1 : murmur3
/*
#define SKETCH_HASH(key) ({          \
    uint64_t out[2];                        \
    MurmurHash3_x64_128(&key, 8, 42, out);   \
    (uint32_t)(out[0] ^ out[1]); \
})
*/
//hash 2 : xxhash
/*
#define SKETCH_HASH(key) ({          \
    XXH64_hash_t hash64 = XXH64(&key, 8, 42);   \
   	(uint32_t) hash64; \
})
*/
//hash 3:
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
  key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
};
// hash 4: Cheap but high-quality 64-bit mix (SplitMix64 finalizer).
static inline uint64_t mix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
};

#define SKETCH_HASH(key) ({   (uint32_t) mix64(key); })

static inline uint32_t GET_SKETCHING_ID(uint64_t v1, uint64_t v2, uint64_t v3 ,uint64_t v4 , uint64_t v5){
	uint64_t test_num = 31415926;
	return SKETCH_HASH(v1) ^ SKETCH_HASH(v2) ^ SKETCH_HASH(v3) ^ SKETCH_HASH(v4) ^ SKETCH_HASH(v5) ^ SKETCH_HASH(test_num) ; 	
}


// ----- make_ctxobj_fast.h ----------------------------------------------------
#pragma once

#if defined(__x86_64__)
  #include <immintrin.h>   // _pext_u64
#endif

// ======================= FASTEST: BMI2 path (if available) ===================
// Compile this function for BMI2 only; caller doesn’t need -mbmi2 globally.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("bmi2"), always_inline))
static inline uint64_t make_ctxobj_bmi2(uint64_t unituple, uint64_t tuplemask, uint64_t ctxmask, uint8_t  nobjbits)
{
    const uint64_t ctx = _pext_u64(unituple, ctxmask);
    const uint64_t obj = _pext_u64(unituple, tuplemask & ~ctxmask);
    // Ensure object occupies exactly the low nobjbits (cheap mask)
    return (ctx << nobjbits) | (obj);
}
#endif

// ====================== BASELINE: portable fallback path ======================
__attribute__((always_inline))
static inline uint64_t make_ctxobj_baseline(uint64_t unituple){
#if NUM_CODENS > 10
        printf("warning: %s() initialization failed due to NUM_CODENS(%d) is out of range (NUM_CODENS <=11) ", __func__, NUM_CODENS);
#endif
    return uint64kmer2generic_ctxobj(unituple );
}

// ========================= PUBLIC ENTRY (fast + portable) ====================
__attribute__((always_inline))
static inline uint64_t make_ctxobj(uint64_t unituple, uint64_t tuplemask, int64_t ctxmask, uint8_t  nobjbits)
{
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    // Cache the CPU feature check (no per-call overhead).
    static int inited = 0;
    static int has_bmi2 = 0;
    if (!inited) {
        has_bmi2 = __builtin_cpu_supports("bmi2");
        __atomic_store_n(&inited, 1, __ATOMIC_RELAXED);
    }
    if (has_bmi2) {
        return make_ctxobj_bmi2(unituple, tuplemask, ctxmask, nobjbits);
    }
#endif
    return make_ctxobj_baseline(unituple);
}












//
//initialize funs
//void public_vars_init(dim_sketch_stat_t* sketch_stat_raw) ; //initla global vars from sketch subcommand pars before sketch generated

void compute_sketch(sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat);
void combine_lco( sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat);
int merge_comblco (sketch_opt_t * sketch_opt_val);
void gen_inverted_index4comblco(const char* sketchdir);
//sketchuing methods family
// produce sorted sketch
void read_genomes2mem2sortedctxobj64 (sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat, int batch_size);
int reads2sketch64 (char* seqfname,char * outfname, bool abundance, int n );
int seq2ht_sortedctxobj64 (char* seqfname, char * outfname, bool abundance, int n ) ;
int seq2sortedsketch64(char* seqfname,char * outfname, bool abundance, int n );
int opt_seq2sortedsketch64(char* seqfname, char * outfname, bool abundance, int n) ;
int opt2_seq2sortedsketch64(char* seqfname, char * outfname, bool abundance, int n) ; // no SIMD optimization
// produce sketch without sort
//void compute_sketch_splitmfa ( sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat);
void mfa2sortedctxobj64( sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat);
//void print_hash_table(khash_t(kmer_hash) *h);
void write_sketch_stat (const char* outdir, infile_tab_t* infile_stat);
simple_sketch_t* simple_genomes2mem2sortedctxobj64_mem (infile_tab_t *infile_stat, int drfold);
static void sketch_many_files_in_parallel(sketch_opt_t *opt, infile_tab_t *tab, int batch_size);
static void sketch_few_files_with_intrafile_parallel(sketch_opt_t *opt, infile_tab_t *tab, int BATCH_READS);
void mfa2sortedctxobj64_v2 (sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat);



#endif
