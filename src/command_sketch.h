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

/*
#define GID_NBITS 20
//glovbal public vars 
uint32_t FILTER,hash_id;
uint64_t ctxmask,tupmask,ho_mask_len, hc_mask_len, io_mask_len, ho_mask_left,hc_mask_left, io_mask, hc_mask_right, ho_mask_right;
uint8_t iolen, klen,hclen,holen ;
uint32_t gid_mask; bitslen_t Bitslen;// context, gid, and obj bits len
*/
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
#define SKETCH_HASH(key) ({          \
    XXH64_hash_t hash64 = XXH64(&key, 8, 42);   \
   	(uint32_t) hash64; \
})

/* abolished 
#define SKETCH_HASH(key) ({          \
    XXH64_hash_t hash64 = XXH3_64bits_dispatch(&key, 8 );   \
    (uint32_t)(hash64 >> 32) ^ (uint32_t)hash64; \
})
//
#define GET_SKETCHING_ID(inv1,v2,v3,v4,v5) ( { \
	uint64_t v1 = v1, test_num = 31415926; \
  SKETCH_HASH(v1) ^ SKETCH_HASH(v2) ^ SKETCH_HASH(v3) \
  ^ SKETCH_HASH(v4) ^ SKETCH_HASH(v5) ^ SKETCH_HASH(test_num) ; \
})
*/

static inline uint32_t GET_SKETCHING_ID(uint64_t v1, uint64_t v2, uint64_t v3 ,uint64_t v4 , uint64_t v5){
	uint64_t test_num = 31415926;
	return SKETCH_HASH(v1) ^ SKETCH_HASH(v2) ^ SKETCH_HASH(v3) ^ SKETCH_HASH(v4) ^ SKETCH_HASH(v5) ^ SKETCH_HASH(test_num) ; 	
}

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
// define khash type

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

#endif
