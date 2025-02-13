#ifndef SKETCH_H
#define SKETCH_H
#include "global_basic.h"
#include "command_sketch_wrapper.h"
#include "../klib/khash.h"
#include "../klib/kseq.h"
#include "MurmurHash3.h"
//#include "xxhash.h"
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

inline uint32_t GET_SKETCHING_ID(uint64_t v1, uint64_t v2, uint64_t v3 ,uint64_t v4 , uint64_t v5){
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
/* mv to global
typedef struct dim_sketch_stat
{
  unsigned int hash_id; // sketching type coding
	bool koc;
  int klen; // full length of kmer, 8..31
  int hclen; // half context length,
  int holen; // half outer object length, 1..64
  int drfold; //dimension reduction fold 2^n , 0..32
	int infile_num;
} dim_sketch_stat_t;
*/
// define khash type

void compute_sketch(sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat);
int merge_comblco (sketch_opt_t * sketch_opt_val);
int reads2sketch64 (char* seqfname,char * outfname, bool abundance, int n );
void compute_sketch_splitmfa ( sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat);
//void print_hash_table(khash_t(kmer_hash) *h);


#endif
