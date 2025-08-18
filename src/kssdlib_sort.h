#ifndef KSSDLIB_SORT_H
#define KSSDLIB_SORT_H
#include "command_ani.h"
#include "sketch_rearrange.h"
/*0. qsort comparator */
int qsort_comparator_uint32 (const void *a, const void *b);
int qsort_comparator_uint64 (const void *a, const void *b);
/*1. count overlap uint64_t a[N] and b[M]*/
size_t count_overlaps_two_pointers(const uint64_t *a, size_t n, const uint64_t *b, size_t m);
bool binary_search(const uint64_t *arr, size_t size, uint64_t target) ;
size_t count_overlaps_binary_search(const uint64_t *small, size_t small_size, const uint64_t *large, size_t large_size);
size_t count_overlaps(const uint64_t *a, size_t n, const uint64_t *b, size_t m) ;

/*2.custom sort*/
static void custom_sort(uint96_t *arr, size_t n); //void custom_sort(co_t *arr, size_t n);
void sort_uint96_array(uint96_t *arr, size_t n) ;
void ctxgidobj_sort_array(ctxgidobj_t *arr, size_t n) ;


/*3. khash sort*/
#include "../klib/khash.h"
KHASH_MAP_INIT_INT64(sort64, uint32_t)
typedef struct { uint64_t key;  uint32_t value;} kv_pair_t;
typedef struct { uint64_t *keys; uint32_t *values; size_t len;} SortedKV_Arrays_t ;
int cmp_kv_pair(const void *a, const void *b);
SortedKV_Arrays_t sort_khash_u64 (khash_t(sort64) *h);
SortedKV_Arrays_t gpt_sort_khash_u64(khash_t(sort64) *h);
void filter_n_SortedKV_Arrays(SortedKV_Arrays_t *result, uint32_t n);
void remove_ctx_with_conflict_obj (SortedKV_Arrays_t *data, uint32_t n_obj_bits);
void remove_ctx_with_conflict_obj_noabund (uint64_t *keys, size_t *len, uint32_t n_obj_bits) ;
//void free_SortedKV_Arrays (SortedKV_Arrays_t *result);


/* 4. related methods for already sorted array  */
//4.1. sorted array deduplicate 
static inline size_t dedup_sorted_uint64(uint64_t *arr, size_t n) {
    if (n <= 1) return n;  // Handle empty/single-element cases

    size_t j = 0;  // Position for next unique element
    
    for (size_t i = 1; i < n; i++) {
        if (arr[i] != arr[j]) {
            arr[++j] = arr[i];  // Shift unique element forward
        }
    }
    return j + 1;  // New array length
}

/**
 * 4.2. Deduplicates sorted uint64 array and returns occurrence counts
 * 
 * @param arr     Sorted array (modified in-place)
 * @param n       Number of elements in input array
 * @param counts  Output parameter for occurrence counts array
 * @return        New length of deduplicated array
 */
static inline size_t dedup_with_counts(uint64_t *arr, size_t n, uint32_t **counts) {
    if (n == 0) {
        *counts = NULL;
        return 0;
    }
    // Allocate maximum possible size for counts
    uint32_t *cnt = malloc(n * sizeof(uint32_t));
    if (!cnt) {
        *counts = NULL;
        // Fallback: deduplicate without counts
        return dedup_sorted_uint64(arr,n);
    }
    size_t j = 0;
    cnt[j] = 1;
    // Single pass with count tracking
    for (size_t i = 1; i < n; i++) {
        if (arr[i] == arr[j]) {
            cnt[j]++;
        } else {
            arr[++j] = arr[i];
            cnt[j] = 1;
        }
    }
    // Trim counts array to actual size
    uint32_t *tmp = realloc(cnt, (j + 1) * sizeof(uint32_t));
    *counts = tmp ? tmp : cnt;  // Keep original if realloc fails
    return j + 1;
}
/* Usage example: 
  int main() {
    uint64_t arr[] = {1, 1, 2, 2, 2, 3, 4, 4, 4, 4};
    size_t n = sizeof(arr)/sizeof(arr[0]);
    uint32_t *counts;

    size_t new_len = dedup_with_counts(arr, n, &counts);
    
    // arr[] = {1, 2, 3, 4} 
    // counts = {2, 3, 1, 4}
    // new_len = 4

    free(counts);
    return 0;
}
*/

// 8-pass stable LSD radix sort on u64
static inline void radix_sort_u64(uint64_t *keys, size_t n){
    if (n < 2) return;
    uint64_t *tmp = (uint64_t*)malloc(n*sizeof(uint64_t));
    size_t cnt[256];
    for (unsigned pass=0; pass<8; ++pass){
        for (int i=0;i<256;++i) cnt[i]=0;
        unsigned sh = pass*8;
        for (size_t i=0;i<n;++i) ++cnt[(keys[i]>>sh)&0xFFu];
        size_t sum=0;
        for (int i=0;i<256;++i){ size_t c=cnt[i]; cnt[i]=sum; sum+=c; }
        for (size_t i=0;i<n;++i) tmp[cnt[(keys[i]>>sh)&0xFFu]++] = keys[i];
        uint64_t *sw = keys; keys = tmp; tmp = sw;
    }
    free(tmp);
}

// ---- radix sort (keys+counts in lockstep) + in-place merge-equals (sum counts) ----

static inline void radix_sort_kv_u64(uint64_t *k, uint32_t *v, size_t n){
    if (n < 2) return;
    uint64_t *tk = (uint64_t*)malloc(n*sizeof(uint64_t));
    uint32_t *tv = (uint32_t*)malloc(n*sizeof(uint32_t));
    size_t cnt[256];
    for (unsigned pass=0; pass<8; ++pass){
        for (int i=0;i<256;++i) cnt[i]=0;
        unsigned sh = pass*8;
        for (size_t i=0;i<n;++i) ++cnt[(k[i]>>sh)&0xFFu];
        size_t sum=0;
        for (int i=0;i<256;++i){ size_t c=cnt[i]; cnt[i]=sum; sum+=c; }
        for (size_t i=0;i<n;++i){
            unsigned b=(unsigned)((k[i]>>sh)&0xFFu);
            size_t p = cnt[b]++;
            tk[p]=k[i]; tv[p]=v[i];
        }
        uint64_t *swk=k; k=tk; tk=swk;
        uint32_t *swv=v; v=tv; tv=swv;
    }
    free(tk); free(tv);
}

#endif 
