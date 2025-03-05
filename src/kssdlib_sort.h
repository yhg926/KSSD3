#ifndef KSSDLIB_SORT_H
#define KSSDLIB_SORT_H
#include "command_ani.h"

/*0. qsort comparator */
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

/*3. methods for already sortted array */
size_t dedup_sorted_uint64(uint64_t *arr, size_t n);
size_t dedup_with_counts(uint64_t *arr, size_t n, uint32_t **counts) ;

/*4. khash sort*/
#include "../klib/khash.h"
KHASH_MAP_INIT_INT64(sort64, uint32_t)
typedef struct { uint64_t key;  uint32_t value;} kv_pair_t;
typedef struct { int64_t *keys;int *values;size_t len;} SortedKV_Arrays_t ;
int cmp_kv_pair(const void *a, const void *b);
SortedKV_Arrays_t sort_khash_u64 (khash_t(sort64) *h);
void filter_n_SortedKV_Arrays(SortedKV_Arrays_t *result, uint32_t n);
//void free_SortedKV_Arrays (SortedKV_Arrays_t *result);
#endif 
