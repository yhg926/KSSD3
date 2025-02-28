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


#endif 
