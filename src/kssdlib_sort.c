#include "kssdlib_sort.h"
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include <err.h>
#include <errno.h>

#define INSERTION_THRESHOLD 32
#define PARALLEL_THRESHOLD (1 << 24)  // 16M elements per task

/* related methods for already sorted array  */
//1. array deduplicate 
size_t dedup_sorted_uint64(uint64_t *arr, size_t n) {
    if (n <= 1) return n;  // Handle empty/single-element cases

    size_t j = 0;  // Position for next unique element
    
    for (size_t i = 1; i < n; i++) {
        if (arr[i] != arr[j]) {
            j++;
            arr[j] = arr[i];  // Shift unique element forward
        }
    }
    
    return j + 1;  // New array length
}

/**
 * 2. Deduplicates sorted uint64 array and returns occurrence counts
 * 
 * @param arr     Sorted array (modified in-place)
 * @param n       Number of elements in input array
 * @param counts  Output parameter for occurrence counts array
 * @return        New length of deduplicated array
 */
size_t dedup_with_counts(uint64_t *arr, size_t n, uint32_t **counts) {
    if (n == 0) {
        *counts = NULL;
        return 0;
    }

    // Allocate maximum possible size for counts
    uint32_t *cnt = malloc(n * sizeof(uint32_t));
    if (!cnt) {
        // Fallback: deduplicate without counts
        size_t j = 0;
        for (size_t i = 1; i < n; i++) {
            if (arr[i] != arr[j]) arr[++j] = arr[i];
        }
        *counts = NULL;
        return j + 1;
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

//* cumstom array sorting Methods *//
/*0. qsort comparator */
int qsort_comparator_uint32 (const void *a, const void *b){
    const uint32_t val_a = *(const uint32_t *)a;
    const uint32_t val_b = *(const uint32_t *)b;
    return (val_a > val_b) - (val_a < val_b);
}
int qsort_comparator_uint64 (const void *a, const void *b){
	const uint64_t val_a = *(const uint64_t *)a;
	const uint64_t val_b = *(const uint64_t *)b;   
    return (val_a > val_b) - (val_a < val_b);
}

/* Method1: paralle cumstom sort for uint96_t */
// Comparison function for uint96_t
static inline int uint96_less(const uint96_t *a, const uint96_t *b) {
    if (a->part[0] != b->part[0]) return a->part[0] < b->part[0];
    if (a->part[1] != b->part[1]) return a->part[1] < b->part[1];
    return a->part[2] < b->part[2];
}

// Swap function for uint96_t
static inline void swap(uint96_t *a, uint96_t *b) {
    uint32_t tmp[3];
    tmp[0] = a->part[0]; tmp[1] = a->part[1]; tmp[2] = a->part[2];
    a->part[0] = b->part[0]; a->part[1] = b->part[1]; a->part[2] = b->part[2];
    b->part[0] = tmp[0]; b->part[1] = tmp[1]; b->part[2] = tmp[2];
}

// Insertion sort for small arrays
static void insertion_sort(uint96_t *arr, size_t n) {
    for (size_t i = 1; i < n; i++) {
        uint96_t key = arr[i];
        size_t j = i;
        while (j > 0 && uint96_less(&key, &arr[j-1])) {
            arr[j] = arr[j-1];
            j--;
        }
        arr[j] = key;
    }
}

// Partition function for quicksort
static size_t partition(uint96_t *arr, size_t n) {
    uint96_t pivot = arr[n/2];  // Mid-point pivot
    size_t i = 0, j = n - 1;
    while (1) {
        while (uint96_less(&arr[i], &pivot)) i++;
        while (uint96_less(&pivot, &arr[j])) j--;
        if (i >= j) return j;
        swap(&arr[i], &arr[j]);
        i++; j--;
    }
}

// Custom recursive sort
static void custom_sort(uint96_t *arr, size_t n) {
    if (n <= INSERTION_THRESHOLD) {
        insertion_sort(arr, n);
        return;
    }

    // Median-of-three pivot selection
    uint96_t *mid = arr + n/2, *left = arr, *right = arr + n - 1;
    if (uint96_less(mid, left)) swap(left, mid);
    if (uint96_less(right, left)) swap(left, right);
    if (uint96_less(mid, right)) swap(mid, right);
    uint96_t pivot = *right;

    // Hoare partition scheme
    size_t i = 0, j = n - 1;
    while (1) {
        do i++; while (uint96_less(&arr[i], &pivot));
        do j--; while (uint96_less(&pivot, &arr[j]));
        if (i >= j) break;
        swap(&arr[i], &arr[j]);
    }
    swap(&arr[i], &arr[n-1]);

    // Recurse on partitions
    custom_sort(arr, i);
    custom_sort(arr + i + 1, n - i - 1);
}

// Parallel quicksort
void parallel_quicksort(uint96_t *arr, size_t n) {
    if (n <= PARALLEL_THRESHOLD) {
        custom_sort(arr, n);
        return;
    }

    size_t p = partition(arr, n);
    #pragma omp task default(none) firstprivate(arr, p)
    { parallel_quicksort(arr, p + 1); }
    #pragma omp task default(none) firstprivate(arr, p, n)
    { parallel_quicksort(arr + p + 1, n - p - 1); }
    #pragma omp taskwait
}

// Main entry point
void sort_uint96_array(uint96_t *arr, size_t n) {
    #pragma omp parallel num_threads(32)
    #pragma omp single
    parallel_quicksort(arr, n);
}
/*==    Method 1 END ===*/

/*Method 2 for ctxgidobj*/ //typedef struct {uint64_t ctxgid; uint32_t obj;} ctxgidobj_t;

// Swap two ctxgidobj_t elements (optimized for 12-byte struct)
static inline void ctxgidobj_swap(ctxgidobj_t *a, ctxgidobj_t *b) {
    ctxgidobj_t temp = *a;
    *a = *b;
    *b = temp;
}

// Insertion sort for small subarrays
static void ctxgidobj_insertion_sort(ctxgidobj_t *arr, size_t n) {
    for (size_t i = 1; i < n; i++) {
        ctxgidobj_t key = arr[i];
        size_t j = i;
        while (j > 0 && (
            arr[j-1].ctxgid > key.ctxgid ||
            (arr[j-1].ctxgid == key.ctxgid && arr[j-1].obj > key.obj)
        )) {
            arr[j] = arr[j-1];
            j--;
        }
        arr[j] = key;
    }
}

// Unified comparison logic
static inline int ctxgidobj_less(const ctxgidobj_t *a, const ctxgidobj_t *b) {
    return (a->ctxgid < b->ctxgid) || 
          ((a->ctxgid == b->ctxgid) && (a->obj < b->obj));
}

// Partition function for both serial and parallel
static size_t ctxgidobj_partition(ctxgidobj_t *arr, size_t n) {
    size_t i = -1, j = n;
    ctxgidobj_t pivot = arr[n/2];
    
    while (1) {
        do i++; while (ctxgidobj_less(&arr[i], &pivot));
        do j--; while (ctxgidobj_less(&pivot, &arr[j]));
        if (i >= j) return j;
        ctxgidobj_swap(&arr[i], &arr[j]);
    }
}

// Optimized quicksort for ctxgidobj_t
void ctxgidobj_custom_sort(ctxgidobj_t *arr, size_t n) {
    if (n <= INSERTION_THRESHOLD) {
        ctxgidobj_insertion_sort(arr, n);
        return;
    }

    // Correct median-of-three selection
    ctxgidobj_t *left = arr, *mid = arr + n/2, *right = arr + n - 1;
    
    // Order left <= mid <= right
    if (ctxgidobj_less(mid, left)) ctxgidobj_swap(left, mid);
    if (ctxgidobj_less(right, left)) ctxgidobj_swap(left, right);
    if (ctxgidobj_less(right, mid)) ctxgidobj_swap(mid, right);
    
    // Place pivot at position right-1
    ctxgidobj_swap(mid, right - 1);
    ctxgidobj_t pivot = *(right - 1);

    // Hoare partition with safe indices
    size_t i = 0, j = n - 2;
    while (1) {
        do i++; while (ctxgidobj_less(&arr[i], &pivot));
        do j--; while (ctxgidobj_less(&pivot, &arr[j]));
        if (i >= j) break;
        ctxgidobj_swap(&arr[i], &arr[j]);
    }
    
    // Restore pivot to final position
    ctxgidobj_swap(&arr[i], right - 1);

    // Recurse on partitions
    ctxgidobj_custom_sort(arr, i);
    ctxgidobj_custom_sort(arr + i + 1, n - i - 1);
}

// Parallel quicksort with consistent partitioning
void ctxgidobj_parallel_quicksort(ctxgidobj_t *arr, size_t n) {
    if (n <= PARALLEL_THRESHOLD) {
        ctxgidobj_custom_sort(arr, n);
        return;
    }

    size_t p = ctxgidobj_partition(arr, n);
    #pragma omp task default(none) firstprivate(arr, p, n)
    { ctxgidobj_parallel_quicksort(arr, p + 1); }
    #pragma omp task default(none) firstprivate(arr, p, n)
    { ctxgidobj_parallel_quicksort(arr + p + 1, n - p - 1); }
    #pragma omp taskwait
}

// Entry point with thread binding
void ctxgidobj_sort_array(ctxgidobj_t *arr, size_t n) {
    #pragma omp parallel num_threads(32)
    #pragma omp single
    ctxgidobj_parallel_quicksort(arr, n);
}

/*Methods 3 lib for count #overlp intergers between uint64_t a[N]  and b[M]*/
#include <stdint.h>
#include <stdbool.h>

// Two-pointer method for similar-sized arrays
size_t count_overlaps_two_pointers(const uint64_t *a, size_t n, 
                                   const uint64_t *b, size_t m) {
    size_t count = 0, i = 0, j = 0;
    while (i < n && j < m) {
        if (a[i] == b[j]) {
            count++;
            i++;
            j++;
        } else if (a[i] < b[j]) {
            i++;
        } else {
            j++;
        }
    }
    return count;
}

// Binary search helper
bool binary_search(const uint64_t *arr, size_t size, uint64_t target) {
    size_t left = 0, right = size;
    while (left < right) {
//        size_t mid = left + (right - left) * (target - arr[left]) / (arr[right]-arr[left]) ;
        size_t mid = left + (right - left) / 2;
        if (arr[mid] < target) left = mid + 1;
        else right = mid;
    }
    return (left < size && arr[left] == target);
}

// Binary search method for small vs large arrays
size_t count_overlaps_binary_search(const uint64_t *small, size_t small_size,
                                    const uint64_t *large, size_t large_size) {
    size_t count = 0;
    for (size_t i = 0; i < small_size; i++) {
        if (binary_search(large, large_size, small[i])) count++;
    }
    return count;
}

// Hybrid selector (automatically chooses fastest method)
size_t count_overlaps(const uint64_t *a, size_t n, 
                      const uint64_t *b, size_t m) {
    const size_t threshold = 64; // Tune based on cache line size (e.g., 64-256)
    const size_t min_size = (n < m) ? n : m;
    const size_t max_size = (n > m) ? n : m;
    
    // Use binary search if smaller array is below threshold
    if (min_size < threshold || min_size * 20 < max_size) {
        return (n < m) ? count_overlaps_binary_search(a, n, b, m)
                       : count_overlaps_binary_search(b, m, a, n);
    } else {
        return count_overlaps_two_pointers(a, n, b, m);
    }
}

/*Methods 4 sort khash */
int cmp_kv_pair(const void *a, const void *b) {
    const kv_pair_t *pa = (const kv_pair_t *)a;
    const kv_pair_t *pb = (const kv_pair_t *)b;
    if (pa->key < pb->key) return -1;
    else if (pa->key > pb->key) return 1;
    else return 0;
}

// Input: pointer to a khash table of type kmer_hash
// Output: SortedArrays struct containing pointers to the sorted key and value arrays, and the array length.
SortedKV_Arrays_t sort_khash_u64 (khash_t(sort64) *h) {
    SortedKV_Arrays_t result = {0};
    size_t n = kh_size(h);
    result.len = n;

    // Allocate temporary array to hold key-value pairs.
    kv_pair_t *pairs = malloc(n * sizeof(kv_pair_t));
    if (!pairs)  err(EXIT_FAILURE, "%s():Memory allocation failed.",__func__);       
    // Iterate over the hash table to collect valid entries.
    size_t idx = 0;
    khiter_t k;
    for (k = kh_begin(h); k != kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            pairs[idx].key = kh_key(h, k);
            pairs[idx].value = kh_value(h, k);
            idx++;
        }
    }
    // Sort the key-value pairs by key.
    qsort(pairs, n, sizeof(kv_pair_t), cmp_kv_pair);

    // Allocate arrays for the sorted keys and values.
    result.keys = malloc(n * sizeof(uint64_t));
    result.values = malloc(n * sizeof(uint32_t));
    if (!result.keys || !result.values) err(EXIT_FAILURE, "%s():Memory allocation failed.",__func__);
    // Copy the sorted keys and values into the arrays.
    for (size_t i = 0; i < n; i++) {
        result.keys[i] = pairs[i].key;
        result.values[i] = pairs[i].value;
    }
    free(pairs);
    return result;
}

void filter_n_SortedKV_Arrays(SortedKV_Arrays_t *result, uint32_t n){
	uint32_t kmer_ct  = 0;
  for(uint32_t i=0; i< result->len; i++) {
  	if(result->values[i] >= n) {
    	result->keys[kmer_ct] = result->keys[i];
	    result->values[kmer_ct] = result->values[i];
      kmer_ct++;
    }
  }
	result->len	= kmer_ct;
}

void remove_ctx_with_conflict_obj (SortedKV_Arrays_t *data, uint32_t n_obj_bits) {
    if (data->len == 0) return;
    size_t write_idx = 0,i = 0;

    while (i < data->len) {
        size_t count = 1;
        // Efficient high-bit comparison using XOR+shift
        while (i + count < data->len &&
               ((data->keys[i] ^ data->keys[i + count]) >> n_obj_bits) == 0) {
            count++;
        }

        if (count == 1) {
            if (write_idx != i) {
                data->keys[write_idx]   = data->keys[i];
                data->values[write_idx] = data->values[i];
            }
            write_idx++;
        }

        i += count;  // Skip current group
    }

    data->len = write_idx;
}

/*
void free_SortedKV_Arrays (SortedKV_Arrays_t result){
	free(result.keys);
	if(result.values != NULL) free(result->values);
	
}

*/
