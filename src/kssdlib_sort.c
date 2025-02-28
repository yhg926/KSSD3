#include "kssdlib_sort.h"
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>

#define INSERTION_THRESHOLD 32
#define PARALLEL_THRESHOLD (1 << 24)  // 16M elements per task

/*0. qsort comparator */
int qsort_comparator_uint64 (const void *a, const void *b){
	const uint64_t val_a = *(const uint64_t *)a;
	const uint64_t val_b = *(const uint64_t *)b;   
    return (val_a > val_b) - (val_a < val_b);
}

/* Method1: paralle cumstom sort */
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
