#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <emmintrin.h> // For SSE2 intrinsics

#define ARRAY_SIZE (1024L * 1024)  // 1 GB of integers
#define LEN (4096)

// Function to measure the time taken
double time_diff(struct timespec *start, struct timespec *end) {
    return (end->tv_sec - start->tv_sec) * 1E9 + (end->tv_nsec - start->tv_nsec);
}

unsigned short space[LEN];

int main() {

    int min = 4096;
    for (int i = 0; i < LEN; i++) {
       space[i] = i;
    }

    struct timespec start, end;

    // Allocate memory for the array (1 GB of integers)
    short *array = (short *)aligned_alloc(16, ARRAY_SIZE * sizeof(short));
    if (array == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    // Generate random integers
    for (long i = 0; i < ARRAY_SIZE; i++) {
        array[i] = rand() % LEN;
    }

    // Start the timer
    clock_gettime(CLOCK_MONOTONIC, &start);

    __m128i min_vec = _mm_set1_epi16(min); // Initialize SIMD min value to 4096
    for (long i = 0; i < ARRAY_SIZE; i += 8) {
        // Load 8 short elements from the array
        __m128i index_vec = _mm_load_si128((__m128i*)&array[i]);

        // Gather the corresponding values from the `space` array
        __m128i space_vals = _mm_set_epi16(
            space[_mm_extract_epi16(index_vec, 7)],
            space[_mm_extract_epi16(index_vec, 6)],
            space[_mm_extract_epi16(index_vec, 5)],
            space[_mm_extract_epi16(index_vec, 4)],
            space[_mm_extract_epi16(index_vec, 3)],
            space[_mm_extract_epi16(index_vec, 2)],
            space[_mm_extract_epi16(index_vec, 1)],
            space[_mm_extract_epi16(index_vec, 0)]
        );

        // Compare and find the minimum values
        min_vec = _mm_min_epi16(min_vec, space_vals);
    }

    // Reduce min_vec to a single minimum value
    short min_array[8];
    _mm_storeu_si128((__m128i*)min_array, min_vec);

    for (int i = 0; i < 8; i++) {
        if (min_array[i] < min) {
            min = min_array[i];
        }
    }

    // Stop the timer
    clock_gettime(CLOCK_MONOTONIC, &end);

    // Print the time taken
    double elapsed_ns = time_diff(&start, &end);
    printf("Time of comparing %ld array elements to %d: %.0f nanoseconds\n", ARRAY_SIZE, min, elapsed_ns);

    // Free the allocated memory
    free(array);

    return 0;
}

