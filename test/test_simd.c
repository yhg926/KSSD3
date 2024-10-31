#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>  // For AVX2 intrinsics

#define ARRAY_SIZE (1024L * 1024L)  // 1 GB of integers
#define LEN (4096)

// Function to measure the time taken
double time_diff(struct timespec *start, struct timespec *end) {
    return (end->tv_sec - start->tv_sec) * 1E9 + (end->tv_nsec - start->tv_nsec);
}

int space[LEN];

int main() {
    int min = 4096;

    // Fill the space array with values 0 to 4095
    for (int i = 0; i < LEN; i++) {
        space[i] = i;
    }

    struct timespec start, end;

    // Allocate memory for the array (1 GB of integers)
    int *array = (int *)malloc(ARRAY_SIZE * sizeof(int));
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

    // SIMD optimization using AVX2
    __m256i vec_min = _mm256_set1_epi16(min);  // Set all elements of the vector to the initial min value
    for (long i = 0; i < ARRAY_SIZE; i += 16) {
        // Load 16 elements from the array (indices) and gather corresponding values from 'space'
        __m256i indices = _mm256_loadu_si256((__m256i*)&array[i]);
        __m256i values = _mm256_i32gather_epi32((const int*)space, indices, 2);

        // Compare the values with current min
        vec_min = _mm256_min_epi16(vec_min, values);
    }

    // Extract the minimum value from the SIMD register
    short simd_min[16];
    _mm256_storeu_si256((__m256i*)simd_min, vec_min);
    for (int i = 0; i < 16; i++) {
        if (simd_min[i] < min) {
            min = simd_min[i];
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

