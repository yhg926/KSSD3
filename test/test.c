#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#define ARRAY_SIZE ( 1024L * 1024 )  // 1 GB of integers
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
    
    // Allocate memory for the array (10 GB of integers)
    short *array = (short *)malloc(ARRAY_SIZE * sizeof(short));
    if (array == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }
		    // Generate random integers
	  for (long i = 0; i < ARRAY_SIZE; i++) {
        array[i] = rand() % LEN ;
    }

    // Start the timer
    clock_gettime(CLOCK_MONOTONIC, &start);
	#pragma omp parallel for reduction(min:min)
		for (long i = 0; i < ARRAY_SIZE; i++) {
			short	temp =  space[array[i]] ;							
				if (temp < min ){
									min = temp ;
				}
    }

    // Stop the timer
    clock_gettime(CLOCK_MONOTONIC, &end);

    // Print the time taken
    double elapsed_ns = time_diff(&start, &end);
    printf("Time of comparing %ld array element to %d: %.0f nanoseconds\n", ARRAY_SIZE, min, elapsed_ns);

    // Free the allocated memory
    free(array);

    return 0;
}

