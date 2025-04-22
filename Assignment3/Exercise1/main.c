/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 3
Date: 2025-04-22
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main() {
    int num_threads;
#pragma omp parallel
    {
    int thread_id = omp_get_thread_num();
    printf("Hello World from Thread %d!\n", thread_id);
    if (thread_id == 0) {
        num_threads = omp_get_num_threads();
    }
    }
    printf("We have a total of %d threads.\n", num_threads);
    return 0;
}

