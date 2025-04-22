/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 3
Date: 2025-04-22
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define N 10000000 // 10^7
#define NUM_WARMUP 3
#define NUM_RUNS 10

double serial_sum(double *x, size_t size)
{
    double sum_val = 0.0;

    for (size_t i = 0; i < size; i++) {
        sum_val += x[i];
    }

    return sum_val;
}

double omp_sum(double *x, size_t size)
{
    double sum_val = 0.0;
#pragma omp parallel for
    for (size_t i = 0; i < size; i++) {
        sum_val += x[i];
    }

    return sum_val;
}

double omp_critical_sum(double *x, size_t size)
{
    double sum_val = 0.0;
#pragma omp parallel for
    for (size_t i = 0; i < size; i++) {
#pragma omp critical
        {
        sum_val += x[i];
        }
    }

    return sum_val;
}

double omp_local_sum(double *x, size_t size)
{
    double *local_sum;
    size_t num_threads;
#pragma omp parallel shared(local_sum,num_threads)
    {
#pragma omp single
    {
    num_threads = omp_get_num_threads();
    local_sum = (double *) malloc(num_threads * sizeof(double));
    for (size_t i = 0; i < num_threads; i++) {
        local_sum[i] = 0.0;
    }
    }

#pragma omp for
    for (size_t i = 0; i < size; i++) {
        local_sum[omp_get_thread_num()] += x[i];
    }
    }
    double sum_val = 0.0;
    for (size_t i = 0; i < num_threads; i++) {
        sum_val += local_sum[i];
    }
    free(local_sum);
    return sum_val;
}


// TODO: Write omp_local_sum without false sharing. Maybe reduction or shedule(static,chunk) ???

double omp_reduction_sum(double *x, size_t size) {
    double sum_val = 0.0;

#pragma omp parallel for default(shared) reduction(+:sum_val)
    for (size_t i = 0; i < size; i++) {
        sum_val += x[i];
    }
    
    return sum_val;
}

void generate_random(double *input, size_t size)
{
    srand(time(0));
    for (size_t i = 0; i < size; i++) {
        input[i] = rand() / (double)(RAND_MAX);
    }
}

void time_serial_sum(double *arr, size_t size) {
    // Do warmup runs
    for (int i = 0; i < NUM_WARMUP; i++) {
        double s = serial_sum(arr, size);
        printf("Warmup run %d got s=%.4f\n", i, s);
    }
    double start_time, end_time;
    for (int i = 0; i < NUM_RUNS; i++) {
        start_time = omp_get_wtime();
        double s = serial_sum(arr, size);
        end_time = omp_get_wtime();

        printf("Run %d took %.6f seconds and got s=%.4f\n", i, end_time - start_time, s);
    }
}

void time_omp_sum(double *arr, size_t size) {
    // Do warmup runs
    for (int i = 0; i < NUM_WARMUP; i++) {
        double s = omp_sum(arr, size);
        printf("Warmup run %d got s=%.4f\n", i, s);
    }
    double start_time, end_time;
    for (int i = 0; i < NUM_RUNS; i++) {
        start_time = omp_get_wtime();
        double s = omp_sum(arr, size);
        end_time = omp_get_wtime();

        printf("Run %d took %.6f seconds and got s=%.4f\n", i, end_time - start_time, s);
    }
}

void time_omp_critical_sum(double *arr, size_t size) {
    // Do warmup runs
    for (int i = 0; i < NUM_WARMUP; i++) {
        double s = omp_critical_sum(arr, size);
        printf("Warmup run %d got s=%.4f\n", i, s);
    }
    double start_time, end_time;
    for (int i = 0; i < NUM_RUNS; i++) {
        start_time = omp_get_wtime();
        double s = omp_critical_sum(arr, size);
        end_time = omp_get_wtime();

        printf("Run %d took %.6f seconds and got s=%.4f\n", i, end_time - start_time, s);
    }
}

void time_omp_local_sum(double *arr, size_t size) {
    // Do warmup runs
    for (int i = 0; i < NUM_WARMUP; i++) {
        double s = omp_local_sum(arr, size);
        printf("Warmup run %d got s=%.4f\n", i, s);
    }
    double start_time, end_time;
    for (int i = 0; i < NUM_RUNS; i++) {
        start_time = omp_get_wtime();
        double s = omp_local_sum(arr, size);
        end_time = omp_get_wtime();

        printf("Run %d took %.6f seconds and got s=%.4f\n", i, end_time - start_time, s);
    }
}

void time_omp_reduction_sum(double *arr, size_t size) {
    // Do warmup runs
    for (int i = 0; i < NUM_WARMUP; i++) {
        double s = omp_reduction_sum(arr, size);
        printf("Warmup run %d got s=%.4f\n", i, s);
    }
    double start_time, end_time;
    for (int i = 0; i < NUM_RUNS; i++) {
        start_time = omp_get_wtime();
        double s = omp_reduction_sum(arr, size);
        end_time = omp_get_wtime();

        printf("Run %d took %.6f seconds and got s=%.4f\n", i, end_time - start_time, s);
    }
}

int main() {
    double *arr = (double*) malloc(N * sizeof(double));

    generate_random(arr, N);

    printf("===        START SERIAL SUM ===\n");
    time_serial_sum(arr, N);
    printf("===          END SERIAL SUM ===\n");

    printf("===           START OMP SUM ===\n");
    time_omp_sum(arr, N);
    printf("===             END OMP SUM ===\n");

    printf("===  START OMP CRITICAL SUM ===\n");
    //time_omp_critical_sum(arr, N);
    printf("===    END OMP CRITICAL SUM ===\n");

    printf("===     START OMP LOCAL SUM ===\n");
    time_omp_local_sum(arr, N);
    printf("===       END OMP LOCAL SUM ===\n");

    printf("=== START OMP REDUCTION SUM ===\n");
    time_omp_reduction_sum(arr, N);
    printf("===   END OMP REDUCTION SUM ===\n");

    free(arr);
    return 0;
}

