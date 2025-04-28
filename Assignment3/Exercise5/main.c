/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 3
Date: 2025-04-28
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define NUM_WARMUP 3
#define NUM_RUNS 10

#define N 500  // Grid size
#define ITER 1000  // Number of iterations
#define DT 0.01  // Time step
#define DX 1.0   // Grid spacing

double h[N][N], u[N][N], v[N][N];

double calc_sum() {
    double s = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            s += h[i][j];
        }
    }
    return s;
}

void initialize() {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            h[i][j] = 1.0;
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
}

void compute() {
    for (int iter = 0; iter < ITER; iter++) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                double dudx = (u[i+1][j] - u[i-1][j]) / (2.0 * DX);
                double dvdy = (v[i][j+1] - v[i][j-1]) / (2.0 * DX);
                h[i][j] -= DT * (dudx + dvdy);
            }
        }
    }
}

void compute_parallel() {
#pragma omp parallel
    {
    for (int iter = 0; iter < ITER; iter++) {
#pragma omp for collapse(2)
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                double dudx = (u[i+1][j] - u[i-1][j]) / (2.0 * DX);
                double dvdy = (v[i][j+1] - v[i][j-1]) / (2.0 * DX);
                h[i][j] -= DT * (dudx + dvdy);
            }
        }
    }
    }
}

void write_output() {
    FILE *f = fopen("output.txt", "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(f, "%f ", h[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void measure_serial() {
    // Do warmup
    for (int i = 0; i < NUM_WARMUP; i++) {
        initialize();
        compute();
        double s = calc_sum();
        printf("Warmup run %d got s=%.4f\n", i, s);
    }
    double start_time, end_time;
    for (int i = 0; i < NUM_RUNS; i++) {
        initialize();

        start_time = omp_get_wtime();
        compute();
        end_time = omp_get_wtime();

        double s = calc_sum();
        printf("Run %d took %.6f seconds and got s=%.4f\n", i, end_time - start_time, s);
    }
}

void measure_parallel() {
    // Do warmup
    for (int i = 0; i < NUM_WARMUP; i++) {
        initialize();
        compute_parallel();
        double s = calc_sum();
        printf("Warmup run %d got s=%.4f\n", i, s);
    }
    double start_time, end_time;
    for (int i = 0; i < NUM_RUNS; i++) {
        initialize();

        start_time = omp_get_wtime();
        compute_parallel();
        end_time = omp_get_wtime();

        double s = calc_sum();
        printf("Run %d took %.6f seconds and got s=%.4f\n", i, end_time - start_time, s);
    }
}

int main() {
    printf("=== MEASURE SERIAL   ===\n");
    measure_serial();
    printf("===     END SERIAL   ===\n");

    printf("=== MEASURE PARALLEL ===\n");
    measure_parallel();
    printf("===     END PARALLEL ===\n");

    printf("Computation completed.\n");
    return 0;
}

