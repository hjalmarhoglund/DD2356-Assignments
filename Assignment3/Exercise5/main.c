/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 3
Date: 2025-04-28
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 1400           // Grid size
#define ITER 1000       // Number of iterations
#define NUM_WARMUP 3    // Number of warm‐up runs
#define NUM_RUNS 10     // Number of timed runs

double h[N][N], u[N][N], v[N][N];

double calc_sum() {
    double s = 0.0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            s += h[i][j];
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

void compute_serial() {
    for (int iter = 0; iter < ITER; iter++) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                double dudx = (u[i+1][j] - u[i-1][j]) / 2.0;
                double dvdy = (v[i][j+1] - v[i][j-1]) / 2.0;
                h[i][j] -= 0.01 * (dudx + dvdy);
            }
        }
    }
}

void compute_parallel() {
#pragma omp parallel
    {
        for (int iter = 0; iter < ITER; iter++) {
#pragma omp for collapse(2) schedule(runtime)
            for (int i = 1; i < N - 1; i++) {
                for (int j = 1; j < N - 1; j++) {
                    double dudx = (u[i+1][j] - u[i-1][j]) / 2.0;
                    double dvdy = (v[i][j+1] - v[i][j-1]) / 2.0;
                    h[i][j] -= 0.01 * (dudx + dvdy);
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

int main() {
    // disable OpenMP dynamic adjustment for consistency/reproducibility
    omp_set_dynamic(0);

    // serial run
    double serial_times[NUM_RUNS];
    for (int i = 0; i < NUM_WARMUP; i++) {
        initialize();
        compute_serial();
        printf("[warmup serial %d] sum=%.4f\n", i, calc_sum());
    }
    for (int run = 0; run < NUM_RUNS; run++) {
        initialize();
        double t0 = omp_get_wtime();
        compute_serial();
        double t1 = omp_get_wtime();
        serial_times[run] = t1 - t0;
        printf("[serial run %d] time=%.6f  check_sum=%.4f\n",
               run, serial_times[run], calc_sum());
    }
    // best (minimal) serial time
    double t_serial = serial_times[0];
    for (int i = 1; i < NUM_RUNS; i++)
        if (serial_times[i] < t_serial) t_serial = serial_times[i];
    printf("-> Best serial time = %.6f s\n\n", t_serial);

    // schedules and thread counts
    struct {
        omp_sched_t policy;
        const char *name;
    } schedules[] = {
        { omp_sched_static, "static" },
        { omp_sched_dynamic, "dynamic" },
        { omp_sched_guided, "guided" }
    };
    int threads[] = { 1, 2, 4, 8 };
    int nsched = sizeof(schedules)/sizeof(schedules[0]);
    int nth = sizeof(threads)/sizeof(threads[0]);

    for (int s = 0; s < nsched; s++) {
        // set OpenMP schedule
        omp_set_schedule(schedules[s].policy, 1);
        printf("=== Schedule: %s ===\n", schedules[s].name);

        for (int t = 0; t < nth; t++) {
            int nthreads = threads[t];
            omp_set_num_threads(nthreads);

            // warmup
            for (int w = 0; w < NUM_WARMUP; w++) {
                initialize();
                compute_parallel();
                printf("[warmup nthreads=%d] check_sum=%.4f\n",
                       nthreads, calc_sum());
            }

            // timed runs
            double par_times[NUM_RUNS];
            for (int run = 0; run < NUM_RUNS; run++) {
                initialize();
                double t0 = omp_get_wtime();
                compute_parallel();
                double t1 = omp_get_wtime();
                par_times[run] = t1 - t0;
                printf("[parallel run %d nthreads=%d] time=%.6f check_sum=%.4f\n",
                       run, nthreads, par_times[run], calc_sum());
            }
            
	    // best parallel time
            double t_par = par_times[0];
            for (int i = 1; i < NUM_RUNS; i++)
                if (par_times[i] < t_par) t_par = par_times[i];
            double speedup = t_serial / t_par;
            printf("-> Best parallel time (nthreads=%d) = %.6f s, speedup = %.2fx\n\n", nthreads, t_par, speedup);
        }
    }

    // final guided run and output
    omp_set_schedule(omp_sched_guided, 1);
    omp_set_num_threads(threads[nth-1]);
    
    initialize();
    compute_parallel();
    
    write_output();
    printf("Output written to output.txt\n");

    return 0;
}

