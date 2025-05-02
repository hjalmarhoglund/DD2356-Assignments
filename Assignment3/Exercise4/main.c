// the code calculates a DFT of a random complex number input and
// then an IDFT. The IDFT result should be the input vector
// to compile with gcc
// gcc -Wall -O2 -fopenmp -o DFTW DFTW.c
// written by stef

// exercise

#include "stdio.h"  // printf
#include "stdlib.h" // malloc and rand for instance. Rand not thread safe!
#include "time.h"   // time(0) to get random seed
#include "math.h"   // sine and cosine
#include "omp.h"    // openmp library like timing
#include "assert.h"

#define NUM_WARMUP 3
#define NUM_REAL 10

// two pi
#define PI2 6.28318530718
// this for the rounding error, increasing N rounding error increases
// 0.01 precision good for N > 8000
#define R_ERROR 0.01

// main routine to calculate DFT
int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N);
int DFT_omp(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N);
int DFT_omp_schedule(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N);
int DFT_omp_reduction(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N);
int DFT_swap(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N);
int DFT_manual(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N);
void perform_test(int version, double *xr, double *xi, double *Xr_o, double *Xi_o, double *xr_check, double *xi_check, int N);
// set the input array with random number
int fillInput(double *xr, double *xi, int N);
// set to zero the input vector
int setOutputZero(double *Xr_o, double *Xi_o, int N);
// check if x = IDFT(DFT(x))
int checkResults(double *xr, double *xi, double *xr_check, double *xi_check,
        double *Xr_o, double *Xi_r, int N);
// print the results of the DFT
int printResults(double *xr, double *xi, int N);

int main(int argc, char *argv[]) {
    // size of input array
    int N = 10000; // 8,000 is a good number for testing
    printf("DFTW calculation with N = %d \n", N);

    // Allocate array for input vector
    double *xr = (double *)malloc(N * sizeof(double));
    double *xi = (double *)malloc(N * sizeof(double));
    fillInput(xr, xi, N);

    // for checking purposes
    double *xr_check = (double *)malloc(N * sizeof(double));
    double *xi_check = (double *)malloc(N * sizeof(double));

    // Allocate array for output vector
    double *Xr_o = (double *)malloc(N * sizeof(double));
    double *Xi_o = (double *)malloc(N * sizeof(double));

    for (int ver = 0; ver <= 5; ver++) {
        for (int warm = 0; warm < NUM_WARMUP; warm++) {
            perform_test(ver, xr, xi, Xr_o, Xi_o, xr_check, xi_check, N);
        }
        for (int real = 0; real < NUM_REAL; real++) {
            perform_test(ver, xr, xi, Xr_o, Xi_o, xr_check, xi_check, N);
        }
    }

    // print the results of the DFT
#ifdef DEBUG
    printResults(Xr_o, Xi_o, N);
#endif

    // take out the garbage
    free(xr);
    free(xi);
    free(Xi_o);
    free(Xr_o);
    free(xr_check);
    free(xi_check);

    return 0;
}

void perform_test(int version, double *xr, double *xi, double *Xr_o, double *Xi_o, double *xr_check, double *xi_check, int N) {
    assert(0 <= version && version <= 5);
    // Reset
    setOutputZero(xr_check, xi_check, N);
    setOutputZero(Xr_o, Xi_o, N);
    double start_time = omp_get_wtime();

    switch (version) {
        case 0:
            DFT(1, xr, xi, Xr_o, Xi_o, N);
            DFT(-1, Xr_o, Xi_o, xr_check, xi_check, N);
            break;
        case 1:
            DFT_omp(1, xr, xi, Xr_o, Xi_o, N);
            DFT_omp(-1, Xr_o, Xi_o, xr_check, xi_check, N);
            break;
        case 2:
            DFT_omp_schedule(1, xr, xi, Xr_o, Xi_o, N);
            DFT_omp_schedule(-1, Xr_o, Xi_o, xr_check, xi_check, N);
            break;
        case 3:
            DFT_omp_reduction(1, xr, xi, Xr_o, Xi_o, N);
            DFT_omp_reduction(-1, Xr_o, Xi_o, xr_check, xi_check, N);
            break;
        case 4:
            DFT_swap(1, xr, xi, Xr_o, Xi_o, N);
            DFT_swap(-1, Xr_o, Xi_o, xr_check, xi_check, N);
            break;
        case 5:
            DFT_manual(1, xr, xi, Xr_o, Xi_o, N);
            DFT_manual(-1, Xr_o, Xi_o, xr_check, xi_check, N);
            break;
        default:     
            printf("BAD VERSION!\n");
    }
    // stop timer
    double run_time = omp_get_wtime() - start_time;
    char *version_names[] = {"Serial", "OMP", "OMP_SCHEDULE", "OMP_REDUCTION", "SWAP", "MANUAL"};
    printf("DFTW %s computation in %f seconds\n", version_names[version], run_time);
    // check the results: easy to make correctness errors with openMP
    checkResults(xr, xi, xr_check, xi_check, Xr_o, Xi_o, N);
}

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int DFT_manual(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N) {
    int num_threads, thread_num;
    double Xr_o_local, Xi_o_local;
#pragma omp parallel shared(num_threads) private(thread_num,Xr_o_local,Xi_o_local)
    {
    thread_num = omp_get_thread_num();
#pragma single
    {
    num_threads = omp_get_num_threads();
    }
    for (int k = thread_num; k < N; k += num_threads) {
        Xr_o_local = 0.0;
        Xi_o_local = 0.0;
        for (int n = 0; n < N; n++) {
            // Real part of X[k]
            Xr_o_local +=
                xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            // Imaginary part of X[k]
            Xi_o_local +=
                -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
        Xr_o[k] = Xr_o_local;
        Xi_o[k] = Xi_o_local;
    }

    // normalize if you are doing IDFT
    if (idft == -1) {
        for (int n = thread_num; n < N; n += num_threads) {
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }
    }
    return 1;
}

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int DFT_swap(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N) {
#pragma omp parallel
    {
    for (int n = 0; n < N; n++) {
#pragma omp for
        for (int k = 0; k < N; k++) {
            // Real part of X[k]
            Xr_o[k] +=
                xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            // Imaginary part of X[k]
            Xi_o[k] +=
                -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
    }

    // normalize if you are doing IDFT
    if (idft == -1) {
#pragma omp for
        for (int n = 0; n < N; n++) {
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }
    }
    return 1;
}

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int DFT_omp_reduction(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N) {
    for (int k = 0; k < N; k++) {
        double Xro_k = 0.0;
        double Xio_k = 0.0;
#pragma omp parallel for reduction(+:Xro_k,Xio_k)
        for (int n = 0; n < N; n++) {
            // Real part of X[k]
            Xro_k +=
                xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            // Imaginary part of X[k]
            Xio_k +=
                -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
        Xr_o[k] = Xro_k;
        Xi_o[k] = Xio_k;
    }

    // normalize if you are doing IDFT
    if (idft == -1) {
#pragma parallel omp for
        for (int n = 0; n < N; n++) {
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }
    return 1;
}

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int DFT_omp_schedule(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N) {
#pragma omp parallel for schedule(static,128)
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            // Real part of X[k]
            Xr_o[k] +=
                xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            // Imaginary part of X[k]
            Xi_o[k] +=
                -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
    }

    // normalize if you are doing IDFT
    if (idft == -1) {
#pragma omp parallel for schedule(static,128)
        for (int n = 0; n < N; n++) {
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }
    return 1;
}

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int DFT_omp(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N) {
#pragma omp parallel for
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            // Real part of X[k]
            Xr_o[k] +=
                xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            // Imaginary part of X[k]
            Xi_o[k] +=
                -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
    }

    // normalize if you are doing IDFT
    if (idft == -1) {
#pragma omp parallel for
        for (int n = 0; n < N; n++) {
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }
    return 1;
}

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N) {
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            // Real part of X[k]
            Xr_o[k] +=
                xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            // Imaginary part of X[k]
            Xi_o[k] +=
                -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
    }

    // normalize if you are doing IDFT
    if (idft == -1) {
        for (int n = 0; n < N; n++) {
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }
    return 1;
}

// set the initial signal
// be careful with this
// rand() is NOT thread safe in case
int fillInput(double *xr, double *xi, int N) {
    srand(time(0));
    for (int n = 0; n < 100000; n++) // get some random number first
        rand();
    for (int n = 0; n < N; n++) {
        // Generate random discrete-time signal x in range (-1,+1)
        // xr[n] = ((double)(2.0 * rand()) / RAND_MAX) - 1.0;
        // xi[n] = ((double)(2.0 * rand()) / RAND_MAX) - 1.0;
        // constant real signal
        xr[n] = 1.0;
        xi[n] = 0.0;
    }
    return 1;
}

// set to zero the output vector
int setOutputZero(double *Xr_o, double *Xi_o, int N) {
    for (int n = 0; n < N; n++) {
        Xr_o[n] = 0.0;
        Xi_o[n] = 0.0;
    }
    return 1;
}

// check if x = IDFT(DFT(x))
int checkResults(double *xr, double *xi, double *xr_check, double *xi_check,
        double *Xr_o, double *Xi_r, int N) {
    // x[0] and x[1] have typical rounding error problem
    // interesting there might be a theorem on this
    for (int n = 0; n < N; n++) {
        if (fabs(xr[n] - xr_check[n]) > R_ERROR)
            printf("ERROR - x[%d] = %f, inv(X)[%d]=%f \n", n, xr[n], n, xr_check[n]);
        if (fabs(xi[n] - xi_check[n]) > R_ERROR)
            printf("ERROR - x[%d] = %f, inv(X)[%d]=%f \n", n, xi[n], n, xi_check[n]);
    }
    printf("Xre[0] = %f \n", Xr_o[0]);
    return 1;
}

// print the results of the DFT
int printResults(double *xr, double *xi, int N) {
    for (int n = 0; n < N; n++)
        printf("Xre[%d] = %f, Xim[%d] = %f \n", n, xr[n], n, xi[n]);
    return 1;
}
