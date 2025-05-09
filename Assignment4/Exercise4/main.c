/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-09
 */

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#define N 1000  // Matrix size

void initialize(double *matrix, double *vector) {
    for (int i = 0; i < N * N; i++) {
        matrix[i] = (double)(i % 100) / 10.0;
    }
    for (int i = 0; i < N; i++) {
        vector[i] = (double)(i % 50) / 5.0;
    }
}

int main() {
    double *matrix = (double*) malloc(N * N * sizeof(double));
    double *vector = (double*) malloc(N * sizeof(double));
    double *result = (double*) malloc(N * sizeof(double));

    initialize(matrix, vector);

    // Perform matrix-vector multiplication using BLAS
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, matrix, N, vector, 1, 0.0, result, 1);

    // Write output to file
    FILE *f = fopen("blas_output.txt", "w");
    for (int i = 0; i < N; i++) {
        fprintf(f, "%f\n", result[i]);
    }
    fclose(f);

    free(matrix);
    free(vector);
    free(result);
    printf("BLAS matrix-vector multiplication complete.\n");
    return 0;
}
