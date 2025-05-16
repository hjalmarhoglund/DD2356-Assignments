/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-16
 */

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <mpi.h>

#define N 1000  // Matrix size

void initialize(double *matrix, double *vector) {
    for (int i = 0; i < N * N; i++) {
        matrix[i] = (double)(i % 100) / 10.0;
    }
    for (int i = 0; i < N; i++) {
        vector[i] = (double)(i % 50) / 5.0;
    }
}

int main(int argc, char **argv) {
    int rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int nrows = N / nprocs;
    int extraforlast = N % nprocs;

    if (rank == nprocs - 1) nrows += extraforlast;

    double *submatrix = (double*) malloc(N * nrows * sizeof(double));
    double *subresult = (double*) malloc(nrows * sizeof(double));
    double *vector = (double*) malloc(N * sizeof(double));
    double *matrix;
    double *result;
    int *sendcounts;
    int *recvcounts;
    int *displs;
    int *gatherdispls;


    if (rank == 0) {
        matrix = (double*) malloc(N * N * sizeof(double));
        result = (double*) malloc(N * sizeof(double));
        initialize(matrix, vector);

        sendcounts = (int *)malloc(nprocs * sizeof(int));
        recvcounts = (int *)malloc(nprocs * sizeof(int));
        displs = (int *)malloc(nprocs * sizeof(int));
        gatherdispls = (int *)malloc(nprocs * sizeof(int));
        for (int i = 0; i < nprocs; i++) {
            sendcounts[i] = nrows * N;
            recvcounts[i] = nrows;
            displs[i] = nrows * i * N;
            gatherdispls[i] = nrows * i;
        }
        sendcounts[nprocs-1] += extraforlast * N;
        recvcounts[nprocs-1] += extraforlast;
    }
    // Broadcast the vector
    MPI_Bcast(vector, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Send the part of the matrix
    MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, submatrix, N * nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Perform matrix-vector multiplication using BLAS
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nrows, N, 1.0, submatrix, N, vector, 1, 0.0, subresult, 1);

    MPI_Gatherv(subresult, nrows, MPI_DOUBLE, result, recvcounts, gatherdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Finalize();
    free(submatrix);
    free(subresult);
    free(vector);

    // Other processes are done now.
    if (rank != 0) return 0;

    // Write output to file
    FILE *f = fopen("blas_output.txt", "w");
    for (int i = 0; i < N; i++) {
        fprintf(f, "%f\n", result[i]);
    }
    fclose(f);

    free(matrix);
    free(result);

    free(sendcounts);
    free(recvcounts);
    free(displs);
    free(gatherdispls);

    printf("BLAS matrix-vector multiplication complete.\n");
    return 0;
}
