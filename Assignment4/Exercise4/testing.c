/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-16
*/

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <mpi.h>

#define NUM_RUNS 10
#define NUM_WARMUP 3

void initialize(double *matrix, double *vector, int N) {
	for (int i = 0; i < N * N; i++) {
		matrix[i] = (double)(i % 100) / 10.0;
	}
	for (int i = 0; i < N; i++) {
		vector[i] = (double)(i % 50) / 5.0;
	}
}

double experiment(int rank, int nprocs, double *matrix, double *vector, double *result, int N) {
	int nrows = N / nprocs;
	int extraforlast = N % nprocs;

	if (rank == nprocs - 1) nrows += extraforlast;

	double *submatrix = (double*) malloc(N * nrows * sizeof(double));
	double *subresult = (double*) malloc(nrows * sizeof(double));

	int *sendcounts;
	int *recvcounts;
	int *displs;
	int *gatherdispls;

	if (rank == 0) {
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

	free(submatrix);
	free(subresult);
	
	// Other processes are done now.
	if (rank != 0) return 0.0;
	double ret = result[0];

	free(sendcounts);
	free(recvcounts);
	free(displs);
	free(gatherdispls);
	
	return ret;
}

int main(int argc, char **argv) {
	int rank, nprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	double *matrix;
	double *vector;
	double *result;

	int N = 500;
	for (int n = 0; n < 5; n++) {
		N *= 2;

		vector = (double*) malloc(N * sizeof(double));
		if (rank == 0) {
			matrix = (double*) malloc(N * N * sizeof(double));
			result = (double*) malloc(N * sizeof(double));
			initialize(matrix, vector, N);
		}

		for (int i = 0; i < NUM_WARMUP; i++) {
			double tval = experiment(rank, nprocs, matrix, vector, result, N);
			if (rank == 0 && tval < 0.0) printf("tval less than 0!\n");
		}

		double t_start, t;
		for (int i = 0; i < NUM_RUNS; i++) {
			if (rank == 0) t_start = MPI_Wtime();
			double tval = experiment(rank, nprocs, matrix, vector, result, N);
			if (rank == 0) {
				t = MPI_Wtime() - t_start;
				printf("Run %d took %.6f seconds\n", i, t);
			}
			if (rank == 0 && tval < 0.0) printf("tval less than 0!\n");
		}

		free(vector);
		if (rank == 0) {
			free(matrix);
			free(result);
		}
	}
	MPI_Finalize();

	return 0;
}
