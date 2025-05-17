/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-09
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#define NUM_RUNS 10
#define NUM_WARMUP 3

void initialize_matrix(double *matrix, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			matrix[i*N+j] = i + j * 0.01;
		}
	}
}

double compute_row_sums(double *matrix, double *row_sums, int N, int rank, int nprocs) {
	double totalsum = 0.0;
	int nrows = N / nprocs;
	int extraforlast = N % nprocs;
	if (rank == nprocs - 1) nrows += extraforlast;

	double *buff = (double*) malloc(nrows * N * sizeof(double));
	double *lsums = (double*) malloc(nrows * sizeof(double));

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
	MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, buff, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Figure out how many rows should be processed by this process
	double *p = buff;
	double thissum = 0.0;
	for (int i = 0; i < nrows; i++) {
		lsums[i] = 0.0;
		for (int j = 0; j < N; j++) {
			lsums[i] += *p;
			thissum += *p;
			p++;
		}
	}
	MPI_Gatherv(lsums, nrows, MPI_DOUBLE, row_sums, recvcounts, gatherdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Reduce(&thissum, &totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	free(buff);
	free(lsums);

	if (rank == 0) {
		free(sendcounts);
		free(recvcounts);
		free(displs);
		free(gatherdispls);
		return totalsum;
	}
	return 0.0;
}

void write_output(double *row_sums, int N) {
	FILE *f = fopen("row_sums_output.txt", "w");
	for (int i = 0; i < N; i++) {
		fprintf(f, "%f\n", row_sums[i]);
	}
	fclose(f);
}

int main(int argc, char **argv) {
	int rank, nprocs, provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	assert(nprocs > 1);
	int N = 500;
	for (int n = 0; n < 5; n++) {
		N *= 2;
		double *matrix;
		double *row_sums;
		if (rank == 0) {
			matrix = (double*) malloc(N * N * sizeof(double));
			row_sums = (double*) malloc(N * N * sizeof(double));
			initialize_matrix(matrix, N);
			printf("N = %d, nprocs = %d\n", N, nprocs);
		}
		for (int i = 0; i < NUM_WARMUP; i++) {
			double tval = compute_row_sums(matrix, row_sums, N, rank, nprocs);
			if (rank == 0 && tval < 0.0) printf("tval less than 0!\n");
		}

		double t_start, t;
		for (int i = 0; i < NUM_RUNS; i++) {
			if (rank == 0) t_start = MPI_Wtime();
			double tval = compute_row_sums(matrix, row_sums, N, rank, nprocs);
			if (rank == 0) {
				t = MPI_Wtime() - t_start;
				printf("Run %d took %.6f seconds\n", i, t);
			}
			if (rank == 0 && tval < 0.0) printf("tval less than 0!\n");
		}

		if (rank == 0) {
			free(matrix);
			free(row_sums);
		}
	}

	MPI_Finalize();

	return 0;
}


