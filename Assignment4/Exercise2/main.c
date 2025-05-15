/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-09
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#define N 1000 // Matrix size

void initialize_matrix(double matrix[N*N]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			matrix[i*N+j] = i + j * 0.01;
		}
	}
}

void compute_row_sums(double matrix[N*N], double row_sums[N], int rank, int nprocs) {
	double buff[N*N];
	double totalsum = 0.0;
    int nrows = N / nprocs;
    int extraforlast = N % nprocs;

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
	MPI_Scatterv(&matrix[0], sendcounts, displs, MPI_DOUBLE, &buff[0], N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Figure out how many rows should be processed by this process
	if (rank == nprocs - 1) nrows += extraforlast;
	double lsums[N];
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
	MPI_Gatherv(&lsums[0], nrows, MPI_DOUBLE, &row_sums[0], recvcounts, gatherdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Reduce(&thissum, &totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		free(sendcounts);
		free(recvcounts);
		free(displs);
		free(gatherdispls);
		printf("totalsum = %f\n", totalsum);
	}
}

void write_output(double row_sums[N]) {
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
	double matrix[N*N], row_sums[N];
	if (rank == 0) {
		initialize_matrix(matrix);
	}
	compute_row_sums(matrix, row_sums, rank, nprocs);

	MPI_Finalize();

	if (rank == 0) {
		write_output(row_sums);
		printf("Row sum computation complete.\n");
	}
	return 0;
}

