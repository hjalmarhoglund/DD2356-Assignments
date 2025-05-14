/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-09
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

#define N 1000 // Grid size
#define STEPS 100 // Time steps
#define C 1.0   // Wave speed
#define DT 0.01 // Time step
#define DX 1.0  // Grid spacing

double u[N], u_prev[N];

int max(int x, int y) {
	if (x < y) return y;
	return x;
}

int min(int x, int y) {
	if (x < y) return x;
	return y;
}

void initialize(int istart, int iend) {
	for (int i = istart; i < iend; i++) {
		u[i] = sin(2.0 * M_PI * i / N);
		u_prev[i] = u[i];
	}
}

void compute_step(int istart, int iend) {
	double u_next[N];
	for (int i = max(istart, 1); i < min(iend, N - 1); i++) {
		u_next[i] = 2.0 * u[i] - u_prev[i] + C * C * DT * DT / (DX * DX) * (u[i+1] - 2.0 * u[i] + u[i-1]);
	}
	for (int i = istart; i < iend; i++) {
		u_prev[i] = u[i];
		u[i] = u_next[i];
	}
}

void write_output(int step) {
	char filename[50];
	sprintf(filename, "wave_output_%d.txt", step);
	FILE *f = fopen(filename, "w");
	for (int i = 0; i < N; i++) {
		fprintf(f, "%f\n", u[i]);
	}
	fclose(f);
}

void sync_arr(double* arr, int rank, int nprocs, int istart, int iend, int stride) {
	for (int i = 0; i < nprocs; i++) {
		if (i == rank) continue;
		int ostart = i * stride;
		int oend = (i+1) * stride;
		if (i == nprocs - 1) oend = N;
		MPI_Sendrecv(arr+istart, iend-istart, MPI_DOUBLE, i, 0,
			     arr+ostart, oend-ostart, MPI_DOUBLE, i, 0,
			     MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
	}
}

int main(int argc, char* argv[]) {
	int rank, nprocs, provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	printf("MPI rank %d out of %d processes.\n", rank, nprocs);
	assert(nprocs > 1);
	int stride = N / nprocs;
	int istart = rank * stride;
	int iend = (rank + 1) * stride;
	if (rank == nprocs - 1) iend = N;

	initialize(istart, iend);
	sync_arr(&u[0], rank, nprocs, istart, iend, stride);
	memcpy(&u_prev, &u, N);
	for (int step = 0; step < STEPS; step++) {
		compute_step(istart, iend);
		sync_arr(&u_prev[0], rank, nprocs, istart, iend, stride);
		sync_arr(&u[0], rank, nprocs, istart, iend, stride);
		if (rank == 0 && step % 10 == 0) write_output(step);
	}
	MPI_Finalize();

	if(rank==0)
		printf("Simulation complete.\n");
	return 0;
}

