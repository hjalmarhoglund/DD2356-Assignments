/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-09
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#define N 20  // Grid size
#define STEPS 100  // Simulation steps

int **grid;
int **new_grid;
int *neig_u, *neig_d, *neig_l, *neig_r;
int *tc1, *tc2, *tc3, *tc4;

int **malloc_grid(int nrows, int ncols) {
    int **ret = (int **) malloc(nrows * sizeof(int *));
    for (int i = 0; i < nrows; i++) {
        ret[i] = (int *) malloc(ncols * sizeof(int));
    }
    return ret;
}

void malloc_neigs(int sq, int segi, int segj, int nrows, int ncols) {
    // Here we could skip the allocation if we are a edge segment
    neig_u = (int *) malloc(ncols * sizeof(int));
    neig_d = (int *) malloc(ncols * sizeof(int));
    neig_l = (int *) malloc(nrows * sizeof(int));
    neig_r = (int *) malloc(nrows * sizeof(int));
}

void malloc_tmp_cols(int nrows) {
    tc1 = (int *) malloc(nrows * sizeof(int));
    tc2 = (int *) malloc(nrows * sizeof(int));
    tc3 = (int *) malloc(nrows * sizeof(int));
    tc4 = (int *) malloc(nrows * sizeof(int));
}

void initialize(int rank, int sq, int segi, int segj, int nrows, int ncols) {
    grid = malloc_grid(nrows, ncols);
    new_grid = malloc_grid(nrows, ncols);
    malloc_neigs(sq, segi, segj, nrows, ncols);
    malloc_tmp_cols(nrows);
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            grid[i][j] = rand() % 2;  // Random initial state
        }
    }
}

void destroy(int nrows) {
    for (int i = 0; i < nrows; i++) {
        free(grid[i]);
        free(new_grid[i]);
    }
    free(grid);
    free(new_grid);

    free(neig_u);
    free(neig_d);
    free(neig_l);
    free(neig_r);

    free(tc1);
    free(tc2);
    free(tc3);
    free(tc4);
}

int count_neighbors(int x, int y, int nrows, int ncols) {
    int sum = 0;
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0) continue;
            int nx = (x + i + nrows) % nrows;
            int ny = (y + j + ncols) % ncols;
            sum += grid[nx][ny];
        }
    }
    return sum;
}

int ijtorank(int i, int j, int sq) {
    if (i < 0) i += sq;
    if (j < 0) j += sq;
    i = i % sq;
    j = j % sq;
    return i * sq + j;
}

void send_to_neig(int step, int rank, int sq, int segi, int segj, int nrows, int ncols, MPI_Request *req) {
    // Send to above
    int above_rank = ijtorank(segi - 1, segj, sq);
    MPI_Isend(grid[0], ncols, MPI_DOUBLE, above_rank, step, MPI_COMM_WORLD, req[0]);
    // Send to below
    int below_rank = ijtorank(segi + 1, segj, sq);
    MPI_Isend(grid[nrows-1], ncols, MPI_DOUBLE, below_rank, step, MPI_COMM_WORLD, req[1]);
    // To send left and right, we move the columns to tc1 and tc2
    for (int i = 0; i < nrows; i++) {
        tc1[i] = grid[i][0];
        tc2[i] = grid[i][ncols-1];
    }
    // Send left
    int left_rank = ijtorank(segi, segj - 1, sq);
    MPI_Isend(tc1, nrows, MPI_DOUBLE, left_rank, step, MPI_COMM_WORLD, req[2]);
    // Send right
    int right_rank = ijtorank(segi, segj + 1, sq);
    MPI_Isend(tc2, nrows, MPI_DOUBLE, right_rank, step, MPI_COMM_WORLD, req[3]);
}

void recv_from_neig(int step, int rank, int sq, int segi, int segj, int nrows, int ncols, MPI_Request *req) {
    // Get from above
    int above_rank = ijtorank(segi - 1, segj, sq);
    MPI_Irecv(neig_u, ncols, MPI_DOUBLE, above_rank, step, MPI_COMM_WORLD, req[0]);
    // Get from below
    int below_rank = ijtorank(segi + 1, segj, sq);
    MPI_Irecv(neig_d, ncols, MPI_DOUBLE, below_rank, step, MPI_COMM_WORLD, req[1]);
    // Get from left
    int left_rank = ijtorank(segi, segj - 1, sq);
    MPI_Irecv(neig_l, nrows, MPI_DOUBLE, left_rank, step, MPI_COMM_WORLD, req[2]);
    // Get from right
    int right_rank = ijtorank(segi, segj + 1, sq);
    MPI_Irecv(neig_r, nrows, MPI_DOUBLE, right_rank, step, MPI_COMM_WORLD, req[3]);
}

void update(int step, int rank, int sq, int segi, int segj, int nrows, int ncols) {
    // We can at most have 4 neigs
    MPI_Request send_req[4];
    MPI_Request recv_req[4];
    send_to_neig(step, rank, sq, segi, segj, nrows, ncols, &send_req[0]);
    recv_from_neig(step, rank, sq, segi, segj, nrows, ncols, &recv_req[0]);
    // Update internal
    for (int i = 1; i < nrows-1; i++) {
        for (int j = 1; j < ncols-1; j++) {
            int neighbors = count_neighbors(i, j);
            if (grid[i][j] == 1 && (neighbors < 2 || neighbors > 3)) {
                new_grid[i][j] = 0;
            } else if (grid[i][j] == 0 && neighbors == 3) {
                new_grid[i][j] = 1;
            } else {
                new_grid[i][j] = grid[i][j];
            }
        }
    }
    // Check if we sent to neig
    MPI_Waitall(4, send_req, MPI_STATUSES_IGNORE);
    // Check if we recv:ed from neig
    MPI_Waitall(4, recv_req, MPI_STATUSES_IGNORE);
    // Update edges
    // Update grid
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            grid[i][j] = new_grid[i][j];
        }
    }
}

void write_output(int step) {
    char filename[50];
    sprintf(filename, "gol_output_%d.txt", step);
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(f, "%d ", grid[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

int perfect_square_root(int n) {
    for (int i = 2; i < n; i++) {
        if (i * i == n) return i;
    }
    return -1;
}

int main(int argc, char** argv) {
    int rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int sq = perfect_square_root(nprocs);
    assert(sq != -1);
    // Figure out which segment of the grid to work on
    int segi = rank / sq;
    int segj = rank % sq;
    // Figure out segment size
    int N = 20;
    int nrows = N / sq;
    int ncols = N / sq;
    int extraforlast = N % sq;
    if (segi == sq-1) nrows += extraforlast;
    if (segj == sq-1) ncols += extraforlast;

    initialize(rank, sq, segi, segj, nrows, ncols);
    for (int step = 0; step < STEPS; step++) {
        update(step, rank, sq, segi, segj, nrows, ncols);
        if (step % 10 == 0) write_output(step);
    }

    MPI_Finalize();
    destroy(nrows);
    return 0;
}
