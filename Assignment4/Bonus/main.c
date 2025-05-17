/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-17
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#define STEPS 100  // Simulation steps

int **grid;
int **new_grid;
int *neig_l, *neig_r;
int *tc1, *tc2;

/*
To make calculations easier when we are updating
we wrap our grid with an extra row and col on each
side. That explains the +2 below.
*/
int **malloc_grid(int nrows, int ncols) {
    nrows += 2;
    ncols += 2;
    int **ret = (int **) malloc(nrows * sizeof(int *));
    for (int i = 0; i < nrows; i++) {
        ret[i] = (int *) malloc(ncols * sizeof(int));
    }
    return ret;
}

void malloc_neigs(int sq, int segi, int segj, int nrows, int ncols) {
    neig_l = (int *) malloc(nrows * sizeof(int));
    neig_r = (int *) malloc(nrows * sizeof(int));
}

void malloc_tmp_cols(int nrows) {
    tc1 = (int *) malloc(nrows * sizeof(int));
    tc2 = (int *) malloc(nrows * sizeof(int));
}

/*
Since we wrapped our grid in one extra row/column
we start at 1 and end at nrows/ncols
*/
void initialize(int rank, int sq, int segi, int segj, int nrows, int ncols) {
    grid = malloc_grid(nrows, ncols);
    new_grid = malloc_grid(nrows, ncols);
    malloc_neigs(sq, segi, segj, nrows, ncols);
    malloc_tmp_cols(nrows);
    for (int i = 1; i <= nrows; i++) {
        for (int j = 1; j <= ncols; j++) {
            grid[i][j] = rand() % 2;  // Random initial state
        }
    }
}

void destroy(int nrows) {
    for (int i = 0; i < nrows+2; i++) {
        free(grid[i]);
        free(new_grid[i]);
    }
    free(grid);
    free(new_grid);

    free(neig_l);
    free(neig_r);

    free(tc1);
    free(tc2);
}

int count_neighbors(int x, int y) {
    int sum = 0;
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0) continue;
            int nx = x + i;
            int ny = y + j;
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
    MPI_Isend(&grid[1][1], ncols, MPI_INT, above_rank, step, MPI_COMM_WORLD, &req[0]);
    // Send to below
    int below_rank = ijtorank(segi + 1, segj, sq);
    MPI_Isend(&grid[nrows][1], ncols, MPI_INT, below_rank, step, MPI_COMM_WORLD, &req[1]);
    // To send left and right, we move the columns to tc1 and tc2
    for (int i = 1; i <= nrows; i++) {
        tc1[i] = grid[i][1];
        tc2[i] = grid[i][ncols];
    }
    // Send left
    int left_rank = ijtorank(segi, segj - 1, sq);
    MPI_Isend(tc1, nrows, MPI_INT, left_rank, step, MPI_COMM_WORLD, &req[2]);
    // Send right
    int right_rank = ijtorank(segi, segj + 1, sq);
    MPI_Isend(tc2, nrows, MPI_INT, right_rank, step, MPI_COMM_WORLD, &req[3]);

    // Send diagonals
    int top_left_rank = ijtorank(segi - 1, segj - 1, sq);
    MPI_Isend(&grid[1][1], 1, MPI_INT, top_left_rank, step, MPI_COMM_WORLD, &req[4]);
    int top_right_rank = ijtorank(segi - 1, segj + 1, sq);
    MPI_Isend(&grid[1][ncols], 1, MPI_INT, top_right_rank, step, MPI_COMM_WORLD, &req[5]);
    int bottom_left_rank = ijtorank(segi + 1, segj - 1, sq);
    MPI_Isend(&grid[nrows][1], 1, MPI_INT, bottom_left_rank, step, MPI_COMM_WORLD, &req[6]);
    int bottom_right_rank = ijtorank(segi + 1, segj + 1, sq);
    MPI_Isend(&grid[nrows][ncols], 1, MPI_INT, bottom_right_rank, step, MPI_COMM_WORLD, &req[7]);
}

void recv_from_neig(int step, int rank, int sq, int segi, int segj, int nrows, int ncols, MPI_Request *req) {
    // Get from above
    int above_rank = ijtorank(segi - 1, segj, sq);
    MPI_Irecv(&grid[0][1], ncols, MPI_INT, above_rank, step, MPI_COMM_WORLD, &req[0]);
    // Get from below
    int below_rank = ijtorank(segi + 1, segj, sq);
    MPI_Irecv(&grid[nrows][1], ncols, MPI_INT, below_rank, step, MPI_COMM_WORLD, &req[1]);
    // Get from left
    int left_rank = ijtorank(segi, segj - 1, sq);
    MPI_Irecv(neig_l, nrows, MPI_INT, left_rank, step, MPI_COMM_WORLD, &req[2]);
    // Get from right
    int right_rank = ijtorank(segi, segj + 1, sq);
    MPI_Irecv(neig_r, nrows, MPI_INT, right_rank, step, MPI_COMM_WORLD, &req[3]);

    // Get diagonals
    int top_left_rank = ijtorank(segi - 1, segj - 1, sq);
    MPI_Irecv(&grid[0][0], 1, MPI_INT, top_left_rank, step, MPI_COMM_WORLD, &req[4]);
    int top_right_rank = ijtorank(segi - 1, segj + 1, sq);
    MPI_Irecv(&grid[0][ncols+1], 1, MPI_INT, top_right_rank, step, MPI_COMM_WORLD, &req[5]);
    int bottom_left_rank = ijtorank(segi + 1, segj - 1, sq);
    MPI_Irecv(&grid[nrows+1][0], 1, MPI_INT, bottom_left_rank, step, MPI_COMM_WORLD, &req[6]);
    int bottom_right_rank = ijtorank(segi + 1, segj + 1, sq);
    MPI_Irecv(&grid[nrows+1][ncols+1], 1, MPI_INT, bottom_right_rank, step, MPI_COMM_WORLD, &req[7]);
}

void do_update(int i, int j) {
    int neighbors = count_neighbors(i, j);
    if (grid[i][j] == 1 && (neighbors < 2 || neighbors > 3)) {
        new_grid[i][j] = 0;
    } else if (grid[i][j] == 0 && neighbors == 3) {
        new_grid[i][j] = 1;
    } else {
        new_grid[i][j] = grid[i][j];
    }
}

void update(int step, int rank, int sq, int segi, int segj, int nrows, int ncols) {
    // We have 8 neigs
    MPI_Request send_req[8];
    MPI_Request recv_req[8];
    send_to_neig(step, rank, sq, segi, segj, nrows, ncols, &send_req[0]);
    recv_from_neig(step, rank, sq, segi, segj, nrows, ncols, &recv_req[0]);
    // Update internal
    for (int i = 2; i < nrows; i++) {
        for (int j = 2; j < ncols; j++) {
            do_update(i,j);
        }
    }
    // Check if we sent to neig
    MPI_Waitall(8, send_req, MPI_STATUSES_IGNORE);
    // Check if we recv:ed from neig
    MPI_Waitall(8, recv_req, MPI_STATUSES_IGNORE);
    // We now update the edge columns with the information
    // we stored in the temporary buffers (tc1 & tc2)
    for (int i = 0; i < nrows; i++) {
        grid[0][i+1] = neig_l[i];
        grid[ncols+1][i+1] = neig_r[i];
    }
    // Update edges
    for (int j = 1; j <= ncols; j++) {
        do_update(1,j);
        do_update(nrows,j);
    }
    for (int i = 2; i < nrows; i++) {
        do_update(i,1);
        do_update(i,ncols);
    }
    // Update grid
    for (int i = 1; i <= nrows; i++) {
        for (int j = 1; j <= ncols; j++) {
            grid[i][j] = new_grid[i][j];
        }
    }
}

void write_output(int step, int N) {
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
        //if (step % 10 == 0) write_output(step);
    }

    MPI_Finalize();
    destroy(nrows);
    return 0;
}
