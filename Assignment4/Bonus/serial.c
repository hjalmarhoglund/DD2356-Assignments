/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-09
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 20  // Grid size
#define STEPS 100  // Simulation steps

int grid[N][N], new_grid[N][N];

void initialize() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            grid[i][j] = rand() % 2;  // Random initial state
        }
    }
}

int count_neighbors(int x, int y) {
    int sum = 0;
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0) continue;
            int nx = (x + i + N) % N;
            int ny = (y + j + N) % N;
            sum += grid[nx][ny];
        }
    }
    return sum;
}

void update() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
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
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
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

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    initialize();
    for (int step = 0; step < STEPS; step++) {
        update();
        if (step % 10 == 0) write_output(step);
    }

    MPI_Finalize();
    return 0;
}
