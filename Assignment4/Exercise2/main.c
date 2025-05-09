/*
Authors: Christiane Kobalt & Hjalmar Höglund for DD2356 Assignment 4
Date: 2025-05-09
 */

#include <stdio.h>
#include <stdlib.h>

#define N 1000 // Matrix size

void initialize_matrix(double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = i + j * 0.01;
        }
    }
}

void compute_row_sums(double matrix[N][N], double row_sums[N]) {
    for (int i = 0; i < N; i++) {
        row_sums[i] = 0.0;
        for (int j = 0; j < N; j++) {
            row_sums[i] += matrix[i][j];
        }
    }
}

void write_output(double row_sums[N]) {
    FILE *f = fopen("row_sums_output.txt", "w");
    for (int i = 0; i < N; i++) {
        fprintf(f, "%f\n", row_sums[i]);
    }
    fclose(f);
}

int main() {
    double matrix[N][N], row_sums[N];
    initialize_matrix(matrix);
    compute_row_sums(matrix, row_sums);
    write_output(row_sums);
    printf("Row sum computation complete.\n");
    return 0;
}
