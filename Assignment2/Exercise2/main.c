#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 1000000000  // Large number of iterations

int main() {
    double *A, *B, *C;
    A = (double*) malloc(N * sizeof(double));
    B = (double*) malloc(N * sizeof(double));
    C = (double*) malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        A[i] = i * 1.0;
        B[i] = i * 2.0;
    }

    clock_t start = clock();
    for (int i = 0; i < N; i++) {
        C[i] = A[i] + B[i];  // Simple vector addition
    }
    clock_t end = clock();

    double execution_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", execution_time);

    // To avoid smart compilers skipping calculations,
    // we calculate the sum of the C-vector
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += C[i];
    }
    printf("Sum of C: %f\n", sum);

    free(A);
    free(B);
    free(C);
    return 0;
}
