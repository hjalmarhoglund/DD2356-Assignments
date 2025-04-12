#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define BUFSIZE 1024

int main(int argc, char *argv[]) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int data[BUFSIZE];
    for (int i = 0; i < BUFSIZE; i++) data[i] = 42;
    double start, end;
    for (int reruns = 0; reruns < 11; reruns++) {
        for (int size = 1; size <= 1024; size *= 2) {
            if (rank == 0) {
                start = MPI_Wtime();
                MPI_Send(&data, size, MPI_INT, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&data, size, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                end = MPI_Wtime();
                printf("Size: %lu Bytes, RTT: %f microseconds\n", size * sizeof(int), (end - start) * 1e6);
            } else {
                MPI_Recv(&data, size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&data, size, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();
    return 0;
}
