#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int left  = (rank - 1 + size) % size;
    int right = (rank + 1) % size;

    srand((unsigned)time(NULL) + 1000003u * (unsigned)rank);
    int mine = rand() % 100;
    printf("rank %d: drew %d\n", rank, mine);
    fflush(stdout);

    int theirs;
    MPI_Sendrecv(&mine,   1, MPI_INT, right, 0,
                 &theirs, 1, MPI_INT, left,  0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("rank %d: received %d from rank %d\n", rank, theirs, left);
    fflush(stdout);

    int total;
    MPI_Allreduce(&mine, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("total: %d\n", total);
    }

    MPI_Finalize();
    return 0;
}
