#include <stdio.h>
#include <omp.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // OpenMP parallel region
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        printf("MPI rank %d: OpenMP thread %d\n", world_rank, thread_id);
    }

    MPI_Finalize();
    return 0;
}
