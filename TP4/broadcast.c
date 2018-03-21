#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    int i = 0;
    int count = 0;
    int n = 10;
    int array[(int)n];
    for (int i = 0; i < n; i++)
    {
        array[i] = 0;
    }

    if (world_rank == 0)
    {
        for (int i = 0; i < n; i++)
        {
            array[i] = i;
        }
    }
    else
    {
        printf("%d \n", world_rank);
        for (int i = 0; i < n; i++)
        {
            array[i] = 0;
            printf("%i ", i, array[i]);            
        }
        printf("\n");
    }
    printf("--------\n");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(array, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank != 0)
    {
        printf("%d \n", world_rank);
        for (int i = 0; i < n; i++)
        {
            printf("%i ", i, array[i]);            
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}