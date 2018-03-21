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
    double array[(int)n];
    int nr = 0;
    double * ar;
    for (int i = 0; i < n; i++)
    {
        array[i] = i + i * 0.1;
    }
    if (world_rank == 0)
    {
        MPI_Send(&n, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        for (int i = 0; i < n; i++)
        {            
            MPI_Send(&array[i], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            printf(">>> a[%i]= %f \n", i, array[i]);

            MPI_Recv(&array[i], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf(">>> a[%i]= %f \n", i, array[i]);
        }
    }
    else
    {
        
        MPI_Recv(&nr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 1 received number %i from process 0\n", nr);
        ar = (double *) malloc((int)nr * sizeof(double));
        for (int i = 0; i < nr; i++)
        {   
            printf("<<<ar[%i]= ", i);
            MPI_Recv(&ar[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("%f \n", ar[i]);
            ar[i] = - ar[i];
            MPI_Send(&ar[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            printf("*** ar[%i]= %f \n", i, ar[i]);
        }
        
    }
    printf("-------\n");
    // if (world_rank == 0)
    // {
    //     for (int i = 0; i < n; i++)
    //     {            
    //         MPI_Recv(&array[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         printf(">>> a[%i]= %f \n", i, array[i]);
    //     }
    // }
    // else
    // {        
    //     for (int i = 0; i < nr; i++)
    //     {   
    //         ar[i] = - ar[i];
    //         MPI_Send(&ar[i], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    //         printf("*** ar[%i]= %f \n", i, ar[i]);
    //     }
        
    // }
    MPI_Finalize();
    return 0;
}