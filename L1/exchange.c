/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int rank, size;
  double a, b;
  MPI_Status status;

  MPI_Request request;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  
  a = 100.0 + (double) rank;  /* Different a on different processors */

  int recieve_from = (rank == 0) ? size-1: rank-1;

  MPI_Isend(&a, 1, MPI_DOUBLE, (rank+1)%size, 1, MPI_COMM_WORLD, &request);
  MPI_Irecv(&b, 1, MPI_DOUBLE, (rank-1)%size, 1, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, MPI_STATUS_IGNORE);

  printf("Processor %i got %f from processor %i\n",rank, b, recieve_from);







/*
 /// General implementation of blocking
	MPI_Send(&a, 1, MPI_DOUBLE, (rank+1)%size, 1, MPI_COMM_WORLD);
	MPI_Recv(&b, 1, MPI_DOUBLE, (rank-1)%size , 1, MPI_COMM_WORLD, &status);
	printf("Processor %i got %f from processor %i\n",rank, b, size-1);
*/


/* /// For three processes
  if (rank == 0) {
    MPI_Send(&a, 1, MPI_DOUBLE, 1, 111, MPI_COMM_WORLD);
//    MPI_Recv(&b, 1, MPI_DOUBLE, 1, 222, MPI_COMM_WORLD, &status);
    MPI_Recv(&b, 1, MPI_DOUBLE, 2, 333, MPI_COMM_WORLD, &status);

    printf("Processor 0 got %f from processor 2\n", b);
  } else if (rank==1) {
    MPI_Recv(&b, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD, &status);
//    MPI_Send(&a, 1, MPI_DOUBLE, 0, 222, MPI_COMM_WORLD);
    MPI_Send(&a, 1, MPI_DOUBLE, 2, 222, MPI_COMM_WORLD);

    printf("Processor 1 got %f from processor 0\n", b);
  } else if (rank == 2){
    MPI_Recv(&b, 1, MPI_DOUBLE, 1, 222, MPI_COMM_WORLD, &status);

    MPI_Send(&a, 1, MPI_DOUBLE, 0, 333, MPI_COMM_WORLD);

    printf("Processor 2 got %f from processor 1\n", b);

  }

*/



  MPI_Finalize(); 

  return 0;
}
