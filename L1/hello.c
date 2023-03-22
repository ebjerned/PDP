/**********************************************************************
 * A simple "hello world" program for MPI/C
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

  int rank, size;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(rank == 0) {
	printf("Hello World says process %i! Executing on %i processes\n", rank, size );

  } else {
	printf("Hello World says process %i!\n", rank);

  }
               /* Print a message              */

  MPI_Finalize();                       /* Shut down and clean up MPI   */

  return 0;
}
