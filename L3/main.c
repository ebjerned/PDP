#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[]){
	int size;
	int rank;
	
	MPI_Init(&argc, &argv);
	MPI_Datatype contigous;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;

	

	int n = 10;	
	int k = (n/size);

	int* local_array = (int*) malloc(k*n*sizeof(int));
	memset(local_array, rank, k*n*sizeof(int));
	int m = 2;
	int* recv_data = (int*) malloc(m*n*sizeof(int));
	MPI_Type_contiguous(m*n, MPI_INT, &contigous);
	MPI_Type_commit(&contigous);
	printf("Rank %i sending to %i\n", rank, (rank+1)%size);
	
	MPI_Sendrecv(&local_array[(k-m)*n], 1, contigous, (rank+1)%size, 0, &recv_data, 1, contigous, rank, 0, MPI_COMM_WORLD, &status);
	
	MPI_Finalize();
	return 0;

}
