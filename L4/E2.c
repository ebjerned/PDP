#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char* argv[]){
	
	int N = atoi(argv[1]);

	int k 
	int size;
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	int Np = N/size;
	int* x = (int*) malloc(Np*sizeof(int));
	int* y = (int*) malloc(Np*sizeof(int));
	int d = 0;
	for(int i = 0; i < Np; i++){
		x[i] = 1;
		y[i] = 1;
		d += x[i]*y[i];

	}

	MPI_Cart_create(MPI_COMM_WORLD, k
	
	
	MPI_Finalize();

}

