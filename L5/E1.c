#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>



int main(int argc, char* argv[]){
	int rank, size;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Win win;
	
	double times[10*size];
	
	MPI_Win_create(times, 10*size*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
	MPI_Win_fence(0, win);	
	for(int iter = 0; iter < 10; iter++){
		srand(((long int)&rank) * (iter+1));
		int N = ((double)rand()/RAND_MAX)*700000000 + 100000000;
		
		double start = MPI_Wtime();
		long int sum = 0;
		for(int i = 0; i < N; i++){
			sum += (int)(sqrt(i*i)/sqrt(N));
		}
		double end = MPI_Wtime() - start;
		MPI_Put(&end, 1, MPI_DOUBLE, 0, size*iter+rank, 1, MPI_DOUBLE, win);
		
	}
	MPI_Win_fence(0, win);

	if(rank == 0){
		for(int i = 0; i < size*10; i++)
			printf("%.10lf\n", times[i]);
	}
	MPI_Finalize();


}
