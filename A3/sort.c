#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int comparisonFunction(const void* a, const void* b);
int median(int* values, int n_values);
int findPivot(int* pivots, int processes, char method);
int main(int argc, char* argv[]){
	

	int p, rank;
	int* values;
	int* pivots;
	int pivotFromStrategy;
	int numberOfValues = 100;

	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0){
		values = (int*)malloc(numberOfValues*sizeof(int));
		pivots = (int*)malloc(p*sizeof(int));
		for(int i = 0; i < numberOfValues; i++)
			values[i] = rand() % 50;
	}
	int* localValues = (int*)malloc(numberOfValues/p*sizeof(int));

	MPI_Scatter(values, numberOfValues/p, MPI_INT, localValues, numberOfValues/p, MPI_INT, 0, MPI_COMM_WORLD);
	qsort(localValues, numberOfValues/p, sizeof(int), comparisonFunction);

	int localPivot = median(localValues, numberOfValues/p);
	MPI_Gather(&localPivot, 1, MPI_INT, pivots, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank == 0){
		pivotFromStrategy = findPivot(pivots, p, 0);

	}
	MPI_Bcast(&pivotFromStrategy, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Finalize();

}

int comparisonFunction(const void* a, const void* b){
	return (*(int*)a - *(int*)b);
}

int median(int* values, int n_values){
	if(n_values % 2 == 0)
		return (values[n_values/2] + values[n_values/2+1])/2;

	return values[n_values/2];
}

int mean(int* values, int n_values){
	int sum = 0;
	for(int i = 0; i < n_values; i++){
		sum += values[i];
	}
	return sum / n_values;
}

int findPivot(int* pivots, int processes, char method){
	switch(method){
		case 1:
			return pivots[0];
		case 2:
			qsort(pivots, processes, sizeof(int), comparisonFunction);
			return median(pivots, processes);
		case 3:
			return mean(pivots, processes);

	}

}

