#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]){
	
	int N = atoi(argv[1]);

	
	int size;
	int rank;
	int* list;
	int* results;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0){
		list = (int*)malloc(N*sizeof(int));
		results = (int*)malloc(size*sizeof(int));
		srand(0);
		for(int i = 0; i < N; i++){
			list[i] = rand() % 100;
//			printf("%i\n", list[i]);

		}
	}

	int Np = N/size;
	int* local_list = (int*)malloc(Np*sizeof(int));
	MPI_Scatter(list, Np, MPI_INT, local_list, Np, MPI_INT, 0, MPI_COMM_WORLD);

	int found_at;
	for(int i = 0; i < Np; i++){
		if(!local_list[i]){
			found_at = i+Np*rank;
			break;
		}
	}

	free(local_list);
	MPI_Gather(&found_at, 1, MPI_INT, results, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank == 0){
		int min = 0;
		for(int i = 0; i < size; i++){
			if(results[i] != NULL){
				printf("%i\n", results[i]);
				break;
			}
			

		}
		free(list);
		free(results);
	}
	
	MPI_Finalize();

}

