#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[]){

	if(argc != 3){
		perror("USAGE: input_path, output_path");
	}

	int size;
	int rank;
	 
	float* A;
	float* B;
	int* send_count;
	int* disp;
	int num_values;

	MPI_Init(&argc, &argv);
	MPI_Datatype sub_matrix_t;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;

	int side = sqrt((float)size);

	int x_coord = rank%side;
	int y_coord = rank == 0 ? 0 : rank/side;

	float* local_sub_matrix = (float*) malloc(1*1*sizeof(float));
	MPI_Type_vector(1, 1, 3, MPI_FLOAT, &sub_matrix_t);
	MPI_Type_commit(&sub_matrix_t);

	if(rank == 0){
		read_input(argv[1], &A, &B, &num_values);

		int mat_side = sqrt((float)num_values);
		int n_side = mat_side/side;

		send_count = (int*)malloc(size*sizeof(int));
		disp = (int*)malloc(size*sizeof(int));
//		printf("%i, %i, %i\n",num_values, mat_side, n_side);

		memset(send_count, 1, size*sizeof(int));
		int counter = 0;


		for(int i = 0; i < size; i++){
			int y = i == 0 ? 0 : i/side;
			int x = i%side;
			disp[i] = y*n_side*mat_side+x*n_side;
//			printf("%i\n", disp[i]);
		}

/*
		t[0] = vtype;
		t[1] = MPI_UB;
		displs[0] = 0;
		displs[1] = 1 * sizeof(float);
		blklen[0] = 1;
		blklen[1] = 1;
		MPI_Type_struct( 2, blklen, displs, t, &stype );
		MPI_Type_commit( &stype );
		MPI_Scatterv(&A[0], sendcount, disp, stype, &local_sub_matrix[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD );
*/
	}
	
	if(rank==0){
		MPI_Scatterv(A, send_count, disp, sub_matrix_t, local_sub_matrix, 1, sub_matrix_t, 0, MPI_COMM_WORLD);
	} else {
//		MPI_Scatterv( (void *)0, (void *)0, (void *)0, MPI_DATATYPE_NULL, local_sub_matrix, 1, MPI_FLOAT, 0, MPI_COMM_WORLD );
	}

	printf("PE: %i, (%i, %i), %f\n", rank,x_coord, y_coord, local_sub_matrix[0]);



	MPI_Finalize();

}

// TODO: Byt ut till Mattias variant

int read_input(const char *file_name, float** A, float** B, int* num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num;
	if (EOF == fscanf(file, "%d", num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	printf("%i", *num_values);
	if (NULL == (*A = malloc(*num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	if (NULL == (*B = malloc(*num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<*num_values; i++) {
		if (EOF == fscanf(file, "%f", &((*A)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}

	for (int i=0; i<*num_values; i++) {
		if (EOF == fscanf(file, "%f", &((*B)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return *num_values;
}


int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}

