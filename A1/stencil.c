#include "stencil.h"
#include <string.h>

int main(int argc, char **argv) {
	if (4 != argc) {
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];
	int num_steps = atoi(argv[3]);

	// Read input file
	double *input;
	int num_values;
	if (0 > (num_values = read_input(input_name, &input))) {
		return 2;
	}

	// Stencil values
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};
	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int partition = num_values / size;

	double* sub_input = (double*)malloc((partition+4)*sizeof(double));

	if(rank == 0){
		if(size == 1){
			memcpy(&(sub_input[2]), input, num_values*sizeof(double));
			memcpy(sub_input, &(input[num_values-2]), 2*sizeof(double));
			memcpy(&(sub_input[num_values+2]), input, 2*sizeof(double));
			
		}else{
			for(int i = 0; i < size; i++){
				printf("PE: %i\n", i);
				if(i == 0){
					double* padded_input = (double*) malloc((partition+4)*sizeof(double));
					memcpy(sub_input, input, (partition+2)*sizeof(double));
					memcpy(&(sub_input[partition+2]), &(input[num_values-2]), 2*sizeof(double));

					printf("Rank 0 finished\n");


				} else if (i == size-1){
					double* padded_input = (double*) malloc((partition+4)*sizeof(double));

					memcpy(&(padded_input[2]), &(input[i*partition-2]), (partition+2)*sizeof(double));
					memcpy(padded_input, input, 2*sizeof(double));
					MPI_Ssend(padded_input,  partition+4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);


				} else {
					MPI_Ssend(&(input[i*partition-2]),  partition+4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);


				}

			}

		}

	} else {
		MPI_Recv(&sub_input, partition+4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
	}
	for(int i = 0; i < size; i++){
		for(int j = 0; j < partition+4; j++){
			printf("PE %i, %lf \n", i, sub_input[j]);
		}
	}
	// Start timer
	double  start = MPI_Wtime();
	// Allocate data for result
	double *output;
	if (NULL == (output = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for output");
		return 2;
	}

	// Repeatedly apply stencil
	for (int s=0; s<num_steps; s++) {

		// Apply stencil
		// This loop handles left hand side wrapping
/*
		for (int i=0; i<EXTENT; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = (i - EXTENT + j + num_values) % num_values;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}*/
		// This loop handles all values in the central entries with no special cases
		/*for (int i=EXTENT; i<num_values-EXTENT; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = i - EXTENT + j;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}*/
		for (int i=EXTENT; i < partition; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = i - EXTENT + j;
				result += STENCIL[j] * sub_input[index];
			}
			output[i-EXTENT] = result;
		}
		// This loop handles right hand side wrapping
/*		for (int i=num_values-EXTENT; i<num_values; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = (i - EXTENT + j) % num_value;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}*/
		// Swap input and output
		if (s < num_steps-1) {
			double *tmp = input;
			input = output;
			output = tmp;
		}
	}




	free(input);
	// Stop timer
	double my_execution_time = MPI_Wtime() - start;

	// Write result
	printf("%f\n", my_execution_time);




#ifdef PRODUCE_OUTPUT_FILE
	if (0 != write_output(output_name, output, num_values)) {
		return 2;
	}
#endif

	// Clean up
	free(output);
	MPI_Finalize();
	return 0;
}


void applyStencil(double* input_data, const int STENCIL_WIDTH, const int EXTENT, const double STENCIL[]){

}


int read_input(const char *file_name, double **values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
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
