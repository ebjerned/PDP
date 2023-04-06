#include "stencil.h"
#include <string.h>
#include <assert.h>
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
	double* output = NULL;
	double ex_time[size];

	const int partition = num_values / size;
	double* sub_output = (double*)malloc(partition*sizeof(double));
	double* sub_input = (double*)malloc(partition*sizeof(double));

	if (rank == 0){
		output = (double*)malloc(num_values * sizeof(double));
		if(NULL == output) {
			perror("Couldn't allocate memory for output");
			return 2;
		}
	}
	
	double l_padding[2];
	double r_padding[2];

	MPI_Scatter(input, partition, MPI_DOUBLE, sub_input, partition, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double start = MPI_Wtime();


	for (int s = 0; s < num_steps; s++) {
		double l_edge_data[] = {sub_input[0], sub_input[1]};
		double r_edge_data[] = {sub_input[partition - 2], sub_input[partition - 1]};

		int u_dest = (rank+1)%size;
		int l_dest = rank == 0 ? size - 1 : rank-1;

		MPI_Send(l_edge_data, 2, MPI_DOUBLE, u_dest, 0, MPI_COMM_WORLD);
		MPI_Recv(r_padding, 2, MPI_DOUBLE, l_dest, 0, MPI_COMM_WORLD, &status);
		MPI_Send(r_edge_data, 2, MPI_DOUBLE, l_dest, 0, MPI_COMM_WORLD);
		MPI_Recv(l_padding, 2, MPI_DOUBLE, u_dest, 0, MPI_COMM_WORLD, &status);


		for (int i = 0; i < EXTENT; i++) {
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++) {
				int index = (i - EXTENT + j + partition) % partition;
				if (index < partition - 2){
					result += STENCIL[j] * sub_input[index];
				} else if (index == partition - 2) {
					result += STENCIL[j] * l_padding[0];
				} else if (index == partition - 1) {
					result += STENCIL[j] * l_padding[1];
				}
			}
			sub_output[i] = result;
		}

		for (int i = EXTENT; i < partition - EXTENT; i++) {
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++) {
				int index = i - EXTENT + j;
				result += STENCIL[j] * sub_input[index];
			}
			sub_output[i] = result;
		}

		for (int i = partition - EXTENT; i < partition; i++) {
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++) {
				int index = (i - EXTENT + j) % partition;
				if(index > 1){
					result += STENCIL[j] * sub_input[index];
				} else if (index == 0) {
					result += STENCIL[j] * r_padding[0];
				} else if (index == 1) {
					result += STENCIL[j] * r_padding[1];
				}
			}
			sub_output[i] = result;
		}

		if (s < num_steps - 1) {
			double* tmp = sub_input;
			sub_input = sub_output;
			sub_output = tmp;
		}
	}

	double my_execution_time = MPI_Wtime() - start;

	MPI_Gather(sub_output, partition, MPI_DOUBLE, output, partition, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	double max_time = 0;


	MPI_Reduce(&my_execution_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank== 0){
		// Write result
		printf("%f ", max_time);
	}





#ifdef PRODUCE_OUTPUT_FILE
	if(rank == 0){
		if (0 != write_output(output_name, output, num_values)) {
			return 2;
		}
	}
#endif

	// Clean up
	MPI_Finalize();
	free(sub_input);
	free(input);
	free(sub_output);
	free(output);

	return 0;
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
