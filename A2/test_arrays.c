
//#define PRODUCE_OUTPUT_FILE

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, float **values, float **values2, int* num_values);
int write_output(char *file_name, const float *output, int num_values);



int main(int argc, char **argv) {

	if (3 != argc) {
		printf("Usage: matrix-matrixmult input_file output_file \n");
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];

	MPI_Comm Cycle_Communication;
	MPI_Status status;


	int shift;
	int dims[2];
	int periods[2];
	int left,right,up,down;
	double start,end;


	MPI_Init(&argc, &argv);
	int p, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int N;
	float* input_A;
	float* input_B;
	if (rank == 0) {
		int size_of_matrix;
		if (0 > (size_of_matrix = read_input(input_name, &input_A, &input_B, &N))) {
			return 2;
		}
	}

	start=MPI_Wtime();


	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	const int COLS = sqrt(N);
	const int ROWS = sqrt(N);
	const int NPROWS=sqrt(p);
	const int NPCOLS=sqrt(p);
	const int BLOCKROWS = ROWS/NPROWS;
	const int BLOCKCOLS = COLS/NPCOLS;

	float* A = (float*)malloc(sizeof(float)*BLOCKROWS*BLOCKCOLS);
	float* B = (float*)malloc(sizeof(float)*BLOCKROWS*BLOCKCOLS);
	float* local_result = (float*)calloc(BLOCKROWS*BLOCKCOLS,sizeof(float));

	float* global_result;

	/* Setting up datatypes for selecting blocks for sub-matricies */
	MPI_Datatype blocktype;
	MPI_Datatype blocktype2;
	MPI_Type_vector(BLOCKROWS, BLOCKCOLS, COLS, MPI_FLOAT, &blocktype2);
	MPI_Type_create_resized( blocktype2, 0, sizeof(float), &blocktype);
	MPI_Type_commit(&blocktype);


	int disps[NPROWS*NPCOLS];
	int counts[NPROWS*NPCOLS];
	for (int ii=0; ii<NPROWS; ii++) {
		for (int jj=0; jj<NPCOLS; jj++) {
			disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;
			counts[ii*NPCOLS+jj] = 1;
		}
	}

	/* Scattering sub-matricies*/
	MPI_Scatterv(input_A, counts, disps, blocktype, A, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(input_B, counts, disps, blocktype, B, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, 0, MPI_COMM_WORLD);
	if(rank==0){
//		free(input_A);
//		free(input_B);
	}

	dims[0]=sqrt(p);
	dims[1]=sqrt(p);
	periods[0]=1;
	periods[1]=1;
	
	int mycoords[2];

	/* Setting up communication grid */
	MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&Cycle_Communication);
	MPI_Comm_rank(Cycle_Communication, &rank); 
	MPI_Cart_coords(Cycle_Communication, rank, 2, mycoords); 
	MPI_Cart_shift(Cycle_Communication,1,-1,&right,&left);
	MPI_Cart_shift(Cycle_Communication,0,-1,&down,&up);


	/* Skew matricies */
	for(int i = 1; i < p; i++ ){
		if(mycoords[0] == i){
			for(int j = 0; j<i;j++){
				MPI_Sendrecv_replace(A, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, left, 3, right, 3, Cycle_Communication, &status);
			}
		}
		if(mycoords[1] == i){
			for(int j = 0;j<i;j++){
				MPI_Sendrecv_replace(B, BLOCKROWS*BLOCKCOLS, MPI_FLOAT,up, 3, down, 3, Cycle_Communication, &status);
			}
		}
	}

	for(shift=0;shift<dims[0];shift++) {

		/* Sub-matrix multiplication*/
		for(int i=0;i<BLOCKROWS;i++){
			for(int k=0;k<BLOCKROWS;k++){
				for(int j=0;j<BLOCKROWS;j++){
					local_result[i*BLOCKROWS+j] += A[i*BLOCKROWS+k]*B[k*BLOCKROWS+j];
				}
			}
		}
		/* Send new matricies to neighbourgs, shifting in Cannon's algorithm */
		MPI_Sendrecv_replace(A, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, left, 1, right, 1, Cycle_Communication, &status);
		MPI_Sendrecv_replace(B, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, up, 1, down, 1, Cycle_Communication, &status); 
	}

	free(A);
	free(B);
	MPI_Barrier(Cycle_Communication);
	if(rank == 0) {
		free(input_A);
		free(input_B);
		global_result = (float*)malloc(BLOCKROWS*BLOCKCOLS*p*sizeof(float));
	}

	MPI_Gatherv(local_result, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, global_result, counts, disps, blocktype, 0, MPI_COMM_WORLD);
/*
	if (rank == 0) {
		printf("Global result: \n");
		for (int ii=0; ii<ROWS; ii++) {
			for (int jj=0; jj<COLS; jj++) {
				printf("%f ",global_result[ii*COLS+jj]);
			}
			printf("\n");
		}
	}*/

	end=MPI_Wtime();
	if(rank==0){
		#ifdef PRODUCE_OUTPUT_FILE
		write_output(output_name, global_result, N);
		#endif
		printf("%.10f\n",end-start);
		free(global_result);
	}


	free(local_result);

	MPI_Finalize();
	return 0;
}



int read_input(const char *file_name, float **values, float **values2, int* num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num;
	if (EOF == fscanf(file, "%d", &num)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	*num_values = num;
	if (NULL == (*values = malloc(num * sizeof(float)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
 	if (NULL == (*values2 = malloc(num * sizeof(float)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num; i++) {
		if (EOF == fscanf(file, "%f", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
  for (int i=0; i<num; i++) {
		if (EOF == fscanf(file, "%f", &((*values2)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num;
}


int write_output(char *file_name, const float *output, int num_values) {
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
