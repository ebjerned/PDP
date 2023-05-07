
//#define PRODUCE_OUTPUT_FILE

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, double **values, double **values2);
int write_output(char *file_name, const double *output, int num_values);



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
	double* input_A;
	double* input_B;
	if (rank == 0) {
		if (0 > (N = read_input(input_name, &input_A, &input_B))) {
			return 2;
		}
	}


	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	const int COLS = N;
	const int ROWS = N;
	const int NPROWS=sqrt(p);
	const int NPCOLS=sqrt(p);
	const int BLOCKROWS = ROWS/NPROWS;
	const int BLOCKCOLS = COLS/NPCOLS;

	double* A = (double*)malloc(sizeof(double)*BLOCKROWS*BLOCKCOLS);
	double* B = (double*)malloc(sizeof(double)*BLOCKROWS*BLOCKCOLS);
	double* local_result = (double*)calloc(BLOCKROWS*BLOCKCOLS,sizeof(double));

	double* global_result;

	/* Setting up datatypes for selecting blocks for sub-matricies */
	MPI_Datatype blocktype;
	MPI_Datatype blocktype2;
	MPI_Type_vector(BLOCKROWS, BLOCKCOLS, COLS, MPI_DOUBLE, &blocktype2);
	MPI_Type_create_resized( blocktype2, 0, sizeof(double), &blocktype);
	MPI_Type_commit(&blocktype);


	int disps[NPROWS*NPCOLS];
	int counts[NPROWS*NPCOLS];
	for (int ii=0; ii<NPROWS; ii++) {
		for (int jj=0; jj<NPCOLS; jj++) {
			disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;
			counts[ii*NPCOLS+jj] = 1;
		}
	}
	start=MPI_Wtime();

	/* Scattering sub-matricies*/
	MPI_Scatterv(input_A, counts, disps, blocktype, A, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(input_B, counts, disps, blocktype, B, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
				MPI_Sendrecv_replace(A, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, left, 3, right, 3, Cycle_Communication, &status);
			}
		}
		if(mycoords[1] == i){
			for(int j = 0;j<i;j++){
				MPI_Sendrecv_replace(B, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE,up, 3, down, 3, Cycle_Communication, &status);
			}
		}
	}
	
	int CACHEBLOCKSIZE;

	// L1 cache 32kB -> 4000 doubles -> ~60 a side
	int suitable_block = (BLOCKROWS/NPROWS) % 60;
	CACHEBLOCKSIZE = suitable_block == 0 ? 60 : suitable_block;
	if(BLOCKROWS % CACHEBLOCKSIZE != 0){
		MPI_Abort(MPI_COMM_WORLD, 1);
		return -1;
	}
	for(shift=0;shift<dims[0];shift++) {

		/* Sub-matrix multiplication*/
		/*for(int i=0;i<BLOCKROWS;i++){
			for(int k=0;k<BLOCKROWS;k++){
				for(int j=0;j<BLOCKROWS;j++){
					local_result[i*BLOCKROWS+j] += A[i*BLOCKROWS+k]*B[k*BLOCKROWS+j];
				}
			}
		}*/
		/*
		for(int i=0;i<BLOCKROWS;i++){
			for(int k=0;k<BLOCKROWS;k++){
				int a = A[i*BLOCKROWS+k];
				for(int j=0;j<BLOCKROWS;j++){
					local_result[i*BLOCKROWS+j] += a*B[k*BLOCKROWS+j];
				}
			}
		}*/
		
		int bi, bj, bk, i, j, k;
		for(bi=0; bi < BLOCKROWS; bi += CACHEBLOCKSIZE)	
			for(bk=0; bk < BLOCKROWS; bk += CACHEBLOCKSIZE)
				for(bj=0; bj < BLOCKCOLS; bj += CACHEBLOCKSIZE)	
					for(i=bi; i < bi+CACHEBLOCKSIZE; i++)
						for(k=bk; k < bk+CACHEBLOCKSIZE; k++){
						int a = A[i*BLOCKROWS+k];
							for(j=bj; j < bj+CACHEBLOCKSIZE; j++)						
								local_result[i*BLOCKROWS +j] += a*B[k*BLOCKROWS+j];
						}
		/* Send new matricies to neighbourgs, shifting in Cannon's algorithm */
		MPI_Sendrecv_replace(A, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, left, 1, right, 1, Cycle_Communication, &status);
		MPI_Sendrecv_replace(B, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, up, 1, down, 1, Cycle_Communication, &status); 
	}
//	MPI_Barrier(Cycle_Communication);

	free(A);
	free(B);
	if(rank == 0) {
		free(input_A);
		free(input_B);
		global_result = (double*)malloc(BLOCKROWS*BLOCKCOLS*p*sizeof(double));
	}

	MPI_Gatherv(local_result, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, global_result, counts, disps, blocktype, 0, MPI_COMM_WORLD);

	/*if (rank == 0) {
		printf("Global result: \n");
		for (int ii=0; ii<ROWS; ii++) {
			for (int jj=0; jj<COLS; jj++) {
				printf("%f ",global_result[ii*COLS+jj]);
			}
			printf("\n");
		}
	}*/

	end=MPI_Wtime();
//	MPI_Barrier(Cycle_Communication);
	if(rank==0){
		#ifdef PRODUCE_OUTPUT_FILE
		write_output(output_name, global_result, N);
		#endif
		printf("%.10f\n",end-start);
//		free(input_A);
//		free(input_B);
		free(global_result);

	}


	free(local_result);
//	free(A);
//	free(B);
	MPI_Type_free(&blocktype);
	MPI_Type_free(&blocktype2);
	MPI_Finalize();
	return 0;
}



int read_input(const char *file_name, double **values, double **values2) {
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
	num;
	if (NULL == (*values = malloc(num * num * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
 	if (NULL == (*values2 = malloc(num * num * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num*num; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
  for (int i=0; i<num*num; i++) {
		if (EOF == fscanf(file, "%lf", &((*values2)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num;
}


int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values*num_values; i++) {
		if (0 > fprintf(file, "%.6lf ", output[i])) {
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
