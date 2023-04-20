
#define PRODUCE_OUTPUT_FILE

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/**
 * Read function data from an input file and store in an array. The input file
 * is supposed to contain an integer representing the number of function values,
 * followed by the function values. All values should be separated by white
 * spaces.
 * @param file_name Name of input file
 * @param values Pointer to array where the values are to be stored
 * @return 0 on success, -1 on error
 */
int read_input(const char *file_name, float **values, float **values2);

/**
 * Write function data to a file, with 4 decimal places. The values are
 * separated by spaces.
 * @param file_name Name of output file
 * @param output Function values to write
 * @param num_values Number of values to print
 * @return 0 on success, -1 on error
 */
int write_output(char *file_name, const float *output, int num_values);



int main(int argc, char **argv) {

	// Read input file
	if (3 != argc) {
		printf("Usage: matrix-matrixmult input_file output_file \n");
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];
  int N = 16;
  const int COLS = sqrt(N);
  const int ROWS = sqrt(N);
  
  MPI_Comm Cycle_Communication;
  MPI_Status status;

  
  int shift;
  int dims[2];
  int periods[2];
  int left,right,up,down;
  float start,end;
    
  
  MPI_Init(&argc, &argv);
  int p, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  

  const int NPROWS=sqrt(p);  /* number of rows in _decomposition_ */
  const int NPCOLS=sqrt(p);  /* number of cols in _decomposition_ */
  const int BLOCKROWS = ROWS/NPROWS;  /* number of rows in _block_ */
  const int BLOCKCOLS = COLS/NPCOLS; /* number of cols in _block_ */


  ////READING
  float *input;
  float *input2;
  if (rank == 0) {
 
	  int size_of_matrix;
	  if (0 > (size_of_matrix = read_input(input_name, &input, &input2))) {
		  return 2;
	  }
  }

  //ERRORS
  if (p != NPROWS*NPCOLS) {
      fprintf(stderr,"Error: number of PEs %d != %d x %d\n", p, NPROWS, NPCOLS);
      MPI_Finalize();
      exit(-1);
  }
  
  
  //ALLOCATING MEMORY ON EACH PROCESSOR
  float *A = (float*)malloc(sizeof(float)*BLOCKROWS*BLOCKCOLS);

  float *B = (float*)malloc(sizeof(float)*BLOCKROWS*BLOCKCOLS);
  
  float *result = (float*)calloc(BLOCKROWS*BLOCKCOLS,sizeof(float));
  memset(result, 0, BLOCKROWS*BLOCKCOLS*sizeof(float));

  //Creating vectorized datatype
  MPI_Datatype blocktype;
  MPI_Datatype blocktype2;
  MPI_Type_vector(BLOCKROWS, BLOCKCOLS, COLS, MPI_FLOAT, &blocktype2);
  MPI_Type_create_resized( blocktype2, 0, sizeof(float), &blocktype);
  MPI_Type_commit(&blocktype);

 
  //Create distribution thingymabob
  int disps[NPROWS*NPCOLS];
  int counts[NPROWS*NPCOLS];
  for (int ii=0; ii<NPROWS; ii++) {
    for (int jj=0; jj<NPCOLS; jj++) {
      disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;
      counts [ii*NPCOLS+jj] = 1;
    }
  }

  //Distribute data between processors
  MPI_Scatterv(input, counts, disps, blocktype, A, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(input2, counts, disps, blocktype, B, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, 0, MPI_COMM_WORLD);
//  printf("You are done with scattering matrix on proc %d \n ",rank);
 
    //Create the square setup of processors
  dims[0]=0; dims[1]=0;
  periods[0]=1; periods[1]=1;
  MPI_Dims_create(p,2,dims); 
//  printf("dims %d, %d : \n", dims[0], dims[1]);
  if(dims[0]!=dims[1]) {
    if(rank==0) printf("The number of processors must be a square.\n");
    MPI_Finalize();
    return 0;
  }
  //Start the shifting, creating the sprt p x sqrt p size grid for communication
  MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&Cycle_Communication);

  
  
  MPI_Cart_shift(Cycle_Communication,0,1,&left,&right);
  MPI_Cart_shift(Cycle_Communication,1,1,&up,&down);
  
//  printf("You are done with creating cart shift on proc %d \n ",rank);
  start=MPI_Wtime();
  for(shift=0;shift<dims[0];shift++) {
  
    for (int proc=0; proc<p; proc++) {
        if (proc == rank) {
            printf("Rank = %d\n", rank);
            
            
            if (rank == 0) {
                printf("Global matrix a: \n");
                for (int ii=0; ii<ROWS; ii++) {
                    for (int jj=0; jj<COLS; jj++) {
                        printf("%lf ",input[ii*COLS+jj]);
                    }
                    printf("\n");
                }
            }
	    if (rank == 0) {
                printf("Global matrix b: \n");
                for (int ii=0; ii<ROWS; ii++) {
                    for (int jj=0; jj<COLS; jj++) {
                        printf("%lf ",input2[ii*COLS+jj]);
                    }
                    printf("\n");
                }
            }

            printf("Local Matrix a at rank:%d\n",rank);
            for (int ii=0; ii<BLOCKROWS; ii++) {
                for (int jj=0; jj<BLOCKCOLS; jj++) {
                    printf("%lf ",A[ii*BLOCKCOLS+jj]);
                }
                printf("\n");
            }
            
            printf("Local matrix b at rank:%d\n",rank);
            for (int ii=0; ii<BLOCKROWS; ii++) {
                for (int jj=0; jj<BLOCKCOLS; jj++) {
                    printf("%lf ",B[ii*BLOCKCOLS+jj]);
                }
                printf("\n");
            }
            
/*            printf("\n");
            printf("Local result:\n");
            for (int ii=0; ii<BLOCKROWS; ii++) {
                for (int jj=0; jj<BLOCKCOLS; jj++) {
                    printf("%lf ",result[ii*BLOCKCOLS+jj]);
                }
                printf("\n");
            }*/
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
  
  printf("SHIFT: %d, proc %d\n",shift,rank);
  // Matrix multiplication
    for(int i=0;i<BLOCKROWS;i++){
      for(int k=0;k<BLOCKROWS;k++){
        for(int j=0;j<BLOCKROWS;j++){
          result[i*BLOCKROWS+j]+=A[i*BLOCKROWS+k]*B[k*BLOCKROWS+j];
        }
      }
      printf("\n");
    }
      
    
//    printf("test\n");

	printf("Result of rank %i with shift %i\n", rank,shift);
	for (int ii=0; ii<BLOCKROWS; ii++) {
		for (int jj=0; jj<BLOCKCOLS; jj++) {
			printf("%lf ", result[ii*BLOCKCOLS+jj]);
		}
		printf("\n");
	}
    if(shift==dims[0]-1){
	

	 break;
    }
  // Communication
    MPI_Sendrecv_replace(B,BLOCKROWS*BLOCKCOLS,MPI_FLOAT,left,1,right,1,Cycle_Communication,&status);

    MPI_Sendrecv_replace(A,BLOCKROWS*BLOCKCOLS,MPI_FLOAT,up,2,down,2,Cycle_Communication,&status);

    

  }
   //Finalize

  end=MPI_Wtime();  
  
    /* each proc prints it's "b" out, in order */ //Test printof A B
/*
    for (int proc=0; proc<p; proc++) {
        if (proc == rank) {
            printf("Rank = %d\n", rank);
            
            
            if (rank == 0) {
                printf("Global matrix a: \n");
                for (int ii=0; ii<ROWS; ii++) {
                    for (int jj=0; jj<COLS; jj++) {
                        printf("%lf ",input[ii*COLS+jj]);
                    }
                    printf("\n");
                }
            }
	    if (rank == 0) {
                printf("Global matrix b: \n");
                for (int ii=0; ii<ROWS; ii++) {
                    for (int jj=0; jj<COLS; jj++) {
                        printf("%lf ",input2[ii*COLS+jj]);
                    }
                    printf("\n");
                }
            }

            printf("Local Matrix a at rank:%d\n",rank);
            for (int ii=0; ii<BLOCKROWS; ii++) {
                for (int jj=0; jj<BLOCKCOLS; jj++) {
                    printf("%lf ",A[ii*BLOCKCOLS+jj]);
                }
                printf("\n");
            }
            
            printf("Local matrix b at rank:%d\n",rank);
            for (int ii=0; ii<BLOCKROWS; ii++) {
                for (int jj=0; jj<BLOCKCOLS; jj++) {
                    printf("%lf ",B[ii*BLOCKCOLS+jj]);
                }
                printf("\n");
            }
            
            printf("\n");
            printf("Local result:\n");
            for (int ii=0; ii<BLOCKROWS; ii++) {
                for (int jj=0; jj<BLOCKCOLS; jj++) {
                    printf("%lf ",result[ii*BLOCKCOLS+jj]);
                }
                printf("\n");
            }
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

  */
  
  float *rbuf;
  if(rank == 0) { 
    rbuf = (float *)malloc(BLOCKROWS*BLOCKCOLS*p*sizeof(float)); 
  }
  
  
  MPI_Barrier(Cycle_Communication);
  MPI_Gather(result, BLOCKROWS*BLOCKCOLS, MPI_FLOAT, rbuf ,BLOCKROWS*BLOCKCOLS, MPI_FLOAT, 0, Cycle_Communication);

  
  
   printf("Rank = %d\n", rank);
   if (rank == 0) {
     printf("Global result: \n");
     for (int ii=0; ii<ROWS; ii++) {
       for (int jj=0; jj<COLS; jj++) {
         printf("%f ",rbuf[ii*COLS+jj]);
       }
       printf("\n");
     }
   } 
   
  end=MPI_Wtime();
  if(rank==0){
   printf("Time: %.4fs\n",end-start);
  }
  if(rank == 0) { 
    free(rbuf);
  }  
  free(A);
 // printf("testA\n");
  free(B);
 // printf("testB\n");
  //free(buf);
  //printf("testbuf\n");
  //free(tmp);
  //printf("testtmp\n");
  free(result);
  
  //printf("testresult\n");
  MPI_Finalize();
  return 0;
}



int read_input(const char *file_name, float **values,float **values2) {
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
	if (NULL == (*values = malloc(num_values * sizeof(float)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
 	if (NULL == (*values2 = malloc(num_values * sizeof(float)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%f", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
  for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%f", &((*values2)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
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
