#include <stdlib.h>
#include <stdio.h>
//#include <mpi.h>

double* generateMatrix(int n, int sparsity);

int main(int argc, char* argv[]){
	double* matrix = generateMatrix(16, 1);
	for(int i = 0; i < 16; i++){
		for(int j = 0; j < 16; j++){
			printf("%lf ", matrix[i*16+j]);
			if(j==15) printf("\n");

		}
	}
	free(matrix);
}


double* generateMatrix(int n, int sparsity){

	double* matrix = (double*) calloc(n*n,sizeof(double));
	for(int i = 0; i < n; i++){
		for(int j = i; j < n; j++){
			matrix[i*n+j] += 1;
			matrix[j*n+i] += 1;

		}
	}
	return matrix;
}
