
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>

int main(int argc, char* argv[]){
	
	if(argc != 2){
		perror("Usage: ./cg n{side length of mesh}\n");
		return -1;
	}
	

	MPI_Init(&argc, &argv);
	MPI_Comm Cycle_Communication;
	MPI_Status status;
	MPI_Request request1;
	MPI_Request request2;
	MPI_Request request3;
	MPI_Request request4;

	int shift;
	int dims[2];
	int periods[2];
	int left,right,up,down;
	double start,end;
	int n = atoi(argv[1]);
	double h = 1.0/(n+1);
	int p, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n_p = (int)sqrt(p);
	int sideElementsPerProc = n/n_p;
	
	int elementsPerProc = sideElementsPerProc*sideElementsPerProc;
	MPI_Barrier(MPI_COMM_WORLD);
	double* local_d = (double*)calloc(elementsPerProc, sizeof(double));
	double* local_q = (double*)calloc(elementsPerProc,sizeof(double));
	double* local_u = (double*)calloc(elementsPerProc, sizeof(double));
	double* local_g = (double*)calloc(elementsPerProc, sizeof(double));

	double* topDest = (double*)calloc(sideElementsPerProc, sizeof(double));
	double* bottomDest = (double*)calloc(sideElementsPerProc, sizeof(double));
	double* leftDest = (double*)calloc(sideElementsPerProc, sizeof(double));
	double* rightDest = (double*)calloc(sideElementsPerProc, sizeof(double));

	dims[0]=sqrt(p);
	dims[1]=sqrt(p);
	periods[0]=0;
	periods[1]=0;
	
	int mycoords[2];

	/* Setting up communication grid */
	MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&Cycle_Communication);
	MPI_Comm_rank(Cycle_Communication, &rank); 
	MPI_Cart_coords(Cycle_Communication, rank, 2, mycoords); 
	MPI_Cart_shift(Cycle_Communication,1, -1,&left,&right);
	MPI_Cart_shift(Cycle_Communication,0, -1,&down,&up);
	printf("PE %i has coordinate [%i, %i]\n", rank, mycoords[0], mycoords[1]);
	double q0;
	double DQ;
	double q1;
	double q0_sub = 0;
	start = MPI_Wtime();
	for(int i = 0; i < elementsPerProc; i++){
		double x_i = mycoords[1]*sideElementsPerProc+i%sideElementsPerProc;
		double y_i = mycoords[0]*sideElementsPerProc+i/sideElementsPerProc;
		//if(x_i != 0 && x_i != n-1 && y_i != 0 && y_i != n-1){	
			local_d[i] =  y_i;//2*h*h*h*(x_i*(1-h*x_i)+y_i*(1-h*y_i));
		//} else {
		//	local_d[i] = 0;
		//}
		q0_sub += local_d[i]*local_d[i];
		local_g[i] = -local_d[i];
		//if(rank==8) printf("\tPE %i contains (%lf, %lf, %lf)\n",rank, x_i, y_i, local_d[i]);
	}
	

	MPI_Allreduce(&q0_sub, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	printf("First q0 %lf\n", q0);
	
	MPI_Datatype sideValues_t;
	MPI_Datatype sideValues;
	MPI_Type_vector(sideElementsPerProc, 1, sideElementsPerProc, MPI_DOUBLE, &sideValues);
	MPI_Type_create_resized(sideValues, 0, sizeof(double), &sideValues_t);
	MPI_Type_commit(&sideValues_t);
	
	// Iteratation start here
	for(int t = 0; t < 1; t++){

		double g_norm = 0;
		double q1_sub = 0;
		double DQ_sub = 0;

		
		

		//if(mycoords[0] != 0){
			MPI_Isend(&local_d[0], sideElementsPerProc, MPI_DOUBLE, up, 0, Cycle_Communication, &request1);
			printf("Rank %i sends data uppwards to %i \n", rank, up);		
		//}
		

		//if(mycoords[0] != n_p-1){
			MPI_Isend(&local_d[elementsPerProc-sideElementsPerProc], sideElementsPerProc, MPI_DOUBLE, down, 1, Cycle_Communication, &request2);
			printf("Rank %i sends data downwards to %i \n", rank, down);			
		/*if(mycoords[1] != 0)*/ MPI_Isend(&local_d[0], 1, sideValues_t, left, 2, Cycle_Communication, &request3);
			printf("Rank %i sends data left to %i \n", rank, down);
		/*if(mycoords[1] != n_p-1)*/ MPI_Isend(&local_d[sideElementsPerProc-1], 1, sideValues_t, right, 3, Cycle_Communication, &request4);
			printf("Rank %i sends data right to %i \n", rank, down);
		//}
		// Alla center element
		
		for(int i = 1; i < sideElementsPerProc-1; i++){
			for(int j = 1; j < sideElementsPerProc-1;j++){
				int index = i*sideElementsPerProc + j;
				
				local_q[index] = -local_d[index+1] - local_d[index-1] + 4*local_d[index]-local_d[index+sideElementsPerProc] -local_d[index-sideElementsPerProc];
				//local_q[index] = local_d[index];
			}
		}

		// Inner över
		if(mycoords[0] == 0){
			for(int i = 1; i < sideElementsPerProc-1; i++)
				local_q[i] = 0;
		}else{
			MPI_Recv(topDest, sideElementsPerProc, MPI_DOUBLE, down, 0, Cycle_Communication, &status);
			for(int i = 1; i < sideElementsPerProc-1; i++){
				//local_q[i] = local_d[i];
				//printf("%i received from below:%lf\n",rank, topDest[i]);
				local_q[i] = -local_d[i +1] - local_d[i-1] + 4*local_d[i] - local_d[i+sideElementsPerProc] -topDest[i];
				// res = höger vänster center nere uppe
			}
			
		}


		// Inner under
		
		if(mycoords[0] == n_p-1){
			for(int i = elementsPerProc-sideElementsPerProc+1; i < elementsPerProc-1; i++)
				local_q[i] = 0;
		} else {
			MPI_Recv(bottomDest, sideElementsPerProc, MPI_DOUBLE, up, 1, Cycle_Communication, &status);
			for(int i = elementsPerProc-sideElementsPerProc+1; i < elementsPerProc-1; i++){
				//local_q[i] = local_d[i];
				//printf("\t%i received from above:%lf\n",rank, bottomDest[i-(elementsPerProc-sideElementsPerProc)]);
				local_q[i] = -local_d[i+1] - local_d[i-1] + 4*local_d[i] - local_d[i+sideElementsPerProc] -bottomDest[i-(elementsPerProc-sideElementsPerProc)];
				//res = höger vänster center uppe nere
			}
		}
		

		


		// Vänster


		if(mycoords[1] == 0){
			for(int i = 0; i < sideElementsPerProc; i++){
				int index = i*sideElementsPerProc;		
				local_q[index] = 0;
			}
		} else {
			MPI_Recv(&leftDest[0], sideElementsPerProc, MPI_DOUBLE, right, 2, Cycle_Communication, &status);
			for(int i = 0; i < sideElementsPerProc; i++){
				int index = i*sideElementsPerProc;			
				//local_q[index] = rank;
				//printf("\t%i received from right:%lf\n",rank, leftDest[i]);
				local_q[index] = -local_d[index+1] - leftDest[i] + 4*local_d[index];
				
				if(i != 0){
					local_q[index] -= local_d[index-sideElementsPerProc];
				} else {
				
					local_q[index] -= bottomDest[i];
					continue;
				}

				if(i != sideElementsPerProc-1){
					local_q[index] -= local_d[index+sideElementsPerProc];
				} else {
					local_q[index] -= topDest[i];
					continue;
				}
			}
		}
		// Höger i
		if(mycoords[1] == n_p-1){
			for(int i = 0; i < sideElementsPerProc; i++){
				int index = (i+1)*sideElementsPerProc-1;
				local_q[index] = 0;
			}
		}else{
			MPI_Recv(&rightDest[0], sideElementsPerProc, MPI_DOUBLE, left, 3, Cycle_Communication, &status);
			for(int i = 0; i < sideElementsPerProc; i++){
				int index = (i+1)*sideElementsPerProc-1;
				//local_q[index] = rank;
				local_q[index] = -local_d[index-1] - rightDest[i] + 4*local_d[index];	

				//printf("\t%i received from left:%lf\n",rank, rightDest[i]);
				if(i != 0){
					local_q[index] -= local_d[index-sideElementsPerProc];
				} else {
				
					local_q[index] -= bottomDest[i];
					continue;
				}

				if(i != sideElementsPerProc-1){
					local_q[index] -= local_d[index+sideElementsPerProc];
				} else {
					local_q[index] -= topDest[i];
					continue;

				}

			}
		}
		for(int i = 0; i < elementsPerProc; i++){
			DQ_sub += local_q[i]*local_d[i];
			//printf("q %i %i %lf\n",rank, i, local_q[i]);
			//printf("d %lf\n", local_d[i]);
		}
		printf("DQ_sub %.10lf\n", DQ_sub);

		
		MPI_Allreduce(&DQ_sub, &DQ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		printf("DQ %.10lf\n", DQ);
		printf("q0 %lf\n", q0);
		double tau = q0/DQ;
		printf("Tau %lf\n", tau);
		
		
		for(int i = 0; i < elementsPerProc; i++){
			local_u[i] += tau*local_d[i];
			local_g[i] += tau*local_q[i];
			q1_sub += local_g[i]*local_g[i];
			
		}


		MPI_Allreduce(&q1_sub, &q1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		double beta = q1/q0;
		for(int i = 0; i < elementsPerProc; i++){
			//local_d[i] = -local_g[i]+beta*local_d[i];
		}
		q0 = q1;

		g_norm = sqrt(q1);
		if(rank==0)printf(" %lf\n", g_norm);
	}
	end = MPI_Wtime();
	printf("%lf\n", end-start);

	double* res = (double*) malloc(n*n*sizeof(double));
	MPI_Datatype blocktype;
	MPI_Datatype blocktype2;
	MPI_Type_vector(sideElementsPerProc, sideElementsPerProc, sideElementsPerProc*n_p, MPI_DOUBLE, &blocktype2);
	MPI_Type_create_resized(blocktype2, 0, sizeof(double), &blocktype);
	MPI_Type_commit(&blocktype);
	int count[p];
	int displ[p];
	for(int i = 0; i < p;i++){
		count[i] = 1;
		displ[i] = (i/n_p)*elementsPerProc*n_p + (i%n_p)*sideElementsPerProc;

	}
	MPI_Gatherv(local_q, elementsPerProc, MPI_DOUBLE, res, count, displ, blocktype, 0, MPI_COMM_WORLD);
	if(rank==0)
		for(int i = 0; i < n*n; i++)
			printf("%.10lf\n", res[i]);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Type_free(&sideValues_t);
	free(topDest);
	free(bottomDest);
	free(leftDest);
	free(rightDest);
	free(local_g);
	free(local_u);
	free(local_d);
	free(local_q);	
	
	MPI_Finalize();
}


int ValueInA(int j, int i, int N, int m) {
	if (i == j) return 4;
	else if ((i == j + 1 && i % m != 0) || (i == j - 1 && (i + 1) % m != 0)) return -1;
	else if (j >= m && i == j - m) return -1;
	else if (j <= N - m && i == j + m) return -1;
	else return 0;
}




