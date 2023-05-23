
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
	MPI_Request request;


	int dims[2];
	int periods[2];
	int left,right,up,down;
	double start,end_time;
	int n = atoi(argv[1]);
	const double h = 1.0/(n+1);
	int p, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const int n_p = (int)sqrt(p);
	const int sideElementsPerProc = n/n_p;
	
	const int elementsPerProc = sideElementsPerProc*sideElementsPerProc;
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
	MPI_Cart_shift(Cycle_Communication,1, 1,&left,&right);
	MPI_Cart_shift(Cycle_Communication,0, -1,&down,&up);

		
	MPI_Datatype sideValues_t;
	MPI_Datatype sideValues;
	MPI_Type_vector(sideElementsPerProc, 1, sideElementsPerProc, MPI_DOUBLE, &sideValues);
	MPI_Type_create_resized(sideValues, 0, sizeof(double), &sideValues_t);
	MPI_Type_commit(&sideValues_t);

	double q0, q1, q1_sub,
	       DQ, DQ_sub, g_norm,
	       tau, beta;

	double q0_sub = 0;
	start = MPI_Wtime();

	/* INITIALIZATION OF b, g.
	   CALCULATION OF q_0 	   */
	for(int i = 0; i < elementsPerProc; i++){
		double x_i = h*(mycoords[1]*sideElementsPerProc+i%sideElementsPerProc);
		double y_i = h*(mycoords[0]*sideElementsPerProc+i/sideElementsPerProc);
		
		if((mycoords[0] == 0 && i < sideElementsPerProc) ||  (mycoords[0] == n_p-1 && i > elementsPerProc-sideElementsPerProc) || 
					(mycoords[1] == 0 && i % sideElementsPerProc == 0) || (mycoords[1] == n_p-1 && i % sideElementsPerProc == sideElementsPerProc-1)){
			local_d[i] = 0;
		}else{
			local_d[i] =  2*h*h*(x_i*(1-x_i)+y_i*(1-y_i));
		}
		
		local_g[i] = -local_d[i];
		
		q0_sub += local_g[i]*local_g[i];

	}

	MPI_Allreduce(&q0_sub, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	
	// Iteratation start here
	for(int t = 0; t < 1; t++){


		q1_sub = 0;
		DQ_sub = 0;

		/* DISTRIBUTION OF GHOSTROWS */
		if(mycoords[0] != n_p-1) MPI_Isend(&local_d[elementsPerProc-sideElementsPerProc], sideElementsPerProc, MPI_DOUBLE, down, 1, Cycle_Communication, &request);
		if(mycoords[0] != 0)     MPI_Isend(&local_d[0], sideElementsPerProc, MPI_DOUBLE, up, 0, Cycle_Communication, &request);
		if(mycoords[1] != n_p-1) MPI_Isend(&local_d[sideElementsPerProc-1], 1, sideValues_t, right, 3, Cycle_Communication, &request);
		if(mycoords[1] != 0)     MPI_Isend(&local_d[0], 1, sideValues_t, left, 2, Cycle_Communication, &request);
		

		if(p == 1){
			local_q[0] = 0;			
			local_q[sideElementsPerProc-1] = 0;
			local_q[elementsPerProc-sideElementsPerProc] = 0;
			local_q[elementsPerProc-1] = 0;
		
		}else if(mycoords[0] == 0){
			// TOP LEFT CORNER
			if(mycoords[1] == 0){
				local_q[0] = 0;			
				local_q[sideElementsPerProc-1] = 0;
				local_q[elementsPerProc-sideElementsPerProc] = 0;
				
			}

			// TOP RIGHT CORNER 
			else if (mycoords[1] == n_p-1){
				local_q[0] = 0;
				local_q[sideElementsPerProc-1] = 0;
				local_q[elementsPerProc-1] = 0;

			} else {
				// TOP SIDE 
				local_q[0] = 0;
				local_q[sideElementsPerProc-1] = 0;
			}
		} else if (mycoords[0] == n_p-1){
			// BOTTOM LEFT CORNER

			if(mycoords[1]==0){
				local_q[0] = 0;
				local_q[elementsPerProc-sideElementsPerProc] = 0;
				local_q[elementsPerProc-1] = 0;
			}
			// BOTTOM RIGHT CORNER
			else if (mycoords[1] == n_p-1){
				local_q[sideElementsPerProc-1] = 0;
				local_q[elementsPerProc-sideElementsPerProc] = 0;
				local_q[elementsPerProc-1] = 0;
			} else {
				// BOTTOM SIDE

				local_q[elementsPerProc-sideElementsPerProc] = 0;
				local_q[elementsPerProc-1] = 0;
			}
		} else if (mycoords[1] == 0){
			// LEFT SIDE
			local_q[0] = 0;
			local_q[elementsPerProc-sideElementsPerProc] = 0;

		} else if (mycoords[1] == n_p-1){
			// RIGHT SIDE
			local_q[sideElementsPerProc-1] = 0;	
			local_q[elementsPerProc-1] = 0;

		}

		
		/* CENTER ELEMENTS */
		
		for(int i = 1; i < sideElementsPerProc-1; i++){
			for(int j = 1; j < sideElementsPerProc-1;j++){
				int index = i*sideElementsPerProc + j;
				local_q[index] = -local_d[index+1] - local_d[index-1] + 4*local_d[index]-local_d[index+sideElementsPerProc] -local_d[index-sideElementsPerProc];
			}
		}

		/* INNER TOP ELEMENTS*/
		if(mycoords[0] == 0){
			for(int i = 0; i < sideElementsPerProc; i++)
				local_q[i] = 0;
		} else {
			MPI_Recv(topDest, sideElementsPerProc, MPI_DOUBLE, up, 1, Cycle_Communication, &status);
			for(int i = 1; i < sideElementsPerProc-1; i++){
				local_q[i] = -local_d[i +1] - local_d[i-1] + 4*local_d[i] - local_d[i+sideElementsPerProc] -topDest[i];
			}
			
		}


		/* INNER BOTTOM ELEMENTS */
		
		if(mycoords[0] == n_p-1){
			for(int i = elementsPerProc-sideElementsPerProc; i < elementsPerProc; i++)
				local_q[i] = 0;
		} else {
			MPI_Recv(bottomDest, sideElementsPerProc, MPI_DOUBLE, down, 0, Cycle_Communication, &status);
			for(int i = elementsPerProc-sideElementsPerProc+1; i < elementsPerProc-1; i++){
				local_q[i] = -local_d[i+1] - local_d[i-1] + 4*local_d[i] - local_d[i-sideElementsPerProc] -bottomDest[i-(elementsPerProc-sideElementsPerProc)];
			}
		}
		

		


		/* LEFT ELEMENTS */

		if(mycoords[1] == 0){
			for(int i = 0; i < sideElementsPerProc; i++){
				int index = i*sideElementsPerProc;		
				local_q[index] = 0;
			}
		} else {
			MPI_Recv(&leftDest[0], sideElementsPerProc, MPI_DOUBLE, left, 3, Cycle_Communication, &status);
			for(int i = 1; i < sideElementsPerProc-1; i++){
				int index = i*sideElementsPerProc;			
				local_q[index] = -local_d[index+1] - leftDest[i] + 4*local_d[index]-local_d[index-sideElementsPerProc]- local_d[index+sideElementsPerProc];
				if(i == sideElementsPerProc-2) printf("%.10lf\n", local_q[index]);
				//if(i != sideElementsPerProc-2 && mycoords[0] == n_p-1) local_q[index] -= local_d[index+sideElementsPerProc];
				/* TOP AND BOTTOM, EDGE CASES*/
				/*if(i != 0){
					local_q[index] -= local_d[index-sideElementsPerProc];
				} else {
					if(mycoords[0]==0){
						local_q[index] = 0;
						continue;
					} else {
						local_q[index] -= topDest[0];
					}
				}

				if(i != sideElementsPerProc-1){
					local_q[index] -= local_d[index+sideElementsPerProc];
				} else {
					if(mycoords[0]==n_p-1){
						local_q[index] = 0;
					} else {
						local_q[index] -= bottomDest[0];
					}
				}*/
			}
		}

		/* RIGHT ELEMENTS */
		if(mycoords[1] == n_p-1){
			for(int i = 0; i < sideElementsPerProc; i++){
				int index = (i+1)*sideElementsPerProc-1;
				local_q[index] = 0;
			}
		}else{
			MPI_Recv(&rightDest[0], sideElementsPerProc, MPI_DOUBLE, right, 2, Cycle_Communication, &status);
			for(int i = 1; i < sideElementsPerProc-1; i++){
				int index = (i+1)*sideElementsPerProc-1;

				local_q[index] = -local_d[index-1] - rightDest[i] + 4*local_d[index]-local_d[index-sideElementsPerProc]-local_d[index+sideElementsPerProc];	
				//if(i != sideElementsPerProc -2 && mycoords[0] == n_p-1) local_q[index] -= ;
				/* TOP AND BOTTOM, EDGE CASES*/
				/*if(i != 0){
					local_q[index] -= local_d[index-sideElementsPerProc];
				} else {
					if(mycoords[0]==0){
						local_q[index] = 0;
						continue;
					} else {
						local_q[index] -= topDest[sideElementsPerProc-1];
					}
				}

				if(i != sideElementsPerProc-1){
					local_q[index] -= local_d[index+sideElementsPerProc];
				} else {
					if(mycoords[0]==n_p-1){
						local_q[index] = 0;
					} else {
						local_q[index] -= bottomDest[sideElementsPerProc-1];
					}
				}*/
			}
		}
		if(p == 1){
			printf("All corners PE %i\n", rank);
			local_q[0] = 0;			
			local_q[sideElementsPerProc-1] = 0;
			local_q[elementsPerProc-sideElementsPerProc] = 0;
			local_q[elementsPerProc-1] = 0;
		
		}else if(mycoords[0] == 0){
			// TOP LEFT CORNER
			if(mycoords[1] == 0){
				printf("Top left corner PE %i\n", rank);
				local_q[0] = 0;			
				local_q[sideElementsPerProc-1] = 0;
				local_q[elementsPerProc-sideElementsPerProc] = 0;
				local_q[elementsPerProc-1] = -rightDest[sideElementsPerProc-1]+4*local_d[elementsPerProc-1]-bottomDest[sideElementsPerProc-1]-local_d[elementsPerProc-1-sideElementsPerProc]-local_d[elementsPerProc-2];
				
			}

			// TOP RIGHT CORNER 
			else if (mycoords[1] == n_p-1){
				printf("Top right corner PE %i\n", rank);
				local_q[0] = 0;
				local_q[sideElementsPerProc-1] = 0;
				local_q[elementsPerProc-sideElementsPerProc] = -local_d[elementsPerProc-2*sideElementsPerProc]+4*local_d[elementsPerProc-sideElementsPerProc]-local_d[elementsPerProc-sideElementsPerProc+1]-leftDest[sideElementsPerProc-1]-bottomDest[0];
				local_q[elementsPerProc-1] = 0;

			} else {
				// TOP SIDE 
				printf("Top side PE %i\n", rank);
				local_q[0] = 0;
				local_q[sideElementsPerProc-1] = 0;
				local_q[elementsPerProc-sideElementsPerProc] = -local_d[elementsPerProc-2*sideElementsPerProc]+4*local_d[elementsPerProc-sideElementsPerProc]-local_d[elementsPerProc-sideElementsPerProc+1]-leftDest[sideElementsPerProc-1]-bottomDest[0];
				local_q[elementsPerProc-1] = -rightDest[sideElementsPerProc-1]+4*local_d[elementsPerProc-1]-bottomDest[sideElementsPerProc-1]-local_d[elementsPerProc-1-sideElementsPerProc]-local_d[elementsPerProc-2];
			}
		} else if (mycoords[0] == n_p-1){
			// BOTTOM LEFT CORNER

			if(mycoords[1]==0){
				printf("Bottom left corner PE %i\n", rank);
				local_q[0] = 0;
				local_q[sideElementsPerProc-1] = -rightDest[0]+4*local_d[sideElementsPerProc-1]-topDest[sideElementsPerProc-1]-local_d[2*sideElementsPerProc-1]-local_d[sideElementsPerProc-2];
				local_q[elementsPerProc-sideElementsPerProc] = 0;
				local_q[elementsPerProc-1] = 0;
			}
			// BOTTOM RIGHT CORNER
			else if (mycoords[1] == n_p-1){
				printf("Bottom right corner PE %i\n", rank);
				local_q[0] = -local_d[sideElementsPerProc]+4*local_d[0]-local_d[1]-leftDest[0]-topDest[0];
				local_q[sideElementsPerProc-1] = 0;
				local_q[elementsPerProc-sideElementsPerProc] = 0;
				local_q[elementsPerProc-1] = 0;
			} else {
				// BOTTOM SIDE
				printf("Bottom side PE %i\n", rank);
				local_q[0] = -local_d[sideElementsPerProc]+4*local_d[0]-local_d[1]-leftDest[0]-topDest[0];
				local_q[sideElementsPerProc-1] = -rightDest[0]+4*local_d[sideElementsPerProc-1]-topDest[sideElementsPerProc-1]-local_d[2*sideElementsPerProc-1]-local_d[sideElementsPerProc-2];
				local_q[elementsPerProc-sideElementsPerProc] = 0;
				local_q[elementsPerProc-1] = 0;
			}
		} else if (mycoords[1] == 0){
			// LEFT SIDE
			printf("Left side PE %i\n", rank);
			local_q[0] = 0;
			local_q[sideElementsPerProc-1] = -rightDest[0]+4*local_d[sideElementsPerProc-1]-topDest[sideElementsPerProc-1]-local_d[2*sideElementsPerProc-1]-local_d[sideElementsPerProc-2];
			local_q[elementsPerProc-sideElementsPerProc] = 0;
			local_q[elementsPerProc-1] = -rightDest[sideElementsPerProc-1]+4*local_d[elementsPerProc-1]-bottomDest[sideElementsPerProc-1]-local_d[elementsPerProc-1-sideElementsPerProc]-local_d[elementsPerProc-2];
		} else if (mycoords[1] == n_p-1){
			// RIGHT SIDE
			printf("Right side PE %i\n", rank);
			local_q[0] = -local_d[sideElementsPerProc]+4*local_d[0]-local_d[1]-leftDest[0]-topDest[0];
			local_q[sideElementsPerProc-1] = 0;
			local_q[elementsPerProc-sideElementsPerProc] = -local_d[elementsPerProc-2*sideElementsPerProc]+4*local_d[elementsPerProc-sideElementsPerProc]-local_d[elementsPerProc-sideElementsPerProc+1]-leftDest[sideElementsPerProc-1]-bottomDest[0];
			local_q[elementsPerProc-1] = 0;

		}else {
			// NO EDGES
			printf("No edges PE %i\n", rank);


			// res = under + center - höger - vänster - toppen
			local_q[0] = -local_d[sideElementsPerProc]+4*local_d[0]-local_d[1]-leftDest[0]-topDest[0];
			//res = - Höger+ center - toppen - under - vänster
			local_q[sideElementsPerProc-1] = -rightDest[0]+4*local_d[sideElementsPerProc-1]-topDest[sideElementsPerProc-1]-local_d[2*sideElementsPerProc-1]-local_d[sideElementsPerProc-2];
			
			//res = -Toppen + center - höger - vänster -botten
			local_q[elementsPerProc-sideElementsPerProc] = -local_d[elementsPerProc-2*sideElementsPerProc]+4*local_d[elementsPerProc-sideElementsPerProc]-local_d[elementsPerProc-sideElementsPerProc+1]-leftDest[sideElementsPerProc-1]-bottomDest[0];
			

			// res = höger + Center - botten - Toppen - vänster
			local_q[elementsPerProc-1] = -rightDest[sideElementsPerProc-1]+4*local_d[elementsPerProc-1]-bottomDest[sideElementsPerProc-1]-local_d[elementsPerProc-1-sideElementsPerProc]-local_d[elementsPerProc-2];
		}

		
		for(int i = 0; i < elementsPerProc; i++){
			DQ_sub += local_q[i]*local_d[i];

		}
		MPI_Allreduce(&DQ_sub, &DQ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tau = q0/DQ;		
		
		for(int i = 0; i < elementsPerProc; i++){
			local_u[i] += tau*local_d[i];
			local_g[i] += tau*local_q[i];
			q1_sub += local_g[i]*local_g[i];
			
		}


		MPI_Allreduce(&q1_sub, &q1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		beta = q1/q0;
		
		for(int i = 0; i < elementsPerProc; i++){
			//local_d[i] = -local_g[i]+beta*local_d[i];
		}
		
		q0 = q1;

		g_norm = sqrt(q1);
		//if(rank==0)printf("g_norm %lf\n", g_norm);
	}

	end_time = MPI_Wtime()-start;
	double max_time;
	MPI_Reduce(&end_time, &max_time,1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank==0) printf("%lf %lf\n",g_norm, end_time);

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
	MPI_Gatherv(local_d, elementsPerProc, MPI_DOUBLE, res, count, displ, blocktype, 0, MPI_COMM_WORLD);
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




