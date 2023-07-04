#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

const int MAX_TIME = 100; 

//Propensity calculation function
void prop(int *x, double *w);

//SSA Function
int* ssa(int *);



int main(int argc, char **argv){

	if (2 != argc) {
		printf("Usage: malsim N \n");
		return 1;
	}

	int N = atoi(argv[1]); //Total number of simulation runs
    //char *file_name = argv[2];

	
 	MPI_Init(&argc, &argv);
 	int size, rank;
 	MPI_Comm_size(MPI_COMM_WORLD, &size);
 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	if (N%size != 0) { //Splitting the total runs to processing elements
		printf("N must be divisible by the number of processing elements\n");
		return 1;
	}

	int n = N/size;
	
	time_t rand_seed;
	srand((unsigned) time(&rand_seed));

	double start_time;
	double time;
	
	int* recbuf;	//Initializing recievebuffer
	int* x_out;		//Temporary "final x", statevector to use in the loop
	int* X;			//Final "X" with the number of susceptible humans from each run


	int* x_0 = (int*)malloc(7 * sizeof(int));
	x_0[0] = 900;
	x_0[1] = 900;
	x_0[2] = 30;
	x_0[3] = 330;
	x_0[4] = 50;
	x_0[5] = 270;
	x_0[6] = 20;


	if (rank == 0){
		recbuf = (int*)malloc(N * sizeof(int));
		if(NULL == recbuf) {
			perror("Couldn't allocate memory for output");
			return 2;
		}
	}

	X = (int*)malloc(n * sizeof(int));

	MPI_Barrier(MPI_COMM_WORLD);
	
	start_time = MPI_Wtime(); 


	for(int i = 0; i<n; i++){

		x_out = ssa(x_0);

		X[i] = x_out[0];	//Saving final number of susceptible humans to the vector that is to be gathered
		
	}


	MPI_Barrier(MPI_COMM_WORLD);

	time = MPI_Wtime() - start_time;
	double tmp;
	MPI_Reduce(&time, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank==0) printf("%lf\n", tmp);
	MPI_Gather(X, n, MPI_INT, recbuf, n, MPI_INT, 0, MPI_COMM_WORLD);		//Gathering all X vectors into the recieve buffer

	MPI_Barrier(MPI_COMM_WORLD);

	
	//if(rank == 0){
	//    FILE *fp = fopen(file_name, "w"); 
    //	if (fp == NULL) { 
    //	    printf("Failed to open file."); 
    //	    return 1; 
    //	} 
	//	fprintf(fp, "Processing took %fs\n", time);
    //	for (int i = 0; i < N; i++) { 
    //	    fprintf(fp, "%d ", recbuf[i]); 
    //	}
	//
    //fclose(fp); 
	//}



	free(recbuf);
	free(X);
	free(x_out);
	MPI_Finalize();
}

void prop(int *x, double *w) {
	// Birth number, humans
	const double LAMBDA_H = 20;
	// Birth number, mosquitoes
	const double LAMBDA_M = 0.5;
	// Biting rate of mosquitoes
	const double B = 0.075;
	/* Probability that a bite by an infectious mosquito results in transmission
	   of disease to human*/
	const double BETA_H = 0.3;
	/* Probability that a bite results in transmission of parasite to a
	   susceptible mosquito*/
	const double BETA_M = 0.5;
	// Human mortality rate
	const double MU_H = 0.015;
	// Mosquito mortality rate
	const double MU_M = 0.02;
	// Disease induced death rate, humans
	const double DELTA_H = 0.05;
	// Disease induced death rate, mosquitoes
	const double DELTA_M = 0.15;
	// Rate of progression from exposed to infectious state, humans
	const double ALFA_H = 0.6;
	// Rate of progression from exposed to infectious state, mosquitoes
	const double ALFA_M = 0.6;
	// Recovery rate, humans
	const double R = 0.05;
	// Loss of immunity rate, humans
	const double OMEGA = 0.02;
	/* Proportion of an antibody produced by human in response to the incidence
	   of infection caused by mosquito. */
	const double NU_H = 0.5;
	/* Proportion of an antibody produced by mosquito in response to the
	   incidence of infection caused by human. */
	const double NU_M = 0.15;

	w[0] = LAMBDA_H;
	w[1] = MU_H * x[0];
	w[2] = (B * BETA_H * x[0] * x[5]) / (1 + NU_H * x[5]);
	w[3] = LAMBDA_M;
	w[4] = MU_M * x[1];
	w[5] = (B * BETA_M * x[1]*x[4]) / (1 + NU_M * x[4]);
	w[6] = MU_H * x[2];
	w[7] = ALFA_H * x[2];
	w[8] = MU_M * x[3];
	w[9] = ALFA_M * x[3];
	w[10] = (MU_H + DELTA_H) * x[4];
	w[11] = R * x[4];
	w[12] = (MU_M + DELTA_M) * x[5];
	w[13] = OMEGA * x[6];
	w[14] = MU_H * x[6];
}

int* ssa(int *x){		//Stochastic simulation algorithm, in accordance with the project instructions

    int P[15][7] = 
    {{1,0,0,0,0,0,0},
	{-1,0,0,0,0,0,0},
    {-1,0,1,0,0,0,0},
    {0, 1,0,0,0,0,0},
    {0,-1,0,0,0,0,0},
    {0,-1,0,1,0,0,0},
    {0,0,-1,0,0,0,0},
    {0,0,-1,0,1,0,0},
    {0,0,0,-1,0,0,0},
    {0,0,0,-1,0,1,0},
    {0,0,0,0,-1,0,0},
    {0,0,0,0,-1,0,1},
    {0,0,0,0,0,-1,0},
    {1,0,0,0,0,0,-1},
    {0,0,0,0,0,0,-1}};

	double w[15] = {};

    double tau, a_0, b_0, b_1, u_1, u_2;
	int count = 0;
	double t = 0;

	while(t<MAX_TIME){
	    prop(x, w);
	    a_0 = 0, b_0 = 0, b_1 = 0, tau = 0;

	    for (int i = 0; i < 15; i++) {
	    	a_0 += w[i];
	    }

		u_1 = (double)rand()/(double)RAND_MAX;
	    u_2 = (double)rand()/(double)RAND_MAX;

	   	tau = -log(u_1) / a_0;

	   	int r = 0;
	   	while (u_2 > 0) {
	   		u_2 -= w[r] / a_0;
	   		r++;
	   	}
	   	r--;

		for(int i=0; i<7; i++){
	    	x[i] += P[r][i];
	    }
		t+=tau;
	   	count++;
	}
	return x;
}
