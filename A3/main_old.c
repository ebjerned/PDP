
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, double **values);
int write_output(char *file_name, const double *output, int num_values);

int comparisonFunction(const void* a, const void* b);
int compare (const void * a, const void * b);
int median(int* values, int n_values);
int findPivot(int* pivots, int processes, char method);
double findMedian(int *, int);
void partition(int *, int, int, int *, int *);
double* Parallel_Qsort(MPI_Comm curr_group, int rank, int size, double* local_arr, int chunk_size);
double pivot(int pivot_strat,int len, double *rcv_buffer, int rank, int size, MPI_Comm group);
double* merge(double *list1, int n1, double* list2, int n2 );
int median(int* values, int n_values);
double find_pivot(double* data, int len);
//Global vars
int STRAT;
int currentsize;


int main(int argc, char **argv) {


	if (4 != argc) {
		printf("Usage: qsort input_file output_file pivotstrategy \n");
		return 1;
	}
  char *input_name = argv[1];
	char *output_name = argv[2];
  STRAT = atoi(argv[3]);
  MPI_Status status;

  
  MPI_Init(&argc, &argv);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  int SIZE;
  int n;
  int* pivots;
  int pivotFromStrategy;
  double *input;
  if (rank == 0) {
	  if (0 > (SIZE = read_input(input_name, &input))) {
		  return 2;
	  }
    // Computing chunk size
    n = SIZE /size;
   }
  // BroadCast the Size to all the
  // process from root process
  
  MPI_Bcast(&SIZE, 1, MPI_INT, 0,MPI_COMM_WORLD);
  n = SIZE /size; 
  MPI_Barrier(MPI_COMM_WORLD);  
  double *local_arr = (double *)malloc(n * sizeof(double));
  MPI_Scatter(input, n, MPI_DOUBLE, local_arr, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  qsort(local_arr, n , sizeof(double), compare );
  

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm current_comm = MPI_COMM_WORLD;
  
  double *rcv_buffer = (double*)malloc(n*sizeof(double));
  
  rcv_buffer = Parallel_Qsort(MPI_COMM_WORLD,rank, size, local_arr,n);
  
  for(int i = 0; i< n; i++){
    printf("val %lf \n", rcv_buffer[i] );  
  }

  
  MPI_Finalize();

  return 0;
}


double* Parallel_Qsort(MPI_Comm curr_group, int rank, int size, double* local_arr, int chunk_size){
  MPI_Comm_size( curr_group, &size);
  MPI_Comm_rank( curr_group, &rank );
  //recursion
  if (size < 2){
  	return (local_arr);
  }

  //Find Pivot
  double PivotPoint = pivot(STRAT, chunk_size, local_arr, rank, size, curr_group);
  //Find size and split according to pivot
  double *low=(double *)malloc(chunk_size*sizeof(double));
  double *high=(double *)malloc(chunk_size*sizeof(double));
  int l_low ,h_high;
  l_low = 0;
  h_high = 0;
  for (int i = 0; i < chunk_size; i++){
  	if (local_arr[i] <= PivotPoint){
  		low[l_low] = local_arr[i];
  		l_low = l_low+1;
  	}else{
  		high[h_high] = local_arr[i];
   		h_high = h_high+1;
   	}
  }
  printf("Pivot = %lf \n", PivotPoint);
  printf("h_high = %d, l_low = %d \n", h_high,l_low);
  //Group up processes into two group, high and low, every process has a pair-process in other group
  
  int group_size = size /2;
  
  int pair_high;
  int pair_low;
  
  if(rank < group_size){
    pair_high = rank + group_size;
    
  }else{
    pair_low = rank - group_size;
  }
  
  int tag = size;

 
  MPI_Request request1;
  MPI_Request request2;
  MPI_Status status1;
  MPI_Status status2;
  double *otherhigh;
  double *otherlow;
  int recv_size_low = 0;
  int recv_size_high = 0;
  if (rank < group_size){
    	//low process sends its highdata to pairhigh and receive pairhigh's lowdata from pairhigh
    	
  	MPI_Isend(high, h_high, MPI_DOUBLE, pair_high, tag+1, curr_group,&request1);
   
   	MPI_Probe(pair_high, tag, curr_group, &status1);
    
   	MPI_Get_count(&status1, MPI_DOUBLE, &recv_size_low);
     
   	otherlow=(double *)malloc(recv_size_low*sizeof(double));

   	MPI_Irecv(otherlow, recv_size_low, MPI_DOUBLE, pair_high, tag, curr_group, &request2);
    
    printf("Recv_done\n"); 
    printf("recv_size_low = %d \n" , recv_size_low);
    
   	MPI_Wait(&request1, &status1);
    
   	MPI_Wait(&request2, &status2);

  }else{
   	//high process sends its lowdata to pairlow and receive pairlow's highdata from pairlow
   	MPI_Isend(low, l_low, MPI_DOUBLE, pair_low, tag, curr_group,&request1);
    
   	MPI_Probe(pair_low, tag+1, curr_group, &status2);
    
   	MPI_Get_count(&status2, MPI_DOUBLE, &recv_size_high);
    //allocate for recv
   	otherhigh=(double *)malloc(recv_size_high*sizeof(double));

   	MPI_Irecv(otherhigh, recv_size_high, MPI_DOUBLE, pair_low, tag+1, curr_group, &request2);
    printf("Recv_done\n"); 
    
    
   	MPI_Wait(&request1, &status1);
    
   	MPI_Wait(&request2, &status2);	
    
  }

  //Merge data into one array
  double *merged_arr;
  int newlen;
  free(local_arr);
  MPI_Comm nextgroup;
  int newsize,newrank;
  if (rank < group_size){
  
    merged_arr = merge(low,l_low,otherlow,recv_size_low);
    newlen = l_low + recv_size_low;
    currentsize = newlen;


    free(low);
    free(otherlow);
    	// Split group into two. This is color 0.
    MPI_Comm_split( curr_group, 0, 1, &nextgroup );
    MPI_Comm_size( nextgroup, &newsize );
    MPI_Comm_rank( nextgroup, &newrank );
 	  //Recursive call
 	  return (Parallel_Qsort(nextgroup, newrank, newsize, merged_arr, newlen));
   }else{
  	merged_arr = merge(high,h_high,otherhigh,recv_size_high);

    
   	newlen = h_high + recv_size_high;
   	currentsize = newlen;
    
   	//qsort(merged_arr, newlen, sizeof(double), compare);

    free(high);
    free(otherhigh);
    // Split group into two. This is color 1.
    MPI_Comm_split( curr_group, 1, 1, &nextgroup );
    MPI_Comm_size( nextgroup, &newsize );
    MPI_Comm_rank( nextgroup, &newrank );
    //Recursive call
    return (Parallel_Qsort(nextgroup, newrank, newsize, merged_arr, newlen));
    }
}


double* merge(double *v1, int n1, double *v2, int n2 ){
    int i;
    int j = 0;
    int k = 0;
    double * result = (double*)malloc((n1+n2)*sizeof(double));
    
    i=0; j=0; k=0;
    while(i<n1 && j<n2)
        if(v1[i]<v2[j])
        {
            result[k] = v1[i];
            i++; k++;
        }
        else
        {
            result[k] = v2[j];
            j++; k++;
        }
    if(i==n1)
        while(j<n2)
        {
            result[k] = v2[j];
            j++; k++;
        }
    else
        while(i<n1)
        {
            result[k] = v1[i];
            i++; k++;
        }
    return result;
}


int comparisonFunction(const void* a, const void* b){
	return (*(double*)a - *(double*)b);
}

int compare(const void * a, const void * b)
{
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;
}

int median(int* values, int n_values){
	if(n_values % 2 == 0)
		return (values[n_values/2] + values[n_values/2+1])/2;

	return values[n_values/2];
}

int mean(int* values, int n_values){
	int sum = 0;
	for(int i = 0; i < n_values; i++){
		sum += values[i];
	}
	return sum / n_values;
}

double find_pivot(double* data, int len){
    /* If the number of elements is zero the pivot is data[0] */
    if (len==0)
        return 0;
    double pivot = 0;
    /* If number of elements is odd, pivot is the middle element */
    if (len % 2 != 0)
        pivot = data[len/2];
    /* Else it is the average of the two middle elements */
    else
        pivot =  ((double)data[(len / 2) -1] + (double)data[(len / 2)]) / 2.0;
    return pivot;
}

double find_mean(double *data, int len){
    double sum = 0;
    for (int i=0; i<len; ++i)
        sum += data[i];
    sum = sum/(double)len;
    return sum;
}

double pivot(int pivot_strat,int len, double *rcv_buffer, int rank, int size, MPI_Comm group){
    double l_pivot = 0;
    MPI_Comm_size( group, &size);
    MPI_Comm_rank( group, &rank );
    switch (pivot_strat){
        case 1 : /* 1. Select the median in one processor in each group of processors*/
        {
            /* ROOT find the Global PIVOT and broadcast it to all other PROCS */
            if (rank==0){
                l_pivot = find_pivot(rcv_buffer, len);


            }
       
            MPI_Bcast ( &l_pivot, 1, MPI_DOUBLE, 0, group);
            
	    return l_pivot;
        }
        break;
        case 2 : /* Select the median of all medians in each processor group */
        {
            double local_pivot = find_pivot(rcv_buffer, len); /* Each PROC find median of local buffer */
            double all_median[size]; /* # of median values as our #PROCS */
            
            MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, group);
            l_pivot = find_pivot(all_median, size); /* median of all medians */
   
	    return l_pivot;
        }
        break;
        case 3: /* Select the mean value of all medians in each processor group */
        {
            double local_pivot = find_pivot(rcv_buffer, len); /* Each PROC find median of local buffer */
            double all_median[size]; /* # of median values as our #PROCS */
            MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, group);
            l_pivot = find_mean( all_median, size);
	    return l_pivot;
        }
        break;
    }
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
		if (0 > fprintf(file, "%lf ", output[i])) {
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

