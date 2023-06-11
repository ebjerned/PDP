#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//Readwrite
int read_input(const char *file_name, int **values);
int write_output(char *file_name, const int *output, int num_values);


//Comparison function for c.qsort
int compare (const void * a, const void * b);
//For pivot strat
int findMedian(int *, int);
int find_mean(int *data, int len);
int pivot(int pivot_strat,int len, int *rcv_buffer, int rank, int size, MPI_Comm group);
int find_pivot(int* data, int len);
//Recusive funciton
int* Parallel_Qsort(MPI_Comm* groups, int rank, int size,    int* local_arr, int chunk_size, int depth);
//Merge two arrays
int* merge( int *list1, int n1,    int* list2, int n2, int rank);


void splitGroup(MPI_Comm curr_group, int rank, int size, int n, MPI_Comm* groups);

//Global Vars
int STRAT;
MPI_Status status;
int padding;
int currentsize = 0;


int main(int argc, char **argv){

	if (4 != argc) {
		printf("Usage: qsort input_file output_file pivotstrategy \n");
		return 1;
	}

	char *input_name = argv[1];
	char *output_name = argv[2];
	STRAT = atoi(argv[3]);    
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

	int SIZE;
	int n;

	int *input;
	int *padded_input;
	int paddedSize;
	if (rank == 0) {

		if (0 > (SIZE = read_input(input_name, &input))) {
			return 2;
		}
		paddedSize =(int)(ceil(SIZE/(double)size)*size);
		padding = paddedSize -SIZE;
		padded_input = (int*) calloc(paddedSize, sizeof(int));
		// Copy the original array into the padded array to make sure that its length is divisble by size
		for (int i = 0; i < SIZE; i++) {
			padded_input[i] = input[i];
		}

	}
	MPI_Bcast(&paddedSize, 1, MPI_INT, 0,MPI_COMM_WORLD);
	n = paddedSize / size;
	//MPI_Barrier(MPI_COMM_WORLD);    
	int *local_arr = (int *)malloc(n * sizeof(int));
	MPI_Scatter(padded_input, n, MPI_INT, local_arr, n, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank==0) free(padded_input);
	double start, end;
	//Sort locally before anything else
	start = MPI_Wtime();
	qsort(local_arr, n , sizeof(int), compare);
	//if(size==1) end = MPI_Wtime()-start;
    
    	MPI_Comm groups[(int)log2(size)];
	splitGroup(MPI_COMM_WORLD, rank, size, 0, &groups);
    
	// Running and then gathering data;
	int* finaldata;
	if (size > 1){
		//start = MPI_Wtime();
		local_arr = Parallel_Qsort(groups,rank, size, local_arr,n, 0);

		//MPI_Barrier(MPI_COMM_WORLD);
		end = MPI_Wtime() - start;
		double tmp;
		MPI_Reduce(&end, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		end = tmp;

		if(rank != 0){
			MPI_Send(local_arr, currentsize, MPI_INT, 0, rank, MPI_COMM_WORLD);
		}else{
			// Zero process gathers data.
			int currentlocation = 0;
			finaldata=(int*) malloc(paddedSize*sizeof(int));
			memcpy(finaldata + currentlocation, local_arr, currentsize*sizeof(int));
			currentlocation = currentlocation + currentsize;
			int* tmp;

			for(int r = 1; r < size; ++r){
				int recive_size = 0;
				MPI_Probe(r, r, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_INT, &recive_size);

				tmp = (int*)malloc(recive_size*sizeof(int));
				MPI_Recv(tmp, recive_size, MPI_INT, r, r, MPI_COMM_WORLD, NULL);

				memcpy(finaldata+currentlocation, tmp, recive_size*sizeof(int));               
				currentlocation = currentlocation+recive_size;
					    
				free(tmp);
			}
		}
	} else{
		end = MPI_Wtime()-start;	
	}

	if (rank==0){    
		int *finald = finaldata;
		if(size != 1){
			write_output(output_name, finald+padding, SIZE);
		} else {
			write_output(output_name, &local_arr[0], SIZE);
		}
		printf("%lf\n", end);
		free(input);
		free(finaldata);
	}

	free(local_arr);
        

	MPI_Finalize();
}

int* Parallel_Qsort(MPI_Comm* groups, int rank, int size,  int* local_arr, int n, int depth){
	MPI_Comm curr_group = groups[depth];
	MPI_Comm_size( curr_group, &size);
	MPI_Comm_rank( curr_group, &rank );

	//Recursion statement 
	if(size < 2){ 
		return(local_arr);
	}

	int PivotPoint = pivot(STRAT, n, local_arr, rank, size, curr_group);

	//Find the size of arrays to share according to 

	int len_top = 0;
	int len_bot = 0;
	for(int i = 0 ; i < n;i++){
		if(local_arr[i] >= PivotPoint){
			len_top = n-len_bot;
			break;
		}else{
			len_bot++;
		}
	}


	//Group up process ids into two groups, high and low, every process has a pair-process in other group

	int group_size = size/2;

	int pair_high;
	int pair_low;

	int recv_size_low = 0;
	int recv_size_high = 0; 
	if(rank < group_size){
		pair_high = rank + group_size;
		recv_size_low = 0;

	}else{
		pair_low = rank - group_size;
		recv_size_high = 0;
	}
	int tag = size;

	//Sending and recieving

	int *recv_low_from_high;
	int *recv_high_from_low;
	MPI_Request request1;
	MPI_Request request2;

	int *merged_arr;
	int newlen;

	MPI_Comm nextgroup;
	int newsize,newrank;


	if (rank < group_size) {
		//low process sends its top to pair in right group and receive pairs bottom
		MPI_Isend(&local_arr[len_bot], len_top, MPI_INT, pair_high, tag+1, curr_group, &request1);
		MPI_Probe(pair_high, tag, curr_group, &status);
		MPI_Get_count(&status, MPI_INT, &recv_size_low);

		recv_low_from_high = (int*) realloc(recv_low_from_high, recv_size_low*sizeof(int));
		MPI_Recv(recv_low_from_high, recv_size_low, MPI_INT, pair_high, tag, curr_group, &status);

		merged_arr = merge(local_arr,len_bot,recv_low_from_high,recv_size_low, rank);
		newlen = len_bot + recv_size_low;
		currentsize = newlen;

		//free(recv_low_from_high);

		// Split group
		/*MPI_Comm_split( curr_group,0, 1, &nextgroup );*/
		//MPI_Comm_size( groups[depth+1], &newsize );
		//MPI_Comm_rank( groups[depth+1], &newrank );
		return Parallel_Qsort(groups, newrank, group_size/2, merged_arr, newlen, depth+1);

	} else {
		//high process sends its bottom to pair in left group and receive pairs top
		MPI_Isend(local_arr, len_bot, MPI_INT, pair_low, tag, curr_group, &request1);
		MPI_Probe(pair_low, tag+1, curr_group, &status);
		MPI_Get_count(&status, MPI_INT, &recv_size_high);

		recv_high_from_low = (int*) realloc(recv_high_from_low, recv_size_high*sizeof(int));
		MPI_Recv(recv_high_from_low, recv_size_high, MPI_INT, pair_low, tag+1, curr_group, &status);

		merged_arr = merge(&local_arr[len_bot],len_top,recv_high_from_low,recv_size_high, rank);
		newlen = len_top + recv_size_high;
		currentsize = newlen;

		//free(recv_high_from_low);

		// Split group
		/*MPI_Comm_split( curr_group, 1, 1, &nextgroup );*/
		//MPI_Comm_size( groups[depth+1], &newsize );
		//MPI_Comm_rank( groups[depth+1], &newrank );
		return Parallel_Qsort(groups, newrank, group_size/2, merged_arr, newlen, depth+1);
	}
	//MPI_Barrier(curr_group);

	free(recv_low_from_high);
	free(recv_high_from_low);
	return(local_arr);
}

void splitGroup(MPI_Comm curr_group, int rank, int size, int n, MPI_Comm* groups){

	MPI_Comm_size( curr_group, &size);
	MPI_Comm_rank( curr_group, &rank );

	//Recursion statement 
	if(size < 2){
		groups[n] = curr_group;
		return;
	}

	//Group up process ids into two groups, high and low, every process has a pair-process in other group

	int group_size = size/2;
	//Sending and recieving

	MPI_Comm nextgroup;
	int newsize,newrank;

	groups[n] = curr_group;
	if (rank < group_size) {

		MPI_Comm_split( curr_group,0, 1, &nextgroup );
		MPI_Comm_size( nextgroup, &newsize );
		MPI_Comm_rank( nextgroup, &newrank );
		return splitGroup(nextgroup, newrank, newsize, n+1, groups);

	} else {

		MPI_Comm_split( curr_group, 1, 1, &nextgroup );
		MPI_Comm_size( nextgroup, &newsize );
		MPI_Comm_rank( nextgroup, &newrank );
		return splitGroup(nextgroup, newrank, newsize, n+1, groups);
	}

}


int pivot(int pivot_strat,int len, int *local_arr, int rank, int size, MPI_Comm group){
	int l_pivot = 0;
	MPI_Comm_size( group, &size);
	MPI_Comm_rank( group, &rank );
	switch (pivot_strat){
		case 1 : /* 1. Select the median in one processor in each group of processors*/
			{
			/* ROOT find the Global PIVOT and broadcast it to all other PROCS */
				if (rank==0){
					l_pivot = find_pivot(local_arr, len);
				}

				MPI_Bcast(&l_pivot, 1, MPI_INT, 0, group);

				return l_pivot;
			}
			break;
		case 2 : /* Select the median of all medians in each processor group */
			{
				int local_pivot = find_pivot(local_arr, len); /* Each PROC find median of local buffer */
				int all_median[size]; /* # of median values as our #PROCS */

				MPI_Allgather(&local_pivot, 1, MPI_INT, all_median, 1, MPI_INT, group);
				l_pivot = find_pivot(all_median, size); /* median of all medians */

				return l_pivot;
			}
			break;
		case 3: /* Select the mean value of all medians in each processor group */
			{
				int local_pivot = find_pivot(local_arr, len); /* Each PROC find median of local buffer */
				int all_median[size]; /* # of median values as our #PROCS */
				MPI_Allgather(&local_pivot, 1, MPI_INT, all_median, 1, MPI_INT, group);
				l_pivot = find_mean(all_median, size);
				return l_pivot;
			}
			break;
	}
	return 0;
}

int compare(const void * a, const void * b){
        return (*(int*)a) - (*(int*)b);
}

int find_mean(int *data, int len){
        int sum = 0;
        for (int i=0; i<len; ++i)
                sum += data[i];
        sum = sum/(int)len;
        return sum;
}

int find_pivot(int* data, int len){
        /* If the number of elements is zero the pivot is data[0] */
        if (len==0)
                return 0;
        int pivot = 0;
        /* If number of elements is odd, pivot is the middle element */
        if (len % 2 != 0)
                pivot = data[len/2];
        /* Else it is the average of the two middle elements */
        else
                pivot = ((int)data[(len / 2) -1] + (int)data[(len / 2)]) / 2.0;
        return pivot;
}


int* merge( int *v1, int n1,    int *v2, int n2 , int rank){

	int *result = (int*)malloc((n1+n2)*sizeof(int));

	for (int k = 0, i=0, j=0; k < n1+n2;k++) {
		if(i > n1-1){ 
			result[k] = v2[j++];
		}else if(j > n2-1){
			result[k] = v1[i++];
		}else{
			if (v1[i] < v2[j]){
				result[k] = v1[i++];
			}else{
				result[k] = v2[j++];
			}
		}
	}
	return result;
}
int read_input(const char *file_name, int **values) {
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

	if (NULL == (*values = malloc((num_values) * sizeof(int)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
 
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%i", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}

int write_output(char *file_name, const int *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%i ", output[i])) {
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

