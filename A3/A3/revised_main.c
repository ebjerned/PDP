#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//Readwrite
int read_input(const char *file_name,/*int rank, int size, int* num_values,*/ double **values);
int write_output(char *file_name, const double *output, int num_values);


//Comparison function for c.qsort
int compare (const void * a, const void * b);
//For pivot strat
double findMedian(int *, int);
double find_mean(double *data, int len);
double pivot(int pivot_strat,int len, double *rcv_buffer, int rank, int size, MPI_Comm group);
double find_pivot(double* data, int len);
//Recusive funciton
double* Parallel_Qsort(MPI_Comm curr_group, int rank, int size, double* local_arr, int chunk_size);
//Merge two arrays
double* merge(double *list1, int n1, double* list2, int n2, double** mergeTo);

double* pad_array(double* arr, int len, int p);
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

        double *input;
        double *padded_input;

	
        if (rank == 0) {
                int addage; 
	        if (0 > (SIZE = read_input(input_name, &input))) {
		        return 2;
	        }
		printf("Size before: %i\n", SIZE);
                padding = SIZE - (SIZE % size); 
                SIZE = SIZE + (SIZE%size);//padding;
		printf("Size after: %i\n", SIZE);
                padded_input = (double*) calloc((SIZE) , sizeof(double));
                // Copy the original array into the padded array to make sure that its length is divisble by size
                for (int i = 0; i < SIZE; i++) {
                        padded_input[i] = input[i];
                }
                //memcpy(padded_input, input, (SIZE)*sizeof(double));
                n = SIZE/size;
 
        }
        MPI_Bcast(&SIZE, 1, MPI_INT, 0,MPI_COMM_WORLD);
        n = SIZE /size;
        MPI_Barrier(MPI_COMM_WORLD);        
        double *local_arr = (double *)malloc(n * sizeof(double));
        MPI_Scatter(padded_input, n, MPI_DOUBLE, local_arr, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i = 0; i<n; i++)
		printf("%i %lf\n", rank, local_arr[i]);
        double start , end;
        //Sort locally before anything else
        if(size==1) start = MPI_Wtime();
        qsort(local_arr, n , sizeof(double), compare );
        if(size==1) end = MPI_Wtime()-start;
	MPI_Barrier(MPI_COMM_WORLD);
        
        
        
        
        // Running and then gathering data;
        double* finaldata;
        if (size > 1){
                start = MPI_Wtime();
                local_arr = Parallel_Qsort(MPI_COMM_WORLD,rank, size, local_arr,n);
                MPI_Barrier(MPI_COMM_WORLD);
                end = MPI_Wtime() - start;

                
                if( rank != 0 ){
                        MPI_Send( local_arr, currentsize, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD );
			for(int i = 0; i < currentsize; i++)
				printf("Recieved %lf \n", local_arr[i]);
                }else{
                // Zero process gathers data.
                int currentlocation = 0;
                finaldata=(double *)malloc(SIZE*sizeof(double));
                memcpy(finaldata+currentlocation,local_arr,currentsize*sizeof(double));
                currentlocation = currentlocation+currentsize;
                double* tmp;

                	for( int r = 1; r < size; ++r ){
		                int recive_size = 0;
		                
		                MPI_Probe(r, r, MPI_COMM_WORLD, &status);
		                
		                MPI_Get_count(&status, MPI_DOUBLE, &recive_size);
		                
		                tmp = (double*)malloc(recive_size*sizeof(double));
		                
		                MPI_Recv(tmp, recive_size, MPI_DOUBLE, r, r, MPI_COMM_WORLD, NULL );
		                for(int i = 0; i < recive_size; i++)
					printf("Recieved %lf \n", tmp[i]);
		                memcpy(finaldata+currentlocation,tmp,recive_size*sizeof(double));
		                
		                currentlocation = currentlocation+recive_size;
		                                                        
		                free(tmp);
                        }
                }
        }
        
        if (rank==0){        
		double *finald = finaldata + padding;
		if(size != 1){
			
			write_output(output_name, finaldata, SIZE);
		} else {
			write_output(output_name, &local_arr[1], SIZE-1);
		}
		printf("%lf\n", end);
		free(input);
		free(finaldata);
        }

        free(local_arr);
                            

        MPI_Finalize();
}

double* Parallel_Qsort(MPI_Comm curr_group, int rank, int size, double* local_arr, int n){

        MPI_Comm_size( curr_group, &size);
        MPI_Comm_rank( curr_group, &rank );

        //Recursion statement 
        if(size < 2){ 
		for(int i = 0; i < n; i++)
			printf("Last Recursion %lf \n", local_arr[i]);
                return(local_arr);
        }
        
        double PivotPoint = pivot(STRAT, n, local_arr, rank, size, curr_group);
        
        //Find the size of arrays to share according to 
        double *top = (double*)malloc(n*sizeof(double));
        double *bot = (double*)malloc(n*sizeof(double));
        int len_top = 0;
        int len_bot = 0;
        for(int i = 0 ; i < n;i++){
                if(local_arr[i] >= PivotPoint){
                        top[len_top++] = local_arr[i];
                }else{
                        bot[len_bot++] = local_arr[i];
                }
        }
        
        free(local_arr);
        
        //Group up process ids into two groups, high and low, every process has a pair-process in other group
        
        int group_size = size /2;
        
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
                                
        double *recv_low_from_high;
        double *recv_high_from_low;
        MPI_Request request1;
        MPI_Request request2;

        MPI_Status status1;
        MPI_Status status2;
        
        if (rank < group_size) {
        //low process sends its top to pair in right group and receive pairs bottom
                MPI_Isend(top, len_top, MPI_DOUBLE, pair_high, tag+1, curr_group, &request1);
		for(int i = 0; i < len_top; i++)
			printf("Top sent %lf \n", top[i]);
                MPI_Probe(pair_high, tag, curr_group, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &recv_size_low);
                recv_low_from_high= (double *)malloc(recv_size_low*sizeof(double));
                MPI_Irecv(recv_low_from_high, recv_size_low, MPI_DOUBLE, pair_high, tag, curr_group, &request2);
		for(int i = 0; i < recv_size_low; i++)
			printf("Bot recv %lf \n", recv_low_from_high[i]);
                MPI_Wait(&request1, &status);
                MPI_Wait(&request2, &status);
        } else {
        //high process sends its bottom to pair in left group and receive pairs top
                MPI_Isend(bot, len_bot, MPI_DOUBLE, pair_low, tag, curr_group, &request1);
		for(int i = 0; i < len_bot; i++)
			printf("Bot sent %lf \n", bot[i]);
                MPI_Probe(pair_low, tag+1, curr_group, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &recv_size_high);
                recv_high_from_low= (double *)malloc(recv_size_high*sizeof(double));
                MPI_Irecv(recv_high_from_low, recv_size_high, MPI_DOUBLE, pair_low, tag+1, curr_group, &request2);
		for(int i = 0; i < recv_size_high; i++)
			printf("Top recv %lf \n", recv_high_from_low[i]);
                MPI_Wait(&request1, &status);
                MPI_Wait(&request2, &status);
        }
        MPI_Barrier(curr_group);
        
        double *merged_arr;
        int newlen;
        
        MPI_Comm nextgroup;
        int newsize,newrank;
        
        if (rank < group_size){
        	merged_arr = (double*) malloc((len_bot+recv_size_low)*sizeof(double));
                merge(bot,len_bot,recv_low_from_high,recv_size_low, &merged_arr);
		for(int i = 0; i < len_bot + recv_size_low; i++)
			printf("Merge 1 %lf \n", merged_arr[i]);
                newlen = len_bot + recv_size_low;
                currentsize = newlen;

                //qsort(merged_arr, newlen, sizeof(double), compare);
                
                free(bot);
                free(recv_low_from_high);
                

                
                	// Split group
                MPI_Comm_split( curr_group,0, 1, &nextgroup );
                MPI_Comm_size( nextgroup, &newsize );
                MPI_Comm_rank( nextgroup, &newrank );


                
                 //Recursive call

 	        return (Parallel_Qsort(nextgroup, newrank, newsize, merged_arr, newlen));
                                        
         }else{
         	merged_arr = (double*) malloc((len_top+recv_size_high)*sizeof(double));
        	merge(top,len_top,recv_high_from_low,recv_size_high, &merged_arr);
		for(int i = 0; i < len_top + recv_size_high; i++)
			printf("Merge 2 %lf \n", merged_arr[i]);
                                
         	newlen = len_top + recv_size_high;
         	currentsize = newlen;
                                
         	//qsort(merged_arr, newlen, sizeof(double), compare);

                free(top);
                free(recv_high_from_low);

         // printf("rank:%d \n")
                                
                // Split group
                MPI_Comm_split( curr_group, 1, 1, &nextgroup );
                MPI_Comm_size( nextgroup, &newsize );
                MPI_Comm_rank( nextgroup, &newrank );


                //Recursive call
                return (Parallel_Qsort(nextgroup, newrank, newsize, merged_arr, newlen));

                }

        return(local_arr);
}

double pivot(int pivot_strat,int len, double *local_arr, int rank, int size, MPI_Comm group){
        double l_pivot = 0;
        MPI_Comm_size( group, &size);
        MPI_Comm_rank( group, &rank );
        switch (pivot_strat){
        case 1 : /* 1. Select the median in one processor in each group of processors*/
        {
                /* ROOT find the Global PIVOT and broadcast it to all other PROCS */
                if (rank==0){
                                                l_pivot = find_pivot(local_arr, len);
                }

                MPI_Bcast ( &l_pivot, 1, MPI_DOUBLE, 0, group);
                
	return l_pivot;
        }
        break;
        case 2 : /* Select the median of all medians in each processor group */
        {
                double local_pivot = find_pivot(local_arr, len); /* Each PROC find median of local buffer */
                double all_median[size]; /* # of median values as our #PROCS */
                
                MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, group);
                l_pivot = find_pivot(all_median, size); /* median of all medians */

	return l_pivot;
        }
        break;
        case 3: /* Select the mean value of all medians in each processor group */
        {
                double local_pivot = find_pivot(local_arr, len); /* Each PROC find median of local buffer */
                double all_median[size]; /* # of median values as our #PROCS */
                MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, group);
                l_pivot = find_mean( all_median, size);
		return l_pivot;
        }
        break;
        }
        return 0;
}
//Compare double
int compare(const void * a, const void * b)
{
        if (*(double*)a > *(double*)b) return 1;
        else if (*(double*)a < *(double*)b) return -1;
        else return 0;
}

double find_mean(double *data, int len){
        double sum = 0;
        for (int i=0; i<len; ++i)
                sum += data[i];
        sum = sum/(double)len;
        return sum;
}

double find_pivot(double* data, int len){
        /* If the number of elements is zero the pivot is data[0] */
        if (len==0)
                return 0;
        double pivot = 0;
        /* If number of elements is odd, pivot is the middle element */
        if (len % 2 != 0){
                pivot = data[len/2];
        /* Else it is the average of the two middle elements */
        }else{
                pivot = ((double)data[(len / 2) -1] + (double)data[(len / 2)]) / 2.0;
	}
        return pivot;
}


double* merge(double *v1, int n1, double *v2, int n2, double** mergeTo){
        int i = 0;
        int j = 0;
        int k = 0;

        for(i = 0; i < n1; i++)
		printf("1: %lf\n", v1[i]);
        for(j = 0; j < n2; j++)
		printf("2: %lf\n", v2[j]);


	printf("\n\n");

         while(i<n1 && j<n2){
		if(v1[i]<v2[j])
		{
		        (*mergeTo)[k] = v1[i];
		       i++; k++;
		}
		else
		{
		        (*mergeTo)[k] = v2[j];
		        j++; k++;
		}
	}
	if(i==n1){
		while(j<n2)
		{
			(*mergeTo)[k] = v2[j];
			j++; k++;
		}
	}else{
		while(i<n1)
		{
			(*mergeTo)[k] = v1[i];
			i++; k++;
		}
	}
	return NULL;

/*
	for (int i = 0, j = 0, k = 0; k < n1+n2;) {
                if (v1[i] < v2[j])
                        result[k++] = v1[i++];
                else
                        result[k++] = v2[j++];
        }

        return result;*/
}
int read_input(const char *file_name,/* int rank, int size, int* num_value,*/ double **values) {
	
	/*MPI_File file;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, file_name,MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	
	if(rank==0){
		int temp;
		MPI_File_read_at(file,4, &temp, 1, MPI_INT, &status);
		printf("num vals = %i\n", temp);
	}*/
	/*MPI_Bcast(num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);
	printf("recived num vals = %i\n", num_values);
	int offset = sizeof(double)*rank*(*num_values)/size+sizeof(int);
	MPI_File_read_at(file, offset, values, (*num_values)/size, MPI_DOUBLE, &status);*/
	//MPI_File_close(&file);

	//for(int i = 0; i < (*num_values)/size; i++){
		//printf("Rank %i %lf\n", rank, values[i]);
	//}
	//return (*num_values);


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

	if (NULL == (*values = malloc((num_values) * sizeof(double)))) {
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
		if (0 > fprintf(file, "%.0lf ", output[i])) {
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

