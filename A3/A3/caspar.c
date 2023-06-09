#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int printlist(int *list, int n, int rank){
    printf("\nProcessor %d\n", rank);
    for(int i = 0; i<n; i++)
        printf("%d ", list[i]);
    printf("\n");
    return 0;
}

int printlisty(int *list, int n, int rank){
    printf("[");
    for(int i = 0; i<n; i++)
        printf("%d ", list[i]);
    printf("]");
    return 0;
}

int read_input(const char *file_name, int **values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "r"))) {
        perror("Couldn't open input file");
        return -1;
    }
    int n;
    if (EOF == fscanf(file, "%d", &n)) {
        perror("Couldn't read element count from input file");
        return -1;
    }
    if (NULL == (*values = malloc(n*sizeof(int)))) {
        perror("Couldn't allocate memory for input");
        return -1;
    }
    for (int i=0; i<n; i++) {
        if (EOF == fscanf(file, "%d", &((*values)[i]))) {
            perror("Couldn't read elements from input file");
            return -1;
        }
    }
    if (0 != fclose(file)){
        perror("Warning: couldn't close input file");
    }
    return n;
}


int write_output(char *file_name, const int *output, int num_values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "w"))) {
        perror("Couldn't open output file");
        return -1;
    }
    for (int i = 0; i < num_values; i++) {
        if (0 > fprintf(file, "%d ", output[i])) {
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

int ascend_comp(const void* a, const void* b) {
    // Compare two integers
    return *(int*)a - *(int*)b;
}

void mergeLists(int* result, int* a, int* b, int len_a, int len_b) {
    // Merge the sorted lists a and b into a single sorted list result
    int a_index = 0;
    int b_index = 0;
    int total_index = 0;
    while ((a_index < len_a) && (b_index < len_b)) {
        if (a[a_index] < b[b_index]) {
            result[total_index] = a[a_index];
            a_index += 1;
            total_index += 1;
        } else {
            result[total_index] = b[b_index];
            b_index += 1;
            total_index += 1;
        }
    }

    for (; a_index < len_a; a_index++) {
        result[total_index] = a[a_index];
        total_index += 1;
    }

    for (; b_index < len_b; b_index++) {
        result[total_index] = b[b_index];
        total_index += 1;
    }
    return;
}

int main(int argc, char **argv){
    int rank, p;
    double max_time;

	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);       

    int p_exp = log2(p);
    
    int n;
    int *input;
    int pivot_strategy;

    if(rank == 0) {
        char *input_name = argv[1];
        
        // Check if user input is of correct format.
        if (4 != argc) {
            printf("Usage: input_file, output_file, pivot strategy (1-3)\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        
        pivot_strategy = atoi(argv[3]);
        
        if(pivot_strategy < 1 || pivot_strategy > 3) {
            printf("Error: Only 3 available pivot stratergies. %d is invalid.\n", pivot_strategy);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        if (0 > (n = read_input(input_name, &input))) {
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

    }

    // Broadcast number of total elements
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast the pivot strategy
    MPI_Request pivot_strategy_request;
    MPI_Ibcast(&pivot_strategy, 1, MPI_INT, 0, MPI_COMM_WORLD, &pivot_strategy_request);

    int block_size = n/p;
    int excess_size = n % p;

    int* data = malloc((block_size+1)*sizeof(int));

    // Divide the excess elements between the first excess_size nodes.
    int send_lens[p];
    int send_displ[p];
    if (rank == 0) {
        for (int i = 0; i < p; i++) {
            if (i < excess_size) {
                send_lens[i] = block_size + 1;
            } else {
                send_lens[i] = block_size;
            }
        }
        send_displ[0] = 0;
        for (int i = 1; i < p; i++) {
            if (send_lens[i] > 0) {
                send_displ[i] = send_displ[i-1] + send_lens[i-1];
            } else {
                send_displ[i] = send_displ[i-1];
            }
        }
    }
    
    // Scatter the block sizes for each node
    MPI_Scatter(send_lens, 1, MPI_INT, &block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Scatter the data for each node
    MPI_Scatterv(input, send_lens, send_displ, MPI_INT, data, block_size, MPI_INT, 0, MPI_COMM_WORLD);
    


    //Sort locally:
    qsort(data, block_size, sizeof(int), ascend_comp);

    // Create all communicators before the actual computations start
    int dims[p_exp];
    int period[p_exp];
    int group_size[p_exp];

    for (int i = 0; i < p_exp; i++) {
        dims[i] = 2;
        period[i] = 1;
    }

    int reorder = 0;

    MPI_Comm comm;
    // Create hypercube cartesian communication 
    MPI_Cart_create(MPI_COMM_WORLD, p_exp, dims, period, reorder, &comm);

    MPI_Comm groups[p_exp];

    // The targets[i] is the rank of the node the current node will exchange elements with in iteration i.
    int targets[p_exp];

    // This node's coordinates in the hypercube.
    int my_coords[p_exp];
    MPI_Cart_coords(comm, rank, p_exp, my_coords);
    
    groups[0] = comm;
    group_size[0] = p;
    MPI_Comm last = comm;
    for (int dim = 0; dim < p_exp; dim++) {
        if (dim != p_exp-1) { // if not last iteration
            // Split group in half depending on which side of the hypercube they are on
            MPI_Comm_split(last, my_coords[dim], rank, &groups[dim + 1]);
            last = groups[dim + 1];
            group_size[dim+1] = group_size[dim]/2;
        }
        int rank_tmp = rank;
        // Calculate neighbour in dimension dim
        MPI_Cart_shift(comm, dim, 1, &rank_tmp, &targets[dim]);
    }

    MPI_Wait(&pivot_strategy_request, MPI_STATUS_IGNORE);

    // The buffer used to receive the elements from other nodes.
    int* recv_buffer = malloc(sizeof(int)*block_size);
    int recv_buffer_allocated = block_size;
    double start = MPI_Wtime();
    for (int i = 0; i < p_exp; i++) {
        int local_median = data[block_size/2];
        int group_pivot;
        
        if(pivot_strategy == 1) {
            // Use median of the middlemost node in group.
            int groups_middle = group_size[i]/2;
            group_pivot = local_median;
            MPI_Bcast(&group_pivot, 1, MPI_INT, groups_middle, groups[i]);
            
        }
        
        else if(pivot_strategy == 2) {
            // Median of medians
            int medians[group_size[i]];
            MPI_Allgather(&local_median, 1, MPI_INT, medians, 1, MPI_INT, groups[i]);
            // sort the medians
            qsort(medians, group_size[i], sizeof(int), ascend_comp);
            group_pivot = medians[group_size[i]/2];
        }
        
        else {
            // Mean of medians
            int medians[group_size[i]];
            MPI_Allgather(&local_median, 1, MPI_INT, medians, 1, MPI_INT, groups[i]);
            group_pivot = 0;
            for(int j = 0; j<group_size[i]; j++)
                group_pivot += medians[j];
            group_pivot = group_pivot/group_size[i];
        }
                
        int medIndex = 0; // index where `larger` array starts
        int completedLoop = 0;
        for (int j = 0; j < block_size; j++) {
            medIndex = j;
            if (data[medIndex] > group_pivot) {
                break;
            }
            if (j == block_size - 1) { // if we never break the loop
                completedLoop = 1;
            }
        }

        int lower_len = medIndex;
        if (completedLoop && (block_size > 0)) {
            lower_len = block_size;
        }
        int larger_len = block_size - lower_len;
        int* lower = data;
        int* larger = data + medIndex;

        int recv_len;
        int breakpoint;
        if (my_coords[i] == 0) { // I receive lower and send higher
            MPI_Sendrecv(&larger_len, 1, MPI_INT, targets[i], 222, &recv_len, 1, MPI_INT, targets[i], 222, comm, MPI_STATUS_IGNORE);
            int did_realloc = 0;
            if (lower_len + recv_len > recv_buffer_allocated) {
                // Realloc more than we need to not have to realloc if we only add a few more elements next time
                recv_buffer = realloc(recv_buffer, 2*(lower_len + recv_len)*sizeof(int));
                recv_buffer_allocated = 2*(lower_len + recv_len);
                did_realloc = 1;
            }
            block_size = lower_len + recv_len;
            // copy our kept data into recv_buffer
            memcpy(recv_buffer, lower, sizeof(int)*lower_len);
            MPI_Sendrecv(larger, larger_len, MPI_INT, targets[i], 333, recv_buffer + lower_len, recv_len, MPI_INT, targets[i], 333, comm, MPI_STATUS_IGNORE);
            if (did_realloc) {
                data = realloc(data, 2*(lower_len + recv_len)*sizeof(int));
            }
            breakpoint = lower_len;
        } else { // I receive higher and send lower
            MPI_Sendrecv(&lower_len, 1, MPI_INT, targets[i], 222, &recv_len, 1, MPI_INT, targets[i], 222, comm, MPI_STATUS_IGNORE);
            int did_realloc = 0;
            if (larger_len + recv_len > recv_buffer_allocated) {
                // Realloc more than we need to not have to realloc if we only add a few more elements next time
                recv_buffer = realloc(recv_buffer, 2*(larger_len + recv_len)*sizeof(int));
                recv_buffer_allocated = 2*(larger_len + recv_len);
                did_realloc = 1;
            }
            block_size = larger_len + recv_len;
            memcpy(recv_buffer, larger, sizeof(int)*larger_len);
            MPI_Sendrecv(lower, lower_len, MPI_INT, targets[i], 333, recv_buffer + larger_len, recv_len, MPI_INT, targets[i], 333, comm, MPI_STATUS_IGNORE);
            if (did_realloc) {
                data = realloc(data, 2*(larger_len + recv_len)*sizeof(int));
            }
            breakpoint = larger_len;
        }

        // Merge kept and received data
        mergeLists(data, recv_buffer, recv_buffer + breakpoint, breakpoint, recv_len);
    }
    double time = MPI_Wtime() - start;

    // Gather lengths of each node's list
    int lengths[p];
    int displs[p];
    MPI_Gather(&block_size, 1, MPI_INT, lengths, 1, MPI_INT, 0, comm);
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    displs[0] = 0;
    for (int i = 1; i < p; i++) {
        displs[i] = displs[i-1] + lengths[i-1];
    }

    // Gather the data froma all nodes
    MPI_Gatherv(data, block_size, MPI_INT, input, lengths, displs, MPI_INT, 0, MPI_COMM_WORLD);

    // Free memory
    for (int i = 0; i < p_exp; i++)
        MPI_Comm_free(&groups[i]);
    free(data);
    free(recv_buffer);

    // Write to output file
    if (rank == 0){
        printf("%f\n", max_time);
        char *output_name = argv[2];
        /*if (0 != write_output(output_name, input, n))
            return 2;*/
        free(input);
    }
    MPI_Finalize();

    return 0;
}

