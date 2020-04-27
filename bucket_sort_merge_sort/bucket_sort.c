
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"


// decides buckets based on number of processors
int bucket_search(double num, int p_num){
  int x;
  for(x=1; x < p_num+1; x++){
	double bucket_range =(double) x / (double)p_num;
	if(num <= bucket_range){
	  return x - 1; //return bucket number
	}
  }
return 0;
}
// Comparison function for sorting
int compare_sort(const void* arg1, const void* arg2){
 double a1 = *(double *) arg1;
 double a2 = *(double *) arg2;
 if (a1 < a2) return -1;
 else if (a1 == a2) return 0;
 else return 1;
}

// Sort array
void qsort_array(double *array, int array_len){
 qsort(array, (size_t)array_len, sizeof(double), compare_sort);
} 
int main(int argc, char *argv[]){

	if(argc != 2){
		printf("\n enter problem size \n");
		return 0;
	}
    
    int sub_count[1];
    
    	int p_rank, num_procs;

	//Allocate Arrays	
	int array_nums = strtol(argv[1], NULL, 10);

	//Init MPI, get process
	MPI_Init(&argc, &argv);
    //and ranks
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
	MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);	

	double t1 = MPI_Wtime();
	
	//elements on each processor
    
	int eachproc = array_nums/num_procs;
    
	double *list = (double*)malloc(eachproc*sizeof(double));
    
	int *count = (int*)calloc(num_procs, sizeof(int));
    
	int i, bucket;
	double r;
	for(i = 0; i < eachproc; i++){

		r = (double)rand() / (double) RAND_MAX;

		list[i] = r;
		bucket = bucket_search(r, num_procs);
		count[bucket]++;
	}

	int *bucket_count = (int*)malloc(num_procs*sizeof(int));
    
	MPI_Alltoall(count, 1, MPI_INT, bucket_count, 1, MPI_INT, MPI_COMM_WORLD);

	int loc_bcount = 0;
	//Add together counts
	for(i = 0; i < num_procs; i++){
		loc_bcount+= bucket_count[i]; 
	}

	//Allocate arrays based on counts
	double *bucket_list = (double*)malloc(loc_bcount*sizeof(double));

	//Distribute list to other processes
		int *displs = (int*)malloc(num_procs*sizeof(int));
	double *dist_list = (double*)malloc(eachproc*sizeof(double));
	int *index = (int*)calloc(num_procs,sizeof(int));

	displs[0] = 0;
	for(i = 1; i < num_procs; i++){
		displs[i] = count[i-1] + displs[i-1];
	}

	int *rdispls = (int*)malloc(num_procs*sizeof(int));
    
	rdispls[0] = 0;
	for(i = 1; i < num_procs; i++){
		rdispls[i] = bucket_count[i-1] + rdispls[i-1];
	}

	for(i = 0; i < eachproc; i++){
		//Find bucket for double
		bucket = bucket_search(list[i], num_procs);
		dist_list[displs[bucket] + index[bucket]] = list[i];
		index[bucket]++;
	}
	free(list);

	MPI_Alltoallv(dist_list, count, displs, MPI_DOUBLE, bucket_list, bucket_count, rdispls, MPI_DOUBLE, MPI_COMM_WORLD); 	

	//Quicksort on each list locally
	qsort_array(bucket_list, loc_bcount);

	//Gather counts of each bucket
	int gathercounts[1];
	gathercounts[0] = loc_bcount;
	MPI_Gather(gathercounts, 1, MPI_INT, count, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(p_rank==0){
	displs[0] = 0;
	for(i = 1; i < num_procs; i++){
		displs[i] = count[i-1] + displs[i-1];
	}
	}
	
	double* final_list = (double*)malloc(array_nums*sizeof(double));
	//Gather all lists at root 
	MPI_Gatherv(bucket_list,loc_bcount, MPI_DOUBLE, final_list, count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	double t2 = MPI_Wtime();

	//Check Result
	if(p_rank == 0){
        
		int sorted = 1;
		int k;
		for(k = 0; k < array_nums - 2; k++){
			if(final_list[k] > final_list[k+1]){
				sorted = 0;
			}
		}

		printf("\nV2- N: %d P: %d  Execution Time: %f\n",array_nums, num_procs, t2-t1);
	}

	//Free allocated Arrays
	free(index);
	free(displs);
	free(rdispls);
	free(count);
	free(bucket_count);
	free(bucket_list);
	free(final_list);
	MPI_Finalize();
  

}