compile: bucket_sort.o merge_sort.o

bucket_sort.o: bucket_sort.c
	mpicc -o bucketsort bucket_sort.c
merge_sort.o: merge_sort.c
	mpicc -o mergesort merge_sort.c
run_bucket_sort: bucket_sort.o
	mpiexec -np 4 --oversubscribe ./bucketsort 10000
run_merge_sort: merge_sort.o
	mpiexec -np 4 --oversubscribe ./mergesort 10000
clean:
	rm -rf bucketsort mergesort
