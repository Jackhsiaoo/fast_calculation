#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

int main()
{
	int i, j,k, N=100000000;
	double a = 1.234, b = 2.345;
	double ot1, ot2;
	clock_t t1, t2;
	ot1 = omp_get_wtime();
	#pragma omp parallel num_threads(2)  //parallel calculation for 4 cpu
	{
		printf("Hello world (%d,%d,%f) \n",omp_get_thread_num(),omp_get_num_threads(),omp_get_wtime()-ot1);
		printf("Hello Program\n");
	}
	#pragma omp parallel for  //add for means give ten circulation to all parallel
	                           //no for are each parallel has ten circulation
	for(i=0;i<10;++i)
	{   
		
		j = i;                 //j=0,j=1....j=9
		printf("%d %d\n",i,omp_get_thread_num());
	}
	printf("j = %d\n",j);
	
	
	j = 0;
	#pragma omp parallel for reduction(+:j)  //reduction for j is to make each cpu not get other's memory
	for(i=0;i<=10;++i)
	{
		j += i;
		//printf("i=%d,j=%d,threat()=%d\n",i,j,omp_get_thread_num());
	}
	printf("sum(1..10) = %d\n",j);

	j = 0;
	#pragma omp parallel for private(i) reduction(+:j)
	for(i=0;i<=100;++i)
	{
		j += i;
		
	}
	printf("sum(1..100) = %d\n",j);

	
	return 0;
}
