#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG  0
int quicksort1(int *x, int left, int right);
int partition(int *x, int low, int high) ;
void swap(int *x,int a, int b);
int main(){
	clock_t t1, t2;				// variables for computing clocks 
	int *x, *y, s, p;
	double T1;
	int i, j, N=9;
	x = (int *) malloc( N * sizeof(int) );
	y = (int *) malloc( N * sizeof(int) );
	srand( time(NULL) );
	for(i=0;i<N;++i)
	{
		y[i] = x[i] = rand() % N;
	}
	for(i=0;i<N;++i)
	{
		printf("x[%d]=%d\n",i,x[i]);
	}
	printf("******partition********\n");
	printf("partiton= %d\n",partition(x, 0, N));
	for(i=0;i<N;++i)
	{
		printf("x[%d]=%d\n",i,x[i]);
	}
	return 0;
}
int quicksort1(int *x, int left, int right)
{
	int i, j, k, pivot, pivot_loc, N = right-left; 
	int *y;
	int middle=N/2;
	if(left < right-1)
	{
		y = (int *) malloc(N*sizeof(int));
	    /*
	  	pivot_loc = left+(rand() % N);
		pivot = x[pivot_loc];
		x[pivot_loc] = x[left];
		x[left] = pivot;
		*/
		pivot=x[left];
		i=0;j=N-1;
		for(k=1;k<N;++k)
		{
			if(x[left+k]<pivot)
			{
				y[i++]=x[left+k];
				//i=i+1;
			}
			else
			{
				y[j--]=x[left+k];
				//j=j-1;
			}
		}
		y[i]=pivot;
		#if DEBUG 
		printf("%d %d %d %d %d %d\n",left,right,i,j,pivot,N);
		for(k=0;k<N;++k)
		{
			printf("y[%d]=%d\n",k,y[k]);
		}
		#endif	
		for(k=0;k<N;++k)
		{
			x[left+k] = y[k];
		}
		free(y);
		
		
		quicksort1(x,left,left+i);
		
	
		quicksort1(x,left+i+1,right);
		
				
	}
	else 
	{
		return 1;
	}
	
}
int partition(int *x, int low, int high)
 {  
		
        int tmp = x[low];  
        int i = low, j = high;  
        while (i < j) {  
            while (i < j && x[j] >= tmp){  
                j--;  
            }  
            while (i < j && x[i] <= tmp){  
                i++;  
            }  
            swap(x, i, j);  
        }  
        x[low] = x[i];  
        x[i] = tmp; 
		 
        return i;  
} 
void swap(int *x, int i, int j) {  
        if(i == j){  
            return;  
        }  
        int tmp = x[i];  
        x[i] = x[j];  
        x[j] = tmp;  
    }

