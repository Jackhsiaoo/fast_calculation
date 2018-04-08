#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
int quicksort(int *x, int left, int right);
void swap(int *a, int *b)
{
int tmp = *a;
*a = *b;
*b = tmp;
}

int main()
{
	int x[]={5,6,9,7,0,3,4,1};
	int i, j, k;
	int pivot, t;
	int left=0;
	int right=8;
	
	if(left < right-1)
	{
		pivot = x[left];
    	i = left;
    	j = right;
    	// 56970341
    	while(i<j)
		{
      		while(pivot >= x[i]) i++; 
      		while(pivot <= x[j]) j--; 
      		printf("%d %d %d\n", i,j,pivot);
      		if(i>=j) break;
      		swap(&x[i],&x[j]);
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			
        }
        swap(&x[left],&x[j]);
        printf("%d\n",j);
		for(k=left;k<right;++k)
		{
			printf("x[%d]=%d\n",k,x[k]);
		}
		
	}
	system("pause");
	return 0;
}


