#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG  0
int quicksort1(int *x, int left, int right);
int partition(int *x, int left, int right) ;
void swap(int *x,int a, int b);
int findK(int *x, int left, int right, int k, int prepart);
int main(){
	clock_t t1, t2;				// variables for computing clocks 
	int *x, *y,*z, s, p;
	double T1;
	int i, j, N;
	printf("please dicide the number of N:");
	scanf("%d", &N);
	printf("N= %d\n",N);
	x = (int *) malloc( N * sizeof(int) );
	y = (int *) malloc( N * sizeof(int) );
	z = (int *) malloc( N * sizeof(int) );
	srand( time(NULL) );
	for(i=0;i<N;++i)
	{
		z[i]=y[i]= x[i] = rand() % N;
	}
	for(i=0;i<N;++i)
	{
		printf("x[%d]=%d\n",i,x[i]);
	}
	int prepart=partition(y,0,N);
	for(i=0;i<N;++i) 
		{
			for(j=i+1;j<N;++j)
			{
				if(y[i]>y[j]) 
				{
					s = y[i];
					y[i] = y[j];
					y[j] = s;
				}
			}
		}
	printf("*********sorting********\n");
	for(i=0;i<N;++i)
	{
		printf("y[%d]=%d\n",i,y[i]);
	}
	int kv=N/2+1;
	double medium = 0;
	double medium1 = 0;
	double medium2 = 0;
	if(N%2==1)
	{
		medium=findK(x, 0,N,kv,prepart);
		printf("medium= %f\n", medium);
	}
	else if(N%2==0)
	{
		medium1=findK(x, 0,N,kv,prepart);
		printf("medium1= %f\n", medium1);
		medium2=findK(z, 0,N,kv-1,prepart);
		printf("medium2= %f\n", medium2);
		medium=(medium1+medium2)/2;
		printf("medium= %f\n", medium);
		
		
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
int partition(int *x, int left, int right)
 {  
    int i, j, k;
	int pivot, t;
	
	//  
	if(left < right-1)
	{
		pivot = x[left];
    	i = left+1;
    	j = right-1;
    	// 
    	while(1)
		{
			while(i < right && pivot >= x[i]) i++; // 往右邊找到第一個  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // 往左邊找到第一個  pivot >= x[j] 
      		#if DEBUG
			printf("%d %d %d\n", i,j,pivot);
			#endif
      		if(i>=j) break;
      		swap(x,i,j);
      		#if DEBUG
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			#endif
		} 
        x[left] = x[j];
        x[j] = pivot;
		 
        return j;  
	} 
}
int findK(int *x, int left, int right, int k, int prepart)
{  
		int i;
        if(left > right-1)
		{  
            return x[prepart];  
        }  
        int pos = partition(x, left, right);  
        int leftnumber = pos -left  + 1;//左??? 
		//printf("pos=%d,prepart=%d,leftnumber=%d,k=%d\n",pos,prepart,leftnumber,k);
		//system("pause"); 
        if(k > leftnumber)
		{  
		   //中位數在pivot右邊 找 pivot 右邊 
            findK(x, pos +1, right, k - leftnumber, pos);  
        }  
        else if(k < leftnumber)
		{   
			//中位數在 pivot左邊 找pivot 左邊    前一步再做一次  
            findK(x, left, pos -1, k, prepart);  
        } 
        else
		{
			return x[pos];
		}
}
void swap(int *x, int i, int j) 
{  
        if(i == j){  	
            return;  
        }  
        int tmp = x[i];  
        x[i] = x[j];  
        x[j] = tmp;  
}

