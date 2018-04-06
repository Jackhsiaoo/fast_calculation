#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
int function(int n)
{
	if(n==1)
	{
		return 1;  // it must have be written 
	} else
	{
		return n*(function(n-1));
	}
}

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	double **A, *x, *b, *c, T1;		// x = 1234, a location in memory, x[0] the value at memory 1234, 1235, 1236, ... 1041
	int **C;
	int i, j, k, L, M, N=4; //
	j=0;
	

	////******************///////////////
	A = (double **) malloc( N * sizeof(double*) ); // 在記憶體中拿N個double*記憶體, 
													   // 並且把一開始指標放在A, A[0], A[1], A[2], A[N-1]
													   // 但是  A[0], A[1], A[2], A[N-1] 的值是未給定的 
	A[0] = (double *) malloc( N*N*sizeof(double)); // 在記憶體中拿N*N個double記憶體, 
													   // 並且把一開始指標放在A[0] 
													   // A[0][0], A[0][1], .... A[0][N*N-1] 尚未給定 
	for(i=1;i<N;++i) A[i] = A[i-1] + N; // A[1] = A[0]+N ==> A[1][0] = A[0][N]
	x = (double *) malloc( N * sizeof(double) );
	b = (double *) malloc( N * sizeof(double) );
	C = (int **) malloc(N*sizeof(int*));
	C[0] = (int *) malloc(N*N*sizeof(int));
	for(i=1;i<N;++i) C[i] = C[0]+i*N; 	
	
	M = N/4;
		#pragma omp parallel num_threads(4) private(k)
		{
			k = omp_get_thread_num();
			printf("the seed at thread %d is : %d\n",k,time(NULL)>>k);
			srand(time(NULL)>>k);			// 在每一個 thread 中設定起始值 

			#pragma omp parallel for // 再往下做取亂數 
			for(i=k*M;i<(k+1)*M;++i)
			{
				//L = omp_get_thread_num();
				//printf("thread %d, %d\n",k,L);
				for(j=0;j<N;++j)
				{
					A[i][j] = rand() % N*N;
				}
				x[i] = rand() % N*N;
			}
		}
		double t;
		for(i=0;i<N;++i) 
		{
			t = 0.0;
			for(j=0;j<N;++j)
			{
				t += A[i][j]*x[j];
			}
			b[i] = t;
		}

		//output
		printf("A Matrix:\n");
		for(i=0;i<N;++i) 
		{ 
			for(j=0;j<N;++j)
			{
				printf("%f ",A[i][j]);
			}
			printf("\n");
		}
		printf("X Matrix:\n");
		for(j=0;j<N;++j) 
		{ 
		
				printf("%f\n",x[j]);
		}
		M=function(N);
		printf("HELLO %d",function(N));
		

	return 0;
} 
