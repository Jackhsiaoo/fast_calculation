#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	double **A, *x, *b, T1;
	int i, j, N=10;
	char *c;
	int n;
	n = 0x04010302;							// c[0]=(02) c[1]=(03) c[2]=(01) c[3]=(04)
	c = &n;

	//srand(time(NULL));    //produce randon different variables
	srand(0);
	
	
	printf("-----------------");

	
		A = (double **) malloc( N * sizeof(double*) );
		A[0] = (double *) malloc( N*N*sizeof(double));
		for(i=1;i<N;++i) A[i] = A[i-1] + N;
		x = (double *) malloc( N * sizeof(double) );
		b = (double *) malloc( N * sizeof(double) );
		
		for(i=0;i<N;++i)
		{
			for(j=0;j<N;++j)
			{
				A[i][j] = rand();
			}
			x[i] = rand();
		}
		t1 = clock();
		for(i=0;i<N;++i) 
		{
			b[i] = 0.0;
			for(j=0;j<N;++j)
			{
				b[i] += A[i][j]*x[j];
			}
		}
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Matrix time vector :%f\n",T1);
		printf("each location of n: %d %d %d %d\n", c[0], c[1], c[2], c[3]); // c[0]: c這個記憶體位置的值
		printf("*c %d\n", *c);	
		float *f = &n;								// sizeof(float) = 4,
	
		// (0111 1111 1000 0000 0000 0000 ...0)
		//     7 F    8    0 
	
		n = 0x7F800000;   // n ( 00 00 80 7F )
		printf("n = %d, f = %e\n", n, *f);

	return 0;
} 
