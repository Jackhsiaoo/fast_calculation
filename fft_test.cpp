#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define DEBUG 0

int main()
{
	int    k, n, N = 1<<15;  // (00000001)=1  (00001000)=8
	double ca, sa, a, theta, c, s, t, *X_re, *X_im, *x_re, *x_im;
	double *ck1,*sk1,*ck,*sk,*temp;
	clock_t t1, t2;
	
	x_re = (double *) malloc(N*sizeof(double));
	x_im = (double *) malloc(N*sizeof(double));
	X_re = (double *) malloc(N*sizeof(double));
	X_im = (double *) malloc(N*sizeof(double));
	ck1  =(double *) malloc(N*sizeof(double));
	sk1  =(double *) malloc(N*sizeof(double));
	ck =(double *) malloc(N*sizeof(double));
	sk =(double *) malloc(N*sizeof(double));
	temp=(double *) malloc(N*sizeof(double));
	
	for(n=0;n<N;++n)
	{
		x_re[n] = n+1;
		x_im[n] = 0.0;
	}
	
	t1 = clock();
	for(k=0;k<2;++k)
	{
		X_re[k] = 0.0;
		X_im[k] = 0.0;
		a  = 2*M_PI*k/N;
		ca = cos(a);
		sa = sin(a);
		c=1.0;
		s=0.0;
		for(n=0;n<N;++n)
		{
			X_re[k] += x_re[n]*c + x_im[n]*s;
			X_im[k] += x_im[n]*c - x_re[n]*s;
			if(k==1)
			{
			ck[n]=ck1[n]=c;
			sk[n]=sk1[n]=s;
			//ck[n]  = ck1[n];   //先丟給k=2要跑得k=1 
			//sk[n] = sk1[n];
			//printf("n[%d],ck1=%fsk1=%f\n",n,ck1[n],sk1[n]);
			}
			
			t =  c;
			c =  c*ca - s*sa;
			s =  s*ca + t*sa; 
			
		}
	}
	for(k=2;k<N;++k)
	{
		X_re[k] = 0.0;
		X_im[k] = 0.0;
		//   cos(2*PI*k*n/N) n = 0, 1, 2, 3...
		//   theta = n * (2*PI*k/N) = n * a
		//   cos(n a) -> cos((n+1)a) = cos(na) cos(a) - sin(na) sin(a)
		//   sin(n a) -> sin((n+1)a) = sin(na) cos(b) + cos(na) sin(b)
		//a  = 2*M_PI*k/N;
		
		/*if(k==2)
		{
			for(n=0;n<N;++n)
			{
			ck[n]  = ck1[n];
			sk[n] = sk1[n];
			}
		}*/
		
		for(n=0;n<N;++n)
		{
			temp[n] =  ck[n];
			ck[n] =  ck[n]*ck1[n] - sk[n]*sk1[n];
			sk[n] =  sk[n]*ck1[n] + temp[n]*sk1[n]; 
			// x_n = x_re + i * x_im
			// x_n (cos(...)-i sin(...)
			// real: x_re * cos(...) + x_im * sin(...)
			// imag: x_im * cos(...) - x_re * sin(...)
			// theta = 2*M_PI*k*n/N;
			// c = cos(theta) = cos(na);
			// s = sin(theta) = sin(na);
			X_re[k] += x_re[n]*ck[n] + x_im[n]*sk[n];
			X_im[k] += x_im[n]*ck[n] - x_re[n]*sk[n];
			//   cos((n+1)a) = cos(na) cos(a) - sin(na) sin(a)
			//   sin((n+1)a) = sin(na) cos(a) + cos(na) sin(a)
		}
		
	}
	t2 = clock();
	
	printf("N = %d, Time: %f\n", N, 1.0*(t2-t1)/CLOCKS_PER_SEC);
	
	#if DEBUG
	for(k=0;k<N;++k)
	{
		if(X_im[k] >=0 )
			printf("%f + %f i\n", X_re[k], X_im[k]);
		else
			printf("%f - %f i\n", X_re[k], -X_im[k]);
			
	}
	#endif
	return 0;
} 
