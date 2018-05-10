#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#define m 24 
int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N);

int main()
{
	int i;
	double y_re[m], y_im[m], x_re[m], x_im[m];
	for(i=0;i<m;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	Fast_Fourier_Transform(y_re, y_im, x_re, x_im, m);
	for(i=0;i<m;++i)
	{
		printf("%f + %f i\n", y_re[i], y_im[i]);
	}
	
	 
}
int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{
	if(N==2) 
	{
		// y, y[0] = x[0]+x[1], y[1] = x[0] - x[1]
		y_re[0] = x_re[0] + x_re[1];
		y_im[0] = x_im[0] + x_im[1];
		y_re[1] = x_re[0] - x_re[1]; 
		y_im[1] = x_im[0] - x_im[1];
	} 
	else if(N==3)
	{
		y_re[0]= x_re[0] + x_re[1] + x_re[2];
		y_im[0]= x_im[0] + x_im[1] + x_im[2];
		y_re[1]= x_re[0] - 0.5*x_re[1] + 0.5*sqrt(3)*x_im[1] - 0.5*x_re[2] - 0.5*sqrt(3)*x_im[2];
		y_im[1]= x_im[0] - 0.5*sqrt(3)*x_re[1] - 0.5*x_im[1] + 0.5*sqrt(3)*x_re[2] - 0.5*x_im[2];
		y_re[2]= x_re[0] - 0.5*x_re[1] - 0.5*sqrt(3)*x_im[1] - 0.5*x_re[2] + 0.5*sqrt(3)*x_im[2];
		y_im[2]= x_im[0] + 0.5*sqrt(3)*x_re[1] - 0.5*x_im[1] - 0.5*sqrt(3)*x_re[2] - 0.5*x_im[2];
	}
	else
	{ 
	  if(N%2==0)
	  {
	  	
			int k;
			double *y_even_re, *y_even_im, *y_odd_re, *y_odd_im;
			double *x_even_re, *x_even_im, *x_odd_re, *x_odd_im;
			double w_re, w_im, w_N_re, w_N_im, a, b, temp;
			y_even_re = (double *) malloc( N/2 * sizeof(double));
			y_even_im = (double *) malloc( N/2 * sizeof(double));
			x_even_re = (double *) malloc( N/2 * sizeof(double));
			x_even_im = (double *) malloc( N/2 * sizeof(double));
			y_odd_re = (double *) malloc( N/2 * sizeof(double));
			y_odd_im = (double *) malloc( N/2 * sizeof(double));
			x_odd_re = (double *) malloc( N/2 * sizeof(double));
			x_odd_im = (double *) malloc( N/2 * sizeof(double));
			for(k=0;k<N/2;++k)
			{
				x_even_re[k] = x_re[2*k];
				x_even_im[k] = x_im[2*k];
				x_odd_re[k]  = x_re[2*k+1];
				x_odd_im[k]  = x_im[2*k+1];
			}
			Fast_Fourier_Transform(y_even_re, y_even_im, x_even_re, x_even_im, N/2);
			Fast_Fourier_Transform(y_odd_re, y_odd_im, x_odd_re, x_odd_im, N/2);
			// y_k = even_k + w_N^k odd_k = even_k + (a + bi)
			w_N_re =  cos(2.0*M_PI/N);
			w_N_im = -sin(2.0*M_PI/N);
			w_re   = 1.0;
			w_im   = 0.0; 
			for(k=0;k<N/2;++k)
			{
				a = w_re*y_odd_re[k] - w_im*y_odd_im[k];
				b = w_re*y_odd_im[k] + w_im*y_odd_re[k];
				y_re[k]     = y_even_re[k] + a;
				y_im[k]     = y_even_im[k] + b;
				y_re[N/2+k] = y_even_re[k] - a;
				y_im[N/2+k] = y_even_im[k] - b;
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
			free(y_even_re);
			free(x_even_re);
			free(y_even_im);
			free(x_even_im);
			free(y_odd_re);
			free(y_odd_im);
			free(x_odd_re);
			free(x_odd_im);
		}
		else if(N%3==0)
		{
			int k;
			double *y_zero_re, *y_zero_im, *y_one_re, *y_one_im, *y_two_re, *y_two_im;
			double *x_zero_re, *x_zero_im, *x_one_re, *x_one_im, *x_two_re, *x_two_im;
			double w_re, w_im, w_N_re, w_N_im, a1, b1, a2, b2, temp;
			double p,q,r;
			y_zero_re = (double *) malloc( N/3 * sizeof(double));
			y_zero_im = (double *) malloc( N/3 * sizeof(double));
			x_zero_re = (double *) malloc( N/3 * sizeof(double));
			x_zero_im = (double *) malloc( N/3 * sizeof(double));
			y_one_re = (double *) malloc( N/3 * sizeof(double));
			y_one_im = (double *) malloc( N/3 * sizeof(double));
			x_one_re = (double *) malloc( N/3 * sizeof(double));
			x_one_im = (double *) malloc( N/3 * sizeof(double));
			y_two_re = (double *) malloc( N/3 * sizeof(double));
			y_two_im = (double *) malloc( N/3 * sizeof(double));
			x_two_re = (double *) malloc( N/3 * sizeof(double));
			x_two_im = (double *) malloc( N/3 * sizeof(double));
			for(k=0;k<N/3;++k)
			{
				x_zero_re[k] = x_re[3*k];
				x_zero_im[k] = x_im[3*k];
				x_one_re[k]  = x_re[3*k+1];
				x_one_im[k]  = x_im[3*k+1];
				x_two_re[k]  = x_re[3*k+2];
				x_two_im[k]  = x_im[3*k+2];
			}
			Fast_Fourier_Transform(y_zero_re, y_zero_im, x_zero_re, x_zero_im, N/3);
			Fast_Fourier_Transform(y_one_re, y_one_im, x_one_re, x_one_im, N/3);
			Fast_Fourier_Transform(y_two_re, y_two_im, x_two_re, x_two_im, N/3);
			// y_k = zero_k + w_N^k one_k + w_N^2k two_k = zero_k + (a1 + b1i)+ (a2+b2i)
			w_N_re =  cos(2.0*M_PI/N);
			w_N_im = -sin(2.0*M_PI/N);
			w_re   = 1.0;
			w_im   = 0.0; 
			for(k=0;k<N/3;++k)
			{   
				
				// y_k = zero_k + w_N^k one_k + w_N^2k two_k = zero_k + (a1 + b1i)+ (a2+b2i)
				// w_N^2k= w_re*w_re - w_im*w_im + 2*w_re*w_im i
				a1 = w_re*y_one_re[k] - w_im*y_one_im[k];
				b1 = w_re*y_one_im[k] + w_im*y_one_re[k];
				a2 = ( w_re*w_re - w_im*w_im )*y_two_re[k] - (2*w_re*w_im )*y_two_im[k];
				b2 = ( w_re*w_re - w_im*w_im )*y_two_im[k] + (2*w_re*w_im )*y_two_re[k];
				//p%3==0
				y_re[k]     = y_zero_re[k] + a1 + a2;   
				y_im[k]     = y_zero_im[k] + b1 + b2;
				//q%3==1
				y_re[k+N/3] = y_zero_re[k] - 0.5*a1 + 0.5*sqrt(3)*b1 - 0.5*a2 - 0.5*sqrt(3)*b2; 
				y_im[k+N/3] = y_zero_im[k] - 0.5*sqrt(3)*a1 - 0.5*b1 + 0.5*sqrt(3)*a2 - 0.5*b2;
				//r%3==2
				y_re[k+2*N/3]= y_zero_re[k] - 0.5*a1 - 0.5*sqrt(3)*b1 - 0.5*a2 + 0.5*sqrt(3)*b2;
				y_im[k+2*N/3]= y_zero_im[k] + 0.5*sqrt(3)*a1 - 0.5*b1 - 0.5*sqrt(3)*a2 - 0.5*b2;
				
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
			free(y_zero_re);
			free(x_zero_re);
			free(y_zero_im);
			free(x_zero_im);
			free(y_one_re);
			free(y_one_im);
			free(x_one_re);
			free(x_one_im);
			free(y_two_re);
			free(y_two_im);
			free(x_two_re);
			free(x_two_im);
		}
	}
}
