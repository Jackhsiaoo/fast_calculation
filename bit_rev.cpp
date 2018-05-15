#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define X 6 
int main()
{
	int i,m,p,q,k,N;
	double x[3], c[3];
	double a[3]={15,5,1};
	double t;
	for(i=0;i<3;++i)
	{
		x[i] = i;
		c[i] =0.0;
	}
	N=30;
    p=29;
    while(p<N & p>5)
    {
    	for(i=0;;i++)
    	{
    		c[i]=p/a[i];
		}
		p=p/a[i];
	}
	for(i=0;i<3;++i)
	{
		printf("c[i]=%f\n",c[i]);
	}
	
}
