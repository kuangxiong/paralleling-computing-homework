/*************************************************************************
	> File Name: main.c
	> Author:kuangxiong 
	> Mail:kuangxiong@lsec.cc.ac.cn 
	> Created Time: 2017年03月29日 星期三 20时34分53秒
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<sys/time.h>
//#include"CG.c"
#include"newCG.c"
#define n 4

int main(int argc, char* argv[])
{
	int i, j;
	double *A, *b, *x, tol = 0.0000001;
    struct timeval start,end;  

	x = malloc(n*sizeof(double));
    A = malloc(n*n*sizeof(double));
	b = malloc(n*sizeof(double));
	for(i=0 ;i< n; i++)
	{	for(j=0 ;j< n; j++)
		{	
			if(i==j)
				A[i*n+j] = i+j+n+3;
			else 
				A[i*n+j] = 1.0;
//			printf("%f\t", A[i*n+j]);
		}
//		printf("\n");
		b[i] = i;
	}
	gettimeofday(&start, NULL );  
	CGmethod(A, n, x, b, tol);
    for(i=0; i< n; i++)
		printf("x:%f\n", x[i]);
    gettimeofday(&end, NULL );  
    long timeuse =1000000 * ( end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;  
	printf("time=%f(s)\n",timeuse /1000000.0); 
    
	free(A);
	free(b);
	free(x);
}
