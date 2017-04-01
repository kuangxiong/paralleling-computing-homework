/*************************************************************************
	> File Name: parallelCG.c
	> Author:kuangxiong 
	> Mail:kuangxiong@lsec.cc.ac.cn 
	> Created Time: 2017年03月30日 星期四 15时12分34秒
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"mpi.h"
#define n 4

int 
main(int argc, char* argv[])
{
   int i, j, k, rank, size, blocksize;
   double *A, *b, *x0, *x, *r, alpha, beta, s, *w, *p, rho_0, rho_1, rho, tols;

   MPI_Init(NULL, NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if(n % size ==0)
	   blocksize = n/size;
   else if(rank < size-1)
	   blocksize = n/(size-1);
   else
	   blocksize = n%(size -1);

   A = malloc(blocksize * n * sizeof(double));
   b = malloc(n * sizeof(double));
   x0 = malloc(n * sizeof(double));
   p = malloc(blocksize * sizeof(double));
   r = malloc(blocksize * sizeof(double));
   w = malloc(blocksize * sizeof(double));
   x = malloc(blocksize * sizeof(double));
//分发数据
   for(i = 0; i< blocksize; i++)
	   for(k = 0; k< n; k++)
	   {   if(k == rank * blocksize + i)
				A[i* n + k] = i + k + n + 3;
           else 
		        A[i * n + k] = 1.0; 
		}
   for(i = 0; i< blocksize; i++)
	{	
		b[i] = rank *blocksize + i;
        x[i] = 1.0;
	}
   for(i = 0; i< n; i++)
    {  
        x0[i] = 1.0;
	}
   for( i =0 ; i< blocksize; i++)
   {	for(j=0 ;j< n; j++)
		   printf("%f\t", A[i* n +j]);
	printf("\n");
   }
   for(i = 0; i< blocksize; i++)
   {   
	   r[i] = 0;
	   for(j = 0; j< n; j++)
	        r[i] += A[i * n +j] *x0[j];
    }
	for(k = 0; k< blocksize; k++)
	    	r[k] = b[k] - r[k];

    rho = 0;
    for( i =0; i< blocksize; i++)
		rho += r[i] * r[i];

	MPI_Allreduce(&rho, &rho_0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(i = 0; i< blocksize; i++)
		p[i] = r[i];

	for(i =0 ; i< blocksize; i++)
	{  w[i] = 0;
	   for(j = 0 ; j< n; j++)
		   w[i] += A[i* n + j]*p[j];
	}
    s = 0;
	for(i=0 ; i< blocksize; i++)
		s += p[i]*w[i];
	MPI_Allreduce(&s, &tols, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    alpha = rho_0/tols;
	for(i=0 ; i< blocksize; i++)
	{
		x[i] = x[i] + alpha * p[i];
	    r[i] = r[i] - alpha * w[i];
	}
	s = 0;
    for(i=0 ; i< blocksize; i++)
       s += r[i]*r[i];
	MPI_Allreduce(&s, &rho_1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	printf("rho_1: %f\n", rho_1);
    while(rho_1 > 0.0000001)
	{
		beta = rho_1/rho_0;
		rho_0 = rho_1;
		for(k = 0; k< n; k++)
			p[k] = r[k] + beta * p[k];
		for(k = 0; k< blocksize; k++)
		{
		   w[k] =0; 
		   for(j=0 ; j< n; j++)
			   w[k] += A[k* n +j] *p[j];
		}
		s = 0;
        for(i=0 ; i< blocksize; i++)
			s += p[i] *w[i];
		MPI_Allreduce(&s, &tols, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		alpha = rho_0/tols;
		rho = 0;
		for(i=0 ; i< blocksize; i++)
		{
			x[i] = x[i] + alpha * p[i];
			r[i] = r[i] - alpha * w[i];
		}
		s = 0;
		for(i=0 ; i< blocksize; i++)
			s += r[i]*r[i];
		MPI_Allreduce(&s, &rho_1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		printf("rho_1:%f\n", rho_1);
	}
	for(i=0; i< blocksize; i++)
		printf("%f\n", x[i]);
   printf("size: %d\t, rank:%d\t, blocksize:%d\n", size, rank, blocksize);
   free(r);
   free(A);
   free(p);
   free(w);
   free(x);
//   free(x0);
   MPI_Finalize();
   return 0;
}
