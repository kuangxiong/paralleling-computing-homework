#include<stdlib.h>
void CGmethod(double *A, int N, double x[], double b[], double tol)
{ /***** 注意 x[]表示初始的值  
  **** n 表示矩阵的维数。
  ****b[]表示右端项  
  *****p[]为任意的矩阵（函数中利用 P[]产生方向向量）。
  *****r[]表示为残差向量 */ 
	int i, j, k;   
	double alpha, beta, s1, *r, *p, rho_0, rho_1, *w;  
    
	rho_0 = 0;
	rho_1 = 0;
	w = malloc(N * sizeof(double));
	r = malloc(N * sizeof(double));
	p = malloc(N * sizeof(double));

	for(i=0; i<N; i++)
	{    x[i] = 1.0;
	     r[i] = 0;
		 w[i] = 0;
	}
	for(k = 0; k < N; k++)  
	  for(i = 0; i < N; i++)      
	  	r[k] += A[k* N +i] * x[i];
	  
	for(k = 0; k < N; k++)  
	{   
		r[k] = b[k] - r[k];
		rho_0  = rho_0 + r[k]*r[k]; 
	    p[k] = r[k];
	}
	for(k = 0; k < N; k++)  
	    for(j=0; j< N; j++)
			w[k] += A[k* N + j] * p[j];
	s1 = 0;
	for(i = 0; i < N; i++) 
			s1 = s1 + p[i] * w[i];  
	printf("s1: %f\n", s1);
	alpha = rho_0/s1;     
	for(i=0; i < N; i++)   
	{	x[i] = x[i] + alpha*p[i]; /*残量 r[N]的更新 */    
		r[i] = r[i]-alpha*w[i];
	    rho_1 = rho_1 + r[i]*r[i];
	} 
	printf("rho_1:%f\n", rho_1);
	//迭代的过程
	while(rho_1 > tol)
	{
		beta = rho_1/rho_0;
		rho_0 = rho_1;
		for(k = 0; k < N; k++)  
		    p[k] = r[k] + beta * p[k];
		for(k = 0; k < N; k++) 
		{	w[k] = 0;
			for(j=0; j< N; j++)
				w[k] += A[k* N + j] * p[j];
		}
		s1 = 0;
		for(i = 0; i < N; i++) 
			s1 = s1 + p[i] * w[i];  

		alpha = rho_0/s1;     
		/*迭代点 x[N]的更新*/    
		rho_1 = 0;
		for(i=0; i < N; i++)   
		{	x[i] = x[i] + alpha * p[i]; /*残量 r[N]的更新 */    
			r[i] = r[i] - alpha * w[i];
		    rho_1 = rho_1 + r[i]*r[i];
		} 
		printf("rho_1:%f\n", rho_1);
	}
	free(r);
	free(p);
	free(w);
}
