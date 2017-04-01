#include<stdlib.h>
void CGmethod(double *A, int N,double x[],double b[])
{ /***** 注意 x[]表示初始的值  
  **** n 表示矩阵的维数。
  ****b[]表示右端项  
  *****p[]为任意的矩阵（函数中利用 P[]产生方向向量）。
  *****r[]表示为残差向量 */ 
	int i,j,k;   
	double alpha, beta, s1, s2, *r, *p;  
    r = malloc(N * sizeof(double));
	p = malloc(N * sizeof(double));

	for(i=0; i<N; i++)
		r[i] = 0;
	for(k = 0; k < N; k++)  
	{  for(i = 0; i < N; i++)      
			r[k] = r[k] + A[k*N+i] * x[i];       
			r[k] = b[k] - r[k];  
	}  
	for(k = 0; k < N; k++)  
	 p[k] = r[k]; 
	
   /*迭代过程*/ 
	for(k=0; k<N; k++)  
	{	 s1 = s2 = 0;  
		for(i = 0; i < N; i++) 
			s1=s1 + r[i] * p[i];    
		for(i=0; i<N; i++)    
			for(j=0; j<N; j++)   
				s2 = s2 + p[i] * A[i*N+j] * p[j];  
		alpha = s1/s2;     
		s1 = s2 = 0;
	/*迭代点 x[N]的更新*/    
		for(i=0; i<N; i++)   
			x[i] = x[i] + alpha*p[i]; /*残量 r[N]的更新 */   
		for(i=0; i<N; i++)       
			for(j=0; j<N; j++)      
				r[i] = r[i]-alpha*A[i*N+j]*p[j];
	 /*若残差=0，说明已找到精确解，退出循环*/ 
	for(i=0;i<N;i++)     
		 for(j=0;j<N;j++)      
			   s2 = s2 + p[i]*A[i*N+j]*p[j];    
	for(i=0; i<N; i++)    
		  for(j=0; j<N; j++)       
			   s1 = s1-p[i]*A[i*N+j]*r[j];    
		 beta = s1/s2;
	 //方向的更新
	 for(i=0; i<N; i++)
		 p[i] = r[i] + beta*p[i];     
	}
	free(r);
	free(p);
}
