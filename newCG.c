#include<stdlib.h>
void CGmethod(double *A, int N, double x[], double b[], double tol)
{ /***** ע�� x[]��ʾ��ʼ��ֵ  
  **** n ��ʾ�����ά����
  ****b[]��ʾ�Ҷ���  
  *****p[]Ϊ����ľ��󣨺��������� P[]����������������
  *****r[]��ʾΪ�в����� */ 
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
	{	x[i] = x[i] + alpha*p[i]; /*���� r[N]�ĸ��� */    
		r[i] = r[i]-alpha*w[i];
	    rho_1 = rho_1 + r[i]*r[i];
	} 
	printf("rho_1:%f\n", rho_1);
	//�����Ĺ���
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
		/*������ x[N]�ĸ���*/    
		rho_1 = 0;
		for(i=0; i < N; i++)   
		{	x[i] = x[i] + alpha * p[i]; /*���� r[N]�ĸ��� */    
			r[i] = r[i] - alpha * w[i];
		    rho_1 = rho_1 + r[i]*r[i];
		} 
		printf("rho_1:%f\n", rho_1);
	}
	free(r);
	free(p);
	free(w);
}
