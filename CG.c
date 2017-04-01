#include<stdlib.h>
void CGmethod(double *A, int N,double x[],double b[])
{ /***** ע�� x[]��ʾ��ʼ��ֵ  
  **** n ��ʾ�����ά����
  ****b[]��ʾ�Ҷ���  
  *****p[]Ϊ����ľ��󣨺��������� P[]����������������
  *****r[]��ʾΪ�в����� */ 
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
	
   /*��������*/ 
	for(k=0; k<N; k++)  
	{	 s1 = s2 = 0;  
		for(i = 0; i < N; i++) 
			s1=s1 + r[i] * p[i];    
		for(i=0; i<N; i++)    
			for(j=0; j<N; j++)   
				s2 = s2 + p[i] * A[i*N+j] * p[j];  
		alpha = s1/s2;     
		s1 = s2 = 0;
	/*������ x[N]�ĸ���*/    
		for(i=0; i<N; i++)   
			x[i] = x[i] + alpha*p[i]; /*���� r[N]�ĸ��� */   
		for(i=0; i<N; i++)       
			for(j=0; j<N; j++)      
				r[i] = r[i]-alpha*A[i*N+j]*p[j];
	 /*���в�=0��˵�����ҵ���ȷ�⣬�˳�ѭ��*/ 
	for(i=0;i<N;i++)     
		 for(j=0;j<N;j++)      
			   s2 = s2 + p[i]*A[i*N+j]*p[j];    
	for(i=0; i<N; i++)    
		  for(j=0; j<N; j++)       
			   s1 = s1-p[i]*A[i*N+j]*r[j];    
		 beta = s1/s2;
	 //����ĸ���
	 for(i=0; i<N; i++)
		 p[i] = r[i] + beta*p[i];     
	}
	free(r);
	free(p);
}
