#include<stdio.h>
#include<stdlib.h> 
#include<math.h> 

#define N 10 /*��̗v�f��*/
#define EPS pow(10.0,-8.0)
#define KMAX 500

void input_matrix(double **a, char c, FILE *fin, FILE *fout);
void input_vector(double *b, char c, FILE *fin, FILE *fout);

double **dmatrix(int nr1, int nr2, int nl1, int nl2);
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
double *dvector(int i, int j);
void free_dvector(double *a, int i);

int double_comp(const void *s1, const void *s2); /*��r�֐��i�����j*/ 
double vector_norm_max(double *a, int m, int n);//�ő�l�m�����̌v�Z
double *gauss_seidel(double **a, double *b, double *x);//Gauss-Seidel�@

int main(void)
{
	FILE *fin, *fout;
	double **a, *b, *x;
	int i;

	/*�����s��, �x�N�g���̗̈�m��*/
	a=dmatrix(1,N,1,N);
	b=dvector(1,N);
	x=dvector(1,N);

	if((fin=fopen("input_matrix.dat","r"))==NULL)
	{
		printf("�t�@�C����������܂���Finput_matrix.dat \n");
		exit(1);
	}

	if((fout=fopen("output_gs.dat","w"))==NULL)
	{
		printf("�t�@�C�����쐬�ł��܂���Foutput_gs.dat \n");
		exit(1);
	}

	input_matrix(a, 'A', fin, fout); /*�s��A�̓��o��*/
	input_vector(b, 'b', fin, fout); /*�x�N�g��b�̓��o��*/
	input_vector(x, 'x', fin, fout); /*�x�N�g��b�̓��o��*/
	x=gauss_seidel(a,b,x);
	
	/*���ʂ̏o��*/
	fprintf(fout, "Ax=b�̉��͎��̒ʂ�ł�\n");
	for(i=1;i<=N;i++)
	{
		fprintf(fout, "%20.16f\n", x[i]);
	}

	fclose(fin); fclose(fout);

	free_dmatrix(a,1,N,1,N);
	free_dvector(b,1);
	free_dvector(x,1);

  return 0;
}

void input_matrix(double **a, char c, FILE *fin, FILE *fout)
{
	int i,j;

	fprintf(fout, "�s��%c�͎��̒ʂ�ł�\n", c);
	for(i=1;i<=N;i++)
	{
		for(j=1;j<=N;j++)
		{
			fscanf(fin,"%lf", &a[i][j]); //���͎���%g�͎g���Ȃ�
			fprintf(fout, "%f\t",a[i][j]); //�S��6�P�^�A�����_�ȉ�2��
		}
		fprintf(fout,"\n");
	}

}

void input_vector(double *b, char c, FILE *fin, FILE *fout)
{
	int i;

	fprintf(fout, "�x�N�g��%c�͎��̒ʂ�ł�\n", c);
	for(i=1;i<=N;i++)
	{
		fscanf(fin,"%lf", &b[i]);
		fprintf(fout, "%6.2f\t",b[i]);
		fprintf(fout,"\n");
	}

}

double **dmatrix(int nr1, int nr2, int nl1, int nl2)
{
	int i, nrow, ncol;
	double **a;

	nrow=nr2-nr1+1; /*�s��*/
	ncol=nl2-nl1+1; /*��*/

	/*�s�̊m��*/
	if((a=malloc(nrow*sizeof(double *))) == NULL)
	{
		printf("���������m�ۂł��܂���(from dvector)\n");
		exit(1);
	}
	a=a-nr1; /*�s�����炷*/

	/*��̊m��*/
	for(i=nr1; i<=nr2; i++) a[i]=malloc(ncol*sizeof(double));

	/*������炷*/
	for(i=nr1; i<=nr2; i++) a[i]=a[i]-nl1;

	return(a);
}

void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2)
{
	int i;

	for(i=nr1; i<=nr2; i++) free((void *)(a[i]+nl1));
	free((void *)(a+nr1));
}

double *dvector(int i, int j)
{
	double *a;

	if((a=malloc((j-i+1)*sizeof(double))) == NULL)
	{
		printf("���������m�ۂł��܂���(from dvector)\n");
		exit(1);
	}

	return(a-i);
}

void free_dvector(double *a, int i)
{
	free((void *)(a+i));
}

int double_comp(const void *s1, const void *s2) /*��r�֐��i�����j*/ 
{
	const double a1=*((double *)s1); /* (double *)�փL���X�g */
	const double a2=*((double *)s2); /* (double *)�փL���X�g */
	
	if(a1<a2)
	{
		return -1;
	}
	else if(a1==a2)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

double vector_norm_max(double *a, int m, int n)//�ő�l�m�����̌v�Z
{
	int i,tmp;
	tmp=n-m+1; /* �S�v�f�� */
	for(i=m;i<=n;i++)a[i]=fabs(a[i]);
	/*���בւ�*/
	qsort(a+m,tmp,sizeof(a[0]),double_comp);
	return a[n];
}

double *gauss_seidel(double **a, double *b, double *x)//�K�E�X�E�U�C�f���@
{
	double eps, *xo, s, t;
	int i,j,k=0;
	
	xo=dvector(1,N);
	
	do
	{
		for(i=1;i<=N;i++) xo[i]=x[i];
		
		/* i=1�̏���*/
		t=0.0;
		for(j=2;j<=N;j++) t+=a[1][j]*xo[j];
		x[1]=(b[1]-t)/a[1][1];
		
		/* i=2...N�̏��� */
		for(i=2;i<=N;i++)
		{
			s=0.0; t=0.0;
			for(j=1;j<i; j++) s+=a[i][j]*x[j];
			for(j=i+1;j<=N; j++) t+=a[i][j]*xo[j];
			x[i]=(b[i]-s-t)/a[i][i];
			
		}
		
		/* �������� */
		for(i=1;i<=N;i++) xo[i]=xo[i]-x[i];
		eps=vector_norm_max(xo,1,N);
		
		k++;
		printf("%d %10.8f \n",k,xo[1]);
	}while(eps>EPS && k<KMAX);
	
	free_dvector(xo,1);
	if(k==KMAX)
	{
		printf("������������܂���ł���\n");
		exit(1);
	}
	else
	{
		printf("Number of iterations : %d\n",k);
		return x;
	}
}
