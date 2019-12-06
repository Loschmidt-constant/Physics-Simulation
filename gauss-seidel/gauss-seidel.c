#include<stdio.h>
#include<stdlib.h> 
#include<math.h> 

#define N 10 /*列の要素数*/
#define EPS pow(10.0,-8.0)
#define KMAX 500

void input_matrix(double **a, char c, FILE *fin, FILE *fout);
void input_vector(double *b, char c, FILE *fin, FILE *fout);

double **dmatrix(int nr1, int nr2, int nl1, int nl2);
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
double *dvector(int i, int j);
void free_dvector(double *a, int i);

int double_comp(const void *s1, const void *s2); /*比較関数（昇順）*/ 
double vector_norm_max(double *a, int m, int n);//最大値ノルムの計算
double *gauss_seidel(double **a, double *b, double *x);//Gauss-Seidel法

int main(void)
{
	FILE *fin, *fout;
	double **a, *b, *x;
	int i;

	/*正方行列, ベクトルの領域確保*/
	a=dmatrix(1,N,1,N);
	b=dvector(1,N);
	x=dvector(1,N);

	if((fin=fopen("input_matrix.dat","r"))==NULL)
	{
		printf("ファイルが見つかりません：input_matrix.dat \n");
		exit(1);
	}

	if((fout=fopen("output_gs.dat","w"))==NULL)
	{
		printf("ファイルが作成できません：output_gs.dat \n");
		exit(1);
	}

	input_matrix(a, 'A', fin, fout); /*行列Aの入出力*/
	input_vector(b, 'b', fin, fout); /*ベクトルbの入出力*/
	input_vector(x, 'x', fin, fout); /*ベクトルbの入出力*/
	x=gauss_seidel(a,b,x);
	
	/*結果の出力*/
	fprintf(fout, "Ax=bの解は次の通りです\n");
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

	fprintf(fout, "行列%cは次の通りです\n", c);
	for(i=1;i<=N;i++)
	{
		for(j=1;j<=N;j++)
		{
			fscanf(fin,"%lf", &a[i][j]); //入力時に%gは使えない
			fprintf(fout, "%f\t",a[i][j]); //全体6ケタ、小数点以下2桁
		}
		fprintf(fout,"\n");
	}

}

void input_vector(double *b, char c, FILE *fin, FILE *fout)
{
	int i;

	fprintf(fout, "ベクトル%cは次の通りです\n", c);
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

	nrow=nr2-nr1+1; /*行数*/
	ncol=nl2-nl1+1; /*列数*/

	/*行の確保*/
	if((a=malloc(nrow*sizeof(double *))) == NULL)
	{
		printf("メモリが確保できません(from dvector)\n");
		exit(1);
	}
	a=a-nr1; /*行をずらす*/

	/*列の確保*/
	for(i=nr1; i<=nr2; i++) a[i]=malloc(ncol*sizeof(double));

	/*列をずらす*/
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
		printf("メモリが確保できません(from dvector)\n");
		exit(1);
	}

	return(a-i);
}

void free_dvector(double *a, int i)
{
	free((void *)(a+i));
}

int double_comp(const void *s1, const void *s2) /*比較関数（昇順）*/ 
{
	const double a1=*((double *)s1); /* (double *)へキャスト */
	const double a2=*((double *)s2); /* (double *)へキャスト */
	
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

double vector_norm_max(double *a, int m, int n)//最大値ノルムの計算
{
	int i,tmp;
	tmp=n-m+1; /* 全要素数 */
	for(i=m;i<=n;i++)a[i]=fabs(a[i]);
	/*並べ替え*/
	qsort(a+m,tmp,sizeof(a[0]),double_comp);
	return a[n];
}

double *gauss_seidel(double **a, double *b, double *x)//ガウス・ザイデル法
{
	double eps, *xo, s, t;
	int i,j,k=0;
	
	xo=dvector(1,N);
	
	do
	{
		for(i=1;i<=N;i++) xo[i]=x[i];
		
		/* i=1の処理*/
		t=0.0;
		for(j=2;j<=N;j++) t+=a[1][j]*xo[j];
		x[1]=(b[1]-t)/a[1][1];
		
		/* i=2...Nの処理 */
		for(i=2;i<=N;i++)
		{
			s=0.0; t=0.0;
			for(j=1;j<i; j++) s+=a[i][j]*x[j];
			for(j=i+1;j<=N; j++) t+=a[i][j]*xo[j];
			x[i]=(b[i]-s-t)/a[i][i];
			
		}
		
		/* 収束判定 */
		for(i=1;i<=N;i++) xo[i]=xo[i]-x[i];
		eps=vector_norm_max(xo,1,N);
		
		k++;
		printf("%d %10.8f \n",k,xo[1]);
	}while(eps>EPS && k<KMAX);
	
	free_dvector(xo,1);
	if(k==KMAX)
	{
		printf("答えが見つかりませんでした\n");
		exit(1);
	}
	else
	{
		printf("Number of iterations : %d\n",k);
		return x;
	}
}
