#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 3   /* N次正方行列 */

/* 行列の入力 */
void input_matrix( double **a, char c, FILE *fin, FILE *fout);
/* ベクトルの入力 */
void input_vector( double *b, char c, FILE *fin, FILE *fout);
/* 行列の領域確保 */
double **dmatrix(int nr1, int nr2, int nl1, int nl2);
/* 行列の領域解放 */
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
/* ベクトル領域の確保 */
double *dvector(int i, int j);  
/* 領域の解放 */
void free_dvector(double *a, int i); 
/* LU分解 */
void lu_decomp( double **a, int *p);
/* LU分解を利用して連立1次方程式を解く */
double *lu_solve( double **a, double *b, int *p );

int main(void)
{
  FILE *fin, *fout;
  double **a, *b; 
  int i, p[N]; /* p[1...N-1]を利用, p[0]は未使用 */

  /* 行列およびベクトルの領域確保 */
  a = dmatrix(1, N, 1, N); /* 行列 a[1...N][1...N] */
  b = dvector(1,N); /* b[1...N] */

  /* ファイルのオープン */
  if ( (fin = fopen( "input_matrix.dat", "r")) == NULL )
  {
      printf("ファイルが見つかりません : input_matrix.dat \n");
      exit(1);
  }
  if( (fout = fopen( "output_LU.dat", "w")) == NULL )
  {
    printf("ファイルが作成できません : output_LU.dat \n");
    exit(1);
   }

  input_matrix( a, 'A', fin, fout );    /* 行列Aの入出力 */  
  input_vector( b, 'b', fin, fout );    /* ベクトルbの入出力 */
  lu_decomp( a, p);                 /* LU分解 */
  b = lu_solve( a, b, p );         /* 前進代入・後退代入 */

  /* 結果の出力 */
  fprintf( fout, "Ax=bの解は次の通りです\n");
  for( i = 1 ; i <= N ; i++)
  {
    fprintf(fout, "%f\n", b[i]);
  }

  fclose(fin); fclose(fout);  /* ファイルのクローズ */

  /* 領域の解放 */
  free_dmatrix( a, 1, N, 1, N );  free_dvector( b, 1 );

  return 0;
}

void lu_decomp( double **a, int *p )
{
  int i, j, k, ip;
  double alpha, tmp;
  double amax, eps = pow(2.0, -50.0); /* eps = 2^{-50}とする */

  for( k = 1; k <= N-1; k++)
  {
    /* ピボットの選択 */
    amax = fabs(a[k][k]); ip = k;
    for( i = k+1; i <= N; i++)
    {
      if ( fabs(a[i][k]) > amax )
      {
        amax = fabs(a[i][k]); ip = i;
      }
    }
    /* 正則性の判定 */
    if ( amax < eps ) printf("入力した行列は正則ではない!!\n");
    /* ipを配列pに保存 */
    p[k]=ip;
    /* 行交換 */
    if ( ip != k)
    {
      for( j = k; j <= N; j++)
      {
        tmp = a[k][j]; a[k][j]=a[ip][j]; a[ip][j]=tmp;
      }
    }
      /* 前進消去 */
    for( i = k+1; i <= N; i++)
    {
      alpha = - a[i][k]/a[k][k];
      a[i][k] = alpha;
      for( j = k+1; j <= N; j++)
      {
        a[i][j] = a[i][j] + alpha * a[k][j];
      }
    }
  }
}

double *lu_solve( double **a, double *b, int *p )
{
  int i, j, k;
  double tmp;

  /* 右辺の行交換 */
  for( k = 1; k <= N-1; k++)
  {
    tmp=b[k]; b[k]=b[p[k]]; b[p[k]]=tmp;
    /* 前進代入 */
    for( i=k+1; i <= N; i++)
    {
      b[i] = b[i] + a[i][k]*b[k];
    }
  }

  /* 後退代入 */
  b[N] = b[N]/a[N][N];
  for( k = N-1; k >= 1; k--)
  {
    tmp = b[k];
    for( j = k+1; j <= N; j++)
    {
      tmp = tmp - a[k][j] * b[j];
    }
    b[k] = tmp/a[k][k];
  }

  return b;
}

/* a[1...N][1...N]の入力 */
void input_matrix( double **a, char c, FILE *fin, FILE *fout)
{
  int i,j;

  fprintf( fout, "行列%cは次の通りです\n", c);
  for( i = 1 ; i <= N ; i++)
  {
    for( j = 1 ; j <= N ; j++)
    {
      fscanf(fin, "%lf", &a[i][j]);
      fprintf(fout, "%5.2f\t", a[i][j]);
    }
    fprintf( fout, "\n");
  }
}

/* b[1...N]の入力 */
void input_vector( double *b, char c, FILE *fin, FILE *fout)
{
  int i;

  fprintf( fout, "ベクトル%cは次の通りです\n", c);
  for( i = 1 ; i <= N ; i++)
  {
    fscanf(fin, "%lf", &b[i]);
    fprintf(fout, "%5.2f\t", b[i]);
    fprintf( fout, "\n");
  }
}

double **dmatrix(int nr1, int nr2, int nl1, int nl2)
{
  int i, nrow, ncol; 
  double **a; 

  nrow = nr2 - nr1 + 1 ; /* 行の数 */
  ncol = nl2 - nl1 + 1 ; /* 列の数 */

  /* 行の確保 */
  if ( ( a=(double **)malloc( nrow*sizeof(double *) ) ) == NULL ) 
  {
    printf("メモリが確保できません(行列a)\n");
    exit(1);
  }
  a = a - nr1; /* 行をずらす */
  /* 列の確保 */
  for( i=nr1; i<=nr2; i++) a[i] = (double *)malloc(ncol*sizeof(double)); 
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;             /* 列をずらす */
  
  return(a);
}

void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  /* メモリの解放  */
  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}

double *dvector(int i, int j) /* a[i]～a[i+j]の領域を確保 */
{
  double *a;

  if ( (a=(double *)malloc( ((j-i+1)*sizeof(double))) ) == NULL )
  {
    printf("メモリが確保できません(from dvector) \n");
    exit(1);
  }  

  return(a-i);
}

void free_dvector(double *a, int i)
{
  free( (void *)(a + i) );  /* (void *)型へのキャストが必要 */
}