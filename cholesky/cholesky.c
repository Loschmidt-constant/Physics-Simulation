#include <stdio.h>
#include <stdlib.h>

#define N 10   /* N次正方行列 */

/* 行列の入力 */
void input_matrix( double **a, char c, FILE *fin, FILE *fout);
/* ･ベクトルの入力 */
void input_vector( double *b, char c, FILE *fin, FILE *fout);
/* 行列のスペース確保 */
double **dmatrix(int nr1, int nr2, int nl1, int nl2);
/* 行列のスペース開放 */
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
/* ベクトルのスペース確保 */
double *dvector(int i, int j);  
/* ベクトルのスペース開放 */
void free_dvector(double *a, int i); 
/* 修正Cholesky分解法 */
double **cholesky_decomp( double **a );
/* 修正Cholesky分解法による連立方程式の解 */
double *cholesky_solve( double **a, double *b );

int main(void)
{
  FILE *fin, *fout;
  double **a, *b; 
  int i; 

  a = dmatrix(1, N, 1, N); /* 行列 a[1...N][1...N] */
  b = dvector(1,N); /* ベクトル b[1...N] */

  /* ファイル出力 */
  if ( (fin = fopen( "input_matrix.dat", "r")) == NULL )
  {
      printf("File not found : input_matrix.dat \n");
      exit(1);
  }
  if( (fout = fopen( "output_matrix.dat", "w")) == NULL )
  {
    printf("Unable to create file : output_matrix.dat \n");
    exit(1);
   }

  input_matrix( a, 'A', fin, fout );  /* ｹﾔﾎ､ﾎﾆｽﾐﾎﾏ */  
  input_vector( b, 'b', fin, fout );  /* ･ﾙ･ｯ･ﾈ･・､ﾎﾆｽﾐﾎﾏ */
  a = cholesky_decomp( a );           /* ｽ､ﾀｵ･ｳ･・ｹ･ｭ｡ｼﾊｬｲ・*/
  b = cholesky_solve( a, b );         /* ﾁｰｿﾊﾂ衄｡ｦｸ蠶狡衄 */

  /* 結果出力 */
  fprintf( fout, "\nThe solution for Ax = b is ;\n");
  for( i = 1 ; i <= N ; i++)
  {
    fprintf(fout, "%f\t", b[i]);
    fprintf( fout, "\n");
  }

  fclose(fin); fclose(fout);  /* ･ﾕ･｡･､･・ﾎ･ｯ･悅ｼ･ｺ */

  /* スペース開放 */
  free_dmatrix( a, 1, N, 1, N );  free_dvector( b, 1 );

  return 0;
}


/* 修正Cholesky分解法 */
double **cholesky_decomp( double **a )
{
  int i, j, k;
  double tmp;

  for( i = 2; i <= N; i++)
  {
    for( j = 1; j <= i-1; j++)
    {
      tmp = 0.0;
      for ( k = 1; k <= j-1; k++)
      {
        tmp += a[i][k]*a[k][k]*a[j][k];
      }
      a[i][j] = (a[i][j] - tmp) / a[j][j]; 
    }
      tmp = 0.0;
      for ( k = 1; k <= j-1; k++)
      {
        tmp += a[i][k]*a[i][k]*a[k][k];
      }
      a[i][i] = a[i][i] - tmp; 
  }
  return a; 
}


/* 修正Cholesky分解法による連立方程式の解 */
double *cholesky_solve( double **a, double *b )
{
  int i, j;
  double tmp;

  /* LDy = b */
  b[1] = b[1]/a[1][1];
  for( i = 2; i <= N; i++)
  {
    tmp = 0.0;
    for( j = 1; j <= i-1; j++)
    {
      tmp += a[j][j]*a[i][j]*b[j];
    }
    b[i] = ( b[i] - tmp ) / a[i][i];
  }

  /* L^t x = y */
  for( i = N-1; i >= 1; i--)
  {
    tmp = 0.0;
    for( j = i+1; j <= N; j++)
    {
      tmp += a[j][i] * b[j];
    }
    b[i] = b[i] - tmp;
  } 

  return b;
}


/* a[1...N][1...N]の入力 */
void input_matrix( double **a, char c, FILE *fin, FILE *fout)
{
  int i,j;

  fprintf( fout, "The matrix %c is as follows;\n", c);
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

  fprintf( fout, "\nThe vector %c is as follows;\n", c);
  for( i = 1 ; i <= N ; i++)
  {
    fscanf(fin, "%lf", &b[i]);
    fprintf(fout, "%5.2f\t", b[i]);
    fprintf( fout, "\n");
  }
}


/* 行列のスペース確保 */
double **dmatrix(int nr1, int nr2, int nl1, int nl2)
{
  int i, nrow, ncol; 
  double **a; 

  nrow = nr2 - nr1 + 1 ; /* 行数 */
  ncol = nl2 - nl1 + 1 ; /* 列数 */

  /* 行確保 */
  if ( ( a=(double **)malloc( nrow*sizeof(double *) ) ) == NULL ) 
  {
    printf("Memory cannot be secured.(from the matrix A)\n");
    exit(1);
  }
  a = a - nr1; /* 行をずらす */
  /* 列確保 */
  for( i=nr1; i<=nr2; i++) a[i] = (double *)malloc(ncol*sizeof(double)); 
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;             /* ﾎｺ､鬢ｹ */
  
  return(a);
}


/* 行列のスペース開放 */
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  /* メモリ開放 */
  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}

double *dvector(int i, int j) /* a[i] ～ a[i+j]のスペース確保 */
{
  double *a;

  if ( (a=(double *)malloc( ((j-i+1)*sizeof(double))) ) == NULL )
  {
    printf("Memory cannot be secured.(from the vector b) \n");
    exit(1);
  }  

  return(a-i);
}

void free_dvector(double *a, int i)
{
  free( (void *)(a + i) ); 
}
