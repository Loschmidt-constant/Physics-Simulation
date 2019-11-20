#include <stdio.h>
#include <stdlib.h>

#define N 10   /* Nｼ｡ﾀｵﾊｹﾔﾎ・*/

/* ｹﾔﾎﾎﾆﾎﾏ */
void input_matrix( double **a, char c, FILE *fin, FILE *fout);
/* ･ﾙ･ｯ･ﾈ･・ﾎﾆﾎﾏ */
void input_vector( double *b, char c, FILE *fin, FILE *fout);
/* ｹﾔﾎﾎﾎﾎｰ雉ﾎﾊﾝ */
double **dmatrix(int nr1, int nr2, int nl1, int nl2);
/* ｹﾔﾎﾎﾎﾎｰ雋・*/
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
/* ･ﾙ･ｯ･ﾈ･・ﾎｰ隍ﾎｳﾎﾊﾝ */
double *dvector(int i, int j);  
/* ﾎﾎｰ隍ﾎｲ・*/
void free_dvector(double *a, int i); 
/* ｽ､ﾀｵ･ｳ･・ｹ･ｭ｡ｼﾊｬｲ・*/
double **cholesky_decomp( double **a );
/* ｽ､ﾀｵ･ｳ･・ｹ･ｭ｡ｼﾊｬｲﾑ､ｷ､ﾆﾏ｢ﾎｩ1ｼ｡ﾊﾄｰ､ｯ */
double *cholesky_solve( double **a, double *b );

int main(void)
{
  FILE *fin, *fout;
  double **a, *b; 
  int i; 

  /* ｹﾔﾎｪ､隍ﾓ･ﾙ･ｯ･ﾈ･・ﾎﾎﾎｰ雉ﾎﾊﾝ */
  a = dmatrix(1, N, 1, N); /* ｹﾔﾎ・a[1...N][1...N] */
  b = dvector(1,N); /* b[1...N] */

  /* ･ﾕ･｡･､･・ﾎ･ｪ｡ｼ･ﾗ･・*/
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

  /* ｷ・ﾌ､ﾎｽﾐﾎﾏ */
  fprintf( fout, "\nThe solution for Ax = b is ;\n");
  for( i = 1 ; i <= N ; i++)
  {
    fprintf(fout, "%f\t", b[i]);
    fprintf( fout, "\n");
  }

  fclose(fin); fclose(fout);  /* ･ﾕ･｡･､･・ﾎ･ｯ･悅ｼ･ｺ */

  /* ﾎﾎｰ隍ﾎｲ・*/
  free_dmatrix( a, 1, N, 1, N );  free_dvector( b, 1 );

  return 0;
}

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

/* a[1...N][1...N]､ﾎﾆﾎﾏ */
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

/* b[1...N]､ﾎﾆﾎﾏ */
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

double **dmatrix(int nr1, int nr2, int nl1, int nl2)
{
  int i, nrow, ncol; 
  double **a; 

  nrow = nr2 - nr1 + 1 ; /* ｹﾔ､ﾎｿ・*/
  ncol = nl2 - nl1 + 1 ; /* ﾎﾎｿ・*/

  /* ｹﾔ､ﾎｳﾎﾊﾝ */
  if ( ( a=(double **)malloc( nrow*sizeof(double *) ) ) == NULL ) 
  {
    printf("Memory cannot be secured.(from the matrix A)\n");
    exit(1);
  }
  a = a - nr1; /* ｹﾔ､ｺ､鬢ｹ */
  /* ﾎﾎｳﾎﾊﾝ */
  for( i=nr1; i<=nr2; i++) a[i] = (double *)malloc(ncol*sizeof(double)); 
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;             /* ﾎｺ､鬢ｹ */
  
  return(a);
}

void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  /* ･皈筵熙ﾎｲ・ */
  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}

double *dvector(int i, int j) /* a[i]｡ﾁa[i+j]､ﾎﾎﾎｰ隍ﾎﾊﾝ */
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
  free( (void *)(a + i) );  /* (void *)ｷｿ､ﾘ､ﾎ･ｭ･罕ｹ･ﾈ､ｬﾉｬﾍﾗ */
}
