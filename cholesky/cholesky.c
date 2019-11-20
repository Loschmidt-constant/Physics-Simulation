#include <stdio.h>
#include <stdlib.h>

#define N 3   /* N���������� */

/* ��������� */
void input_matrix( double **a, char c, FILE *fin, FILE *fout);
/* �٥��ȥ������ */
void input_vector( double *b, char c, FILE *fin, FILE *fout);
/* ������ΰ���� */
double **dmatrix(int nr1, int nr2, int nl1, int nl2);
/* ������ΰ���� */
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
/* �٥��ȥ��ΰ�γ��� */
double *dvector(int i, int j);  
/* �ΰ�β��� */
void free_dvector(double *a, int i); 
/* �������쥹����ʬ�� */
double **cholesky_decomp( double **a );
/* �������쥹����ʬ������Ѥ���ϢΩ1����������� */
double *cholesky_solve( double **a, double *b );

int main(void)
{
  FILE *fin, *fout;
  double **a, *b; 
  int i; 

  /* ���󤪤�ӥ٥��ȥ���ΰ���� */
  a = dmatrix(1, N, 1, N); /* ���� a[1...N][1...N] */
  b = dvector(1,N); /* b[1...N] */

  /* �ե�����Υ����ץ� */
  if ( (fin = fopen( "input_cho.dat", "r")) == NULL )
  {
      printf("File not found : input_cho.dat \n");
      exit(1);
  }
  if( (fout = fopen( "output_cho.dat", "w")) == NULL )
  {
    printf("Unable to create file : output_cho.dat \n");
    exit(1);
   }

  input_matrix( a, 'A', fin, fout );  /* ����A�������� */  
  input_vector( b, 'b', fin, fout );  /* �٥��ȥ�b�������� */
  a = cholesky_decomp( a );           /* �������쥹����ʬ�� */
  b = cholesky_solve( a, b );         /* ������������������ */

  /* ��̤ν��� */
  fprintf( fout, "\nThe solution for Ax = b is ;\n");
  for( i = 1 ; i <= N ; i++)
  {
    fprintf(fout, "%f\t", b[i]);
    fprintf( fout, "\n");
  }

  fclose(fin); fclose(fout);  /* �ե�����Υ����� */

  /* �ΰ�β��� */
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

/* a[1...N][1...N]������ */
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

/* b[1...N]������ */
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

  nrow = nr2 - nr1 + 1 ; /* �Ԥο� */
  ncol = nl2 - nl1 + 1 ; /* ��ο� */

  /* �Ԥγ��� */
  if ( ( a=(double **)malloc( nrow*sizeof(double *) ) ) == NULL ) 
  {
    printf("Memory cannot be secured.(from the matrix A)\n");
    exit(1);
  }
  a = a - nr1; /* �Ԥ򤺤餹 */
  /* ��γ��� */
  for( i=nr1; i<=nr2; i++) a[i] = (double *)malloc(ncol*sizeof(double)); 
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;             /* ��򤺤餹 */
  
  return(a);
}

void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  /* ����β���  */
  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}

double *dvector(int i, int j) /* a[i]��a[i+j]���ΰ����� */
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
  free( (void *)(a + i) );  /* (void *)���ؤΥ��㥹�Ȥ�ɬ�� */
}
