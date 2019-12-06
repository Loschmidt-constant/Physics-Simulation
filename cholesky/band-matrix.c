#include <stdio.h>
#include <stdlib.h>

/* �s�񐬕��̐錾 */
enum {
	alpha = -2,	//�Ίp����
	N = 10			//�T�C�Y
};

int main(void)
{	
	/* �ϐ��Ɗ֐��̐錾 */
	int ix, jy, p, q, k;
	int a[N+1][N+1];
	int b[N+1];
	void print_matrix(int, int matrix[][N+1]);

	/*****************************/
	/* 	p�̃��[�v�ō��ꍇ		    */
	/*****************************/
	/* (1) ���߂́C�S�� 0 �����Ă����D*/
	for (p=1; p<=N; p++){
		for (q=1; q<=N; q++){
			a[p][q] = 0;
		}
	}

	/* (2) 0�łȂ�����������D*/
	for (p=1; p<=N; p++){
		ix=(p-1)/N+1;
		jy=(p-1)%N+1;

		/* �Ίp�u���b�N */
		a[p][p] = alpha;
		if (jy >= 2)      a[p][p-1] = 1;
		if (jy <= N-1) a[p][p+1] = 1;

	}
	
	/********************/
	/* 		�x�N�g��		    */
	/*******************/
	for (k=1; k<=N; k++){
		b[k] = 0;
	}
	b[N] = -(N+1);
	

	/*�t�@�C���o��*/
	 FILE *outputfile;         // �o�̓X�g���[��
	outputfile = fopen("input_matrix.dat", "w");  // �t�@�C�����������ݗp�ɃI�[�v��(�J��)
	if (outputfile == NULL) {          // �I�[�v���Ɏ��s�����ꍇ
		printf("cannot open\n");         // �G���[���b�Z�[�W���o����
		exit(1);                         // �ُ�I��
	}	
	
	fprintf(outputfile,"",N);		
	print_matrix(N, a);
	fprintf(outputfile,"\n");		
	
	for (k=1; k<=N; k++) {
			fprintf(outputfile,"%4d", b[k]);
		}

	fclose(outputfile);          // �t�@�C�����N���[�Y(����)
	
	return 0;
}


void print_matrix(int N, int matrix[][N+1])
{
	FILE *outputfile;         // �o�̓X�g���[��
	
	int p, q;

	for (p=1; p<=N; p++) {
		for (q=1; q<=N; q++) {
			fprintf(outputfile,"%4d", matrix[p][q]);
		}
		fprintf(outputfile,"\n");		
	}
}
