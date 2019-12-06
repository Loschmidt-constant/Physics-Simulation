#include <stdio.h>
#include <stdlib.h>

/* 行列成分の宣言 */
enum {
	alpha = -2,	//対角成分
	N = 10			//サイズ
};

int main(void)
{	
	/* 変数と関数の宣言 */
	int ix, jy, p, q, k;
	int a[N+1][N+1];
	int b[N+1];
	void print_matrix(int, int matrix[][N+1]);

	/*****************************/
	/* 	pのループで作る場合		    */
	/*****************************/
	/* (1) 初めは，全部 0 を入れておく．*/
	for (p=1; p<=N; p++){
		for (q=1; q<=N; q++){
			a[p][q] = 0;
		}
	}

	/* (2) 0でない成分を入れる．*/
	for (p=1; p<=N; p++){
		ix=(p-1)/N+1;
		jy=(p-1)%N+1;

		/* 対角ブロック */
		a[p][p] = alpha;
		if (jy >= 2)      a[p][p-1] = 1;
		if (jy <= N-1) a[p][p+1] = 1;

	}
	
	/********************/
	/* 		ベクトル		    */
	/*******************/
	for (k=1; k<=N; k++){
		b[k] = 0;
	}
	b[N] = -(N+1);
	

	/*ファイル出力*/
	 FILE *outputfile;         // 出力ストリーム
	outputfile = fopen("input_matrix.dat", "w");  // ファイルを書き込み用にオープン(開く)
	if (outputfile == NULL) {          // オープンに失敗した場合
		printf("cannot open\n");         // エラーメッセージを出して
		exit(1);                         // 異常終了
	}	
	
	fprintf(outputfile,"",N);		
	print_matrix(N, a);
	fprintf(outputfile,"\n");		
	
	for (k=1; k<=N; k++) {
			fprintf(outputfile,"%4d", b[k]);
		}

	fclose(outputfile);          // ファイルをクローズ(閉じる)
	
	return 0;
}


void print_matrix(int N, int matrix[][N+1])
{
	FILE *outputfile;         // 出力ストリーム
	
	int p, q;

	for (p=1; p<=N; p++) {
		for (q=1; q<=N; q++) {
			fprintf(outputfile,"%4d", matrix[p][q]);
		}
		fprintf(outputfile,"\n");		
	}
}
