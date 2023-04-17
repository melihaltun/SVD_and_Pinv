/**
* @file svd.cpp
* @author Melih Altun @2015
**/
#include "svd.h"

/*Singular Value Decomposition
Parameters: (outputs) 1st eigen vector set, diagonal singular values, 2nd eigenvector set, (inputs) input matrix, row count, col count */
void svd(float UU[], float S[], float VV[], float X[], int N, int M)
{
	int i, j;
	bool transposed = false;
	float *V = VV, *U = UU;
	if (N < M) { //transpose if rows are less than columns
		int tmp = N;
		float *Xnew;
		Xnew = new float[M*N];
		transpose(Xnew, X, N, M);
		N = M;
		M = tmp;
		copy_matrix(X, Xnew, N, M);
		V = UU;
		U = VV;
		transposed = true;
		delete[] Xnew;
	}

	float *Xtr, *Xtr_X, *D, *V0, *U0;
	Xtr_X = new float[M*M];
	Xtr = new float[N*M];
	D = new float[M];
	U0 = new float[N];
	V0 = new float[M];


	transpose(Xtr, X, N, M);  //X'
	multiply_matrices(Xtr_X, Xtr, X, M, N, M);   //X'*X

	eig_symetric(V, D, Xtr_X, M, 0); //[V,D]=eig(X'*X)
	
	for (i = 0; i < M; i++) {
		if (D[i] >= 0)
			D[i] = sqrt(D[i]);
		else
			D[i] = sqrt(-D[i]);  // Something is wrong if X'*X produces negative eigenvalues.
	}

	vector_to_diagonal(S, D, M);

	for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++)
			V0[j] = V[lin_index(j, i, M)];
		multiply_matrix_with_vector(U0, X, V0, N, M);
		for (j = 0; j < N; j++)
			U[lin_index(j, i, M)] = U0[j] / D[i];
	}

	if (transposed)
	{
		float *Xnew;
		Xnew = new float[M*N];
		transpose(Xnew, X, N, M);
		copy_matrix(X, Xnew, M, N);
		delete[] Xnew;
	}

	delete[] Xtr;
	Xtr = NULL;
	delete[] Xtr_X;
	Xtr_X = NULL;
	delete[] D;
	D = NULL;
	delete[] V0;
	V0 = NULL;
	delete[] U0;
	U0 = NULL;
}

/* Computes Pseudo Inverse of a matrix using SVD
Parameters: (output) inverted matrix, (inputs) matrix, row count, col count*/
void psInv(float Y[], float X[], int N, int M)
{
	int i;

	float *U, *V, *S, *Sinv, *Utr, *V_Sinv;
	bool transposed = false;
	if (N < M) { //transpose if rows are less than columns
		int tmp = N;
		float* Xnew;
		Xnew = new float[M * N];
		transpose(Xnew, X, N, M);
		N = M;
		M = tmp;
		copy_matrix(X, Xnew, N, M);
		transposed = true;
		delete[] Xnew;
	}

	S = new float[M*M];
	Sinv = new float[M*M];
	U = new float[N*M];
	V = new float[M*M];
	Utr = new float[N*M];
	V_Sinv = new float[M*M];

	svd(U, S, V, X, N, M);

	memset(Sinv, 0, M*M*sizeof(float));

	for (i = 0; i < M; i++)
		Sinv[lin_index(i, i, M)] = 1 / S[lin_index(i, i, M)];

	multiply_square_matrices(V_Sinv, V, Sinv, M);

	transpose(Utr, U, N, M);

	multiply_matrices(Y, V_Sinv, Utr, M, M, N);

	if (transposed) {
		float* Xnew;
		Xnew = new float[M * N];
		transpose(Xnew, X, N, M);
		copy_matrix(X, Xnew, N, M);
		transpose(Xnew, Y, M, N );
		copy_matrix(Y, Xnew, M, N);
		delete[] Xnew;
	}

	delete[] S;
	S = NULL;
	delete[] Sinv;
	Sinv = NULL;
	delete[] U;
	U = NULL;
	delete[] V;
	V = NULL;
	delete[] Utr;
	Utr = NULL;
	delete[] V_Sinv;
	V_Sinv = NULL;
}
