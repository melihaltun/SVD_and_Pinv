/**
* @file eigen.cpp
* @author Melih Altun @2015
**/

#include "eigen.h"

//calculates eigenvalues for an input matrix. eigenvectors are also computed if the matrix is symmetric
//parameters: (outputs) eigenvector matrix as column vectors, eigenvalues, (inputs) input matrix - must be square, size of matrix, tolerance for precision. 
void eig_symetric(float V[], float D[], float in[], int size, float tolr) {
	int i, j, k;
	float norm1, residue;

	float *diagVec;
	diagVec = new float[size];
	float *M;
	M = new float[size*size];
	float *Q;
	Q = new float[size*size];
	float *R;
	R = new float[size*size];
	float *Vnew;
	Vnew = new float[size*size];

	copy_matrix(M, in,size,size);

	if (tolr == 0)
		tolr = DEFAULT_TOLR;

	memset(V, 0, size*size*sizeof(float));

#if defined (CHECK_SYMMETRY)
	bool symmetric = true;
	bool init = true;

	for (i = 0; i < size; i++) {
		for (j = i+1; j < size; j++){
			if (M[lin_index(i, j, size)] != M[lin_index(j, i, size)]) {
				symmetric = false;
				break;
			}
		}
		if (!symmetric)
			break;
	}
#endif

	for (k = 0; k < MAX_ITER; k++) {
		diagonal_to_vector(diagVec, M, size);
		norm1 = vector_norm(diagVec, size);
		QR(Q, R, M, size, size);
#if defined(CHECK_SYMMETRY)
		if (symmetric) {
			if (init) {
				copy_matrix(V, Q, size, size);
				init = false;
			} else {
				multiply_square_matrices(Vnew, V, Q, size);
				copy_matrix(V, Vnew, size, size);
			}
		}
#endif
		multiply_square_matrices(M, R, Q, size);
		diagonal_to_vector(diagVec, M, size);
		residue = fabs(norm1 - vector_norm(diagVec, size));

		if (residue < tolr)
			break;
	}
	diagonal_to_vector(D, M, size);

	//TODO: Add Hessenberg reduction for non-symmetric matrices for later use.

	delete[] diagVec;
	diagVec = NULL;
	delete[] M;
	M = NULL;
	delete[] Q;
	Q = NULL;
	delete[] R;
	R = NULL;
	delete[] Vnew;
	Vnew = NULL;
}