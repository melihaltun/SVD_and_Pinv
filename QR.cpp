/**
* @file qr.cpp
* @author Melih Altun @2015
**/

#include "QR.h"

//Givens rotation subroutine
void givensRotation(float *c, float *s, float a, float b)
{
	float r;
	if (b == 0) {
		*c = (float)(1 * sgn(a));
		if (*c == 0)
			*c = 1.0f;
		*s = 0;
	}
	else if (a == 0) {
		*c = 0;
		*s = (float)(-1 * sgn(b));
	}
	else {
		if (fabs(b) > fabs(a)) {
			r = a / b;
			*s = -sgn(b) / sqrtf(1 + r * r);
			*c = -(*s) * r;
		}
		else {
			r = b / a;
			*c = sgn(a) / sqrtf(1 + r * r);
			*s = -(*c) * r;
		}
	}
}

/*QR decomposition using Givens rotation
Parameters: (outputs) NxN Q Matrix, NxM R Matrix, (inputs) NxM input matrix, # of rows, # of cols*/
void QR(float Q[], float R[], float in[], int N, int M)
{
	int i, j;
	float c, s;
	
	float *G, *Gt, *Rnew, *Qnew;
	G = new float[N*N];
	Gt = new float[N*N];
	Rnew = new float[N*M];
	Qnew = new float[N*N];

	set_to_identity(Q, N);
	copy_matrix(R, in, N, M);

	for (j = 0; j < M; j++) {
		for (i = N - 1; i > j; i--) {
			set_to_identity(G, N);
			givensRotation(&c, &s, R[lin_index((i - 1), j, M)], R[lin_index(i, j, M)]);
			G[lin_index((i - 1), (i - 1), N)] = c;
			G[lin_index((i - 1), i, N)] = -s;
			G[lin_index(i, (i - 1), N)] = s;
			G[lin_index(i, i, N)] = c;
			transpose(Gt, G, N, N);
			multiply_matrices(Rnew, G, R, N, N, M);
			multiply_square_matrices(Qnew, Q, Gt, N);
			copy_matrix(R, Rnew, N, M);
			copy_matrix(Q, Qnew, N, N);
		}
	}
	//release dynamic memory
	delete[] G;
	G = NULL;
	delete[] Gt;
	Gt = NULL;
	delete[] Rnew;
	Rnew = NULL;
	delete[] Qnew;
	Qnew = NULL;
}
