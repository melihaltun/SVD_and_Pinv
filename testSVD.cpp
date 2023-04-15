/**
* @file testSVD.cpp
* @author Melih Altun @2023
**/

#include "svd.h"


int main()
{
	float A[] = { 1, 2,
				  3, 0,
				  4, -3 };
	float Ainv[3 * 2], U[6], S[4], V[4];

	svd(U, S, V, A, 3, 2);
	psInv(Ainv, A, 3, 2);

	return 0;
}
