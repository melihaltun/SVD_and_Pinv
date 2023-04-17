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
	float Ainv[2 * 3], U_A[3 * 2], S_A[2 * 2], V_A[2 * 2];

	svd(U_A, S_A, V_A, A, 3, 2);
	psInv(Ainv, A, 3, 2);

	/*
	A_inv = [ 0.1387  0.1639  0.0924
			  0.2605  0.1261 -0.1597]
	*/

	float B[] = { 1, 3, 4,
				  2, 0,-3 };
	float Binv[3 * 2], U_B[2 * 2], S_B[2 * 2], V_B[3 * 2 ];
	svd(U_B, S_B, V_B, B, 2, 3);
	psInv(Binv, B, 2, 3);

	/*
	B_inv  = [0.1387  0.2605
			  0.1639  0.1261
			  0.0924 -0.1597 ]
	*/

	return 0;
}
