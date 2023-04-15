/**
* @file qr.h
* @author Melih Altun @2015
**/

#include "matrixOperations.h"
#include <math.h>

#define sgn(a) ((a) > 0 ? 1 : ((a)< 0 ? -1 : 0))

void QR(float Q[], float R[], float in[], int N, int M);
