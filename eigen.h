/**
* @file eigen.h
* @author Melih Altun @2015
**/

#include "QR.h"

#define MAX_ITER 1000
#define DEFAULT_TOLR 1E-15
#define CHECK_SYMMETRY

void eig_symetric(float V[], float D[], float in[], int size, float tolr);