/**
* @file svd.h
* @author Melih Altun @2015
**/

#include "eigen.h"

/*Singular Value Decomposition*/
void svd(float U[], float S[], float V[], float X[], int N, int M);

/* Computes Pseudo Inverse of a matrix using SVD*/
void psInv(float Y[], float X[], int N, int M);
