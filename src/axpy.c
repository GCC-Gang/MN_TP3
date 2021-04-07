#include "../include/mnblas.h"
#include "../include/complexe2.h"
#include <stdlib.h>
#include <stdio.h>

/*
* 2 FLOP
*/
void mnblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;

    for (i = 0; i < N; i += incX)
    {
        Y[j] = alpha * X[i] + Y[j];
        j += incY;
    }
}

/*
* 2 FLOP
*/
void mnblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;

    for (i = 0; i < N; i += incX)
    {
        Y[j] = alpha * X[i] + Y[j];
        j += incY;
    }
}

/*
* 8 FLOP
*/
void mnblas_caxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;

    for (i = 0; i < N; i += incX)
    {
        complexe_float_t mult_scalaire = mult_complexe_float(*(complexe_float_t *)alpha, ((complexe_float_t *)X)[i]);
        ((complexe_float_t *)Y)[j] = add_complexe_float(mult_scalaire, ((complexe_float_t *)Y)[j]);
        j += incY;
    }
}

/*
* 8 FLOP
*/
void mnblas_zaxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;

    for (i = 0; i < N; i += incX)
    {
        complexe_double_t mult_scalaire = mult_complexe_double(*(complexe_double_t *)alpha, ((complexe_double_t *)X)[i]);
        ((complexe_double_t *)Y)[j] = add_complexe_double(mult_scalaire, ((complexe_double_t *)Y)[j]);
        j += incY;
    }
}