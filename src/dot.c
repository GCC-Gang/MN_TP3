#include "../include/mnblas.h"
#include "../include/complexe2.h"
#include <stdlib.h>
#include <stdio.h>

/*
* 1 FLOP
*/
float mncblas_sdot(const int N, const float *X, const int incX,
                   const float *Y, const int incY)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  float dot = 0.0;

  for (i = 0; i < N; i += incX)
  {
    dot += X[i] * Y[j];
    j += incY;
  }

  return dot;
}

/*
* 1 FLOP
*/
double mncblas_ddot(const int N, const double *X, const int incX,
                    const double *Y, const int incY)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  double dot = 0.0;

  for (i = 0; i < N; i += incX)
  {
    dot += X[i] * Y[j];
    j += incY;
  }

  return dot;
}

/*
* 8 FLOP
*/
void mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  register complexe_float_t dot;
  register complexe_float_t mult;
  for (i = 0; i < N; i += incX)
  {
    mult = mult_complexe_float(((complexe_float_t *)X)[i], ((complexe_float_t *)Y)[j]); // 6 FLOP
    dot = add_complexe_float(dot, mult);  // 2 FLOP
    j += incY;
  }
  *(complexe_float_t *)dotu = dot;
}

/*
* 9 FLOP
*/
void mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  register complexe_float_t dot;
  register complexe_float_t mult, conjuge_X;
  for (i = 0; i < N; i += incX)
  {
    conjuge_X.real = ((complexe_float_t *)X)[i].real;
    conjuge_X.imaginary = -((complexe_float_t *)X)[i].imaginary;
    mult = mult_complexe_float(conjuge_X, ((complexe_float_t *)Y)[j]);
    dot = add_complexe_float(dot, mult);
    j += incY;
  }
  *(complexe_float_t *)dotc = dot;

  return;
}

/*
* 8 FLOP
*/
void mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  register complexe_double_t dot;
  register complexe_double_t mult;
  for (i = 0; i < N; i += incX)
  {
    mult = mult_complexe_double(((complexe_double_t *)X)[i], ((complexe_double_t *)Y)[j]);
    dot = add_complexe_double(dot, mult);
    j += incY;
  }
  *(complexe_double_t *)dotu = dot;
}

/*
* 9 FLOP
*/
void mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0;
  register unsigned int j = 0;
  register complexe_double_t dot;
  register complexe_double_t mult, conjuge_X;
  for (i = 0; i < N; i += incX)
  {
    conjuge_X.real = ((complexe_double_t *)X)[i].real;
    conjuge_X.imaginary = -((complexe_double_t *)X)[i].imaginary;
    mult = mult_complexe_double(conjuge_X, ((complexe_double_t *)Y)[j]);
    dot = add_complexe_double(dot, mult);
    j += incY;
  }
  *(complexe_double_t *)dotc = dot;

  return;
}
