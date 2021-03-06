#include "../include/mnblas.h"
#include "../include/complexe.h"
#include <math.h>

void mncblas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY) 
{
  unsigned int i;
  
  //on part du postulat que incX = incY
  register unsigned int maxIter = N/incX; 
  


  #pragma omp parallel for
  
  for (i=0; i < maxIter ; i += incX)
    {
        Y [i] = X [i] ;      
    }

  return ;
}

void mncblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY)
{
    register unsigned int i = 0 ;
    //on part du postulat que incX = incY

    register unsigned int maxIter = N/incX; 
    #pragma omp parallel for

     for (i=0; i < maxIter ; i += incX)
    {
        Y [i] = X [i] ;      
    }

  return ;
}

void mncblas_ccopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;

    register unsigned int maxIter = N/incX; 
    #pragma omp parallel for

     for (i=0; i < maxIter ; i += incX)
      {
        ((complexe_float_t*) Y) [i].real = ((complexe_float_t*) X) [i].real ;
        ((complexe_float_t*) Y) [i].imaginary = ((complexe_float_t*) X) [i].imaginary ;
      }

    return ;
}

void mncblas_zcopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;
   register unsigned int maxIter = N/incX; 
    #pragma omp parallel for

     for (i = 0; i < maxIter ; i += incX)
      {
          ((complexe_double_t*) Y) [i].real = ((complexe_double_t*) X) [i].real ;
          ((complexe_double_t*) Y) [i].imaginary = ((complexe_double_t*) X) [i].imaginary ;
       }

    return ;
}
