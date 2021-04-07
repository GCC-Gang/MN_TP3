#include "../include/mnblas.h"
#include "../include/complexe.h"

#include <math.h>
/*

The ?nrm2 routines perform a vector reduction operation defined as
res = ||x||,
where:
x is a vector,
res is a value containing the Euclidean norm of the elements of x.

*/



float  mnblas_snrm2(const int N, const float *X, const int incX){
    register unsigned int i = 0 ;
    register double somme = 0;
    for(; i< N; i+=incX){
        somme += X[i] * X[i];
    }
    return (float) sqrt(somme);


}

double mnblas_dnrm2(const int N, const double *X, const int incX){
    register unsigned int i = 0 ;
    register double somme = 0;
    for(; i< N; i+=incX){
        somme += X[i] * X[i];
    }
    return (double) sqrt(somme);

}

float  mnblas_scnrm2(const int N, const void *X, const int incX){
      register unsigned int i = 0 ;
      const complexe_float_t* new_X = X;

     double somme = 0;
    for(; i< N; i+=incX){

        //somme += new_X[i];
        somme += sqrt(new_X[i].real * new_X[i].real + new_X[i].imaginary * new_X[i].imaginary);

    }
    return (float) sqrt(somme);

}

double mnblas_dznrm2(const int N, const void *X, const int incX){
    register unsigned int i = 0 ;
    const complexe_double_t* new_X = X;

     double somme = 0;
    for(; i< N; i+=incX){

        //somme += new_X[i];
        somme += sqrt(new_X[i].real * new_X[i].real + new_X[i].imaginary * new_X[i].imaginary);

    }
    return (double) sqrt(somme);

}
