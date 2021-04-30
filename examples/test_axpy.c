#include "../include/mnblas.h"
#include "../include/complexe2.h"
#include "flop.h"

#include <stdio.h>
#include <x86intrin.h>

#define NB_FOIS 10000
#define VECSIZE 65536

typedef float vfloat[VECSIZE];
typedef complexe_float_t vComplexefloat[VECSIZE];

vfloat vec1, vec2;
vComplexefloat vecComplexef1, vecComplexef2;

//initialisation du vecteur
void vector_init(vfloat V, float x)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        V[i] = x;

    return;
}
void vector_init_complexe(vComplexefloat V, float reel, float imaginaire)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
    {
        V[i].real = reel;
        V[i].imaginary = imaginaire;
    }

    return;
}

void vector_print(vfloat V)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("%f ", V[i]);
    printf("\n");

    return;
}

int main(int argc, char **argv)
{

    //Test float
    unsigned long long int start, end ;
    int i;
    float sum_saxpy = 0;
    float sum_caxpy = 0;
    printf("--------------TEST mnblas_saxpy-------------\n");

    init_flop();
    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init(vec1, 1.0);
        vector_init(vec2, 2.0);
        start = _rdtsc();
        mnblas_saxpy(VECSIZE, 3, vec1, 1, vec2, 1);
        end = _rdtsc();

        printf("mncblas_saxpy %d : res = %3.2f nombre de cycles: %Ld \n", i, vec2[0], end - start);
        sum_saxpy += calcul_flop_return("saxpy ", 2 * VECSIZE, end - start);
    }



    //Test avec complexes float
    printf("--------------TEST mnblas_caxpy-------------\n");

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init_complexe(vecComplexef1, 1.0, 0.0);
        vector_init_complexe(vecComplexef2, 2.0, 0.0);
        start = _rdtsc();
        complexe_float_t alpha = {3.0, 0.0};

        mnblas_caxpy(VECSIZE, &alpha, vecComplexef1, 1, vecComplexef2, 1);
        end = _rdtsc();

        printf("mncblas_caxpy %d :  nombre de cycles: %Ld \n", i, end - start);
        sum_caxpy += calcul_flop_return("caxpy ", 8 * VECSIZE, end - start);
    }


    printf("perf moyenne saxpy : %f GFLOP/s\n", sum_saxpy / NB_FOIS);
    printf("perf moyenne caxpy : %f GFLOP/s \n", sum_caxpy / NB_FOIS);

    return 0;
}
