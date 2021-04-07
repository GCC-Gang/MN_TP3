#include "../include/mnblas.h"
#include "../include/complexe.h"
#include "flop.h"

#include <stdio.h>
#include <x86intrin.h>
#include <math.h>

#define NB_FOIS 41943
#define VECSIZE 65536

typedef double vdouble[VECSIZE];
typedef complexe_double_t vComplexedouble[VECSIZE];

vdouble vec1, vec2;
vComplexedouble vecComplexef1;

void vector_init(vdouble V, float x)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

void vector_init_complexe(vComplexedouble V, float reel, float imaginaire)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
  {
    V[i].real = reel;
    V[i].imaginary = imaginaire;
  }

  return;
}

void vector_print(vdouble V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%f ", V[i]);
  printf("\n");

  return;
}

int main(int argc, char **argv)
{
    unsigned long long start, end;
    float res;
    int i;

    init_flop();
    double somme = 0.f;
    for (i = 0; i < NB_FOIS; i++)
    {
    vector_init(vec1, 1.0);
    res = 0.0;

    start = _rdtsc();
    res =  mnblas_dnrm2(VECSIZE, vec1, 1);
    end = _rdtsc();

    printf("mnblas_dnrm2 %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end - start);
    somme += calcul_flop_return("mnblas_dnrm2", 2 * VECSIZE + 1, end - start);
    }
    printf("\n");
    //complexe_float_t resComplexe ;

    double somme2 = 0.f;
    for (i = 0; i < NB_FOIS; i++)
  {
      vecComplexef1->real = 1.0;
      vecComplexef1->imaginary = 2.0;
    vector_init_complexe(vecComplexef1, 1.0, 2.0);
    res = 0.0;

    start = _rdtsc();
    mnblas_dznrm2(VECSIZE, vecComplexef1,1);
    end = _rdtsc();


    printf("mnblas_dznrm2 %d : res = %3.2f +i %3.2f nombre de cycles: %Ld \n", i, vecComplexef1->real, vecComplexef1->imaginary, end - start);

    somme2 += calcul_flop_return("mnblas_dznrm2 ", 8 * VECSIZE, end - start);
  }

printf("\n=========== Résultats nrm2 double précision ===========\n");
printf("En moyenne on a une performance de %5.3f GFLOP/s\n", somme/NB_FOIS);

printf("=========== Résultats nrm2 complexe double précision ===========\n");
printf("En moyenne on a une performance de %5.3f GFLOP/s\n", somme2/NB_FOIS);





}
