#include <stdio.h>
#include <x86intrin.h>

#include "../include/mnblas.h"
#include "../include/complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef complexe_float_t vcomplexe[VECSIZE];

vfloat vec1;
vcomplexe vec2;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vectorC_init (vcomplexe V, float a, float b)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++) {
      V[i].real = a;
      V[i].imaginary = b;
  }

  return ;
}

void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}

void vectorC_print (vcomplexe V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f %f, ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;

  return ;
}

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 register int i ;
 float res = 0.f;

 init_flop () ;

 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;

     start = _rdtsc () ;
        res = mnblas_sasum (VECSIZE, vec1, 1) ;
     end = _rdtsc () ;

     printf ("mnblas_sasum %d : resultat = %f, nombre de cycles: %Ld \n", i, res, end-start) ;
     calcul_flop ("sasum ", VECSIZE, end-start) ;
   }

   printf("=========================================\n");

   for (i = 0; i < NB_FOIS; i++) {
       vectorC_init(vec2, 1.0, 2.0);

       start = _rdtsc () ;
        res = mnblas_scasum (VECSIZE, vec2, 1) ;
       end = _rdtsc () ;

       printf ("mnblas_scasum %d : resultat = %f, nombre de cycles: %Ld \n", i, res, end-start) ;
       calcul_flop ("scasum ", VECSIZE * 2, end-start) ;
   }

}
