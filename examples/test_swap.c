#include "../include/mnblas.h"
#include "../include/complexe.h"
#include "flop.h"

#include <stdio.h>
#include <x86intrin.h>

#define NB_FOIS 41943
#define VECSIZE  500


typedef float vfloat [VECSIZE] ;
vfloat vec1, vec2 ;

typedef double vdouble [VECSIZE];
vdouble vec12, vec22 ;

typedef complexe_float_t vCompFloat [VECSIZE];
vCompFloat vec13, vec23;

typedef complexe_double_t vCompDouble [VECSIZE];
vCompDouble vec14, vec24;


//initialisation du vecteur avec float
void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

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

void vector_init_double (vdouble V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print_double (vdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}
//======
//initialisation du vecteur avec float
void vector_init_float_comp (vCompFloat V, float x, float y)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    V[i].real = x ;
    V[i].imaginary = y;
  }

  return ;
}

void vector_print_float_comp (vCompFloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    printf ("reel: %f; imaginaire: %f ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;
  }

  return ;
}

void vector_init_double_comp (vCompDouble V, double x, double y)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
     V [i].real = x ;
    V [i].imaginary = y;
  }

  return ;
}

void vector_print_double_comp (vCompDouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("reel: %f; imaginaire: %f ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;

  return ;
}


int main(int argc, char **argv)
{

 init_flop () ;

unsigned long long int start, end;
int i;

  //Echange des deux valeurs

    vector_init (vec1, 1.0) ;
    vector_init (vec2, 2.0) ;
    float val_moy = 0;
    float val_moy2 = 0;
    float val_moy3 = 0;
    float val_moy4 = 0;

  for (i = 0; i < NB_FOIS; i++)
  {
    start = _rdtsc () ;


    mncblas_sswap(VECSIZE, vec1, 1, vec2, 1);
    end = _rdtsc () ;
    val_moy += calcul_byte_swap ("mncblas_sswap ", VECSIZE * sizeof(float) * 3, end-start) ; // le *3 car on a besoin de faire 3 affectation mémoire
    printf("mncblas_sswap %d : nombre de cycles: %Ld \n", i, end - start);


  }


//-----------------------------------------------------------------------------------------------------------

    vector_init_double (vec12, 1.f) ;
    vector_init_double (vec22, 2.f) ;

  for (i = 0; i < NB_FOIS; i++)
  {
    start = _rdtsc ();


    mncblas_dswap(VECSIZE, vec12, 1, vec22, 1);
    //printf ("mncblas_sswap %d : nombre de cycles: %lld \n", i, end-start) ;
    end = _rdtsc () ;
    val_moy2 += calcul_byte_swap ("mncblas_dswap ", VECSIZE * sizeof(double) * 3, end-start) ; // le *3 car on a besoin de faire 3 affectation mémoire
    printf("mncblas_dswap %d : nombre de cycles: %Ld \n", i, end - start);


  }

//-----------------------------------------------------------------------------------------------------------
  vector_init_float_comp(vec13, 1.f, 5.f);
  vector_init_float_comp(vec23, 2.f, 2.f);

  for (i = 0; i < NB_FOIS; i++)
  {
    start = _rdtsc () ;


    mncblas_cswap(VECSIZE, vec13, 1, vec23, 1);
    //printf ("mncblas_sswap %d : nombre de cycles: %lld \n", i, end-start) ;
    end = _rdtsc () ;
    val_moy3 += calcul_byte_swap ("mncblas_cswap ", VECSIZE * sizeof(float) * 6, end-start) ; // le *3 car on a besoin de faire 3 affectation mémoire
    printf("mncblas_cswap %d : nombre de cycles: %Ld \n", i, end - start);


  }


//-----------------------------------------------------------------------------------------------------------
  vector_init_double_comp(vec14, 1.f, 5.f);
  vector_init_double_comp(vec24, 2.f, 2.f);


  for (i = 0; i < NB_FOIS; i++)
  {
    start = _rdtsc () ;


    mncblas_zswap(VECSIZE, vec14, 1, vec24, 1);
    //printf ("mncblas_sswap %d : nombre de cycles: %lld \n", i, end-start) ;
    end = _rdtsc () ;
    val_moy4 += calcul_byte_swap ("mncblas_zswap ", VECSIZE * sizeof(double) * 6, end-start) ; // le *3 car on a besoin de faire 3 affectation mémoire
    printf("mncblas_zswap %d : nombre de cycles: %Ld \n", i, end - start);


  }

//vec size * size fo float
printf("\n=========== Résultats swap simple précision ===========\n");
printf("En moyenne on a une performance de %5.3f GBytes/s\n", val_moy/NB_FOIS);

printf("=========== Résultats swap double précision ===========\n");
printf("En moyenne on a une performance de %5.3f GBytes/s\n", val_moy2/NB_FOIS);

printf("=========== Résultats swap nombre complexe simple précision ===========\n");
printf("En moyenne on a une performance de %5.3f GBytes/s\n", val_moy3/NB_FOIS);

printf("=========== Résultats swap nombre complexe double précision ===========\n");
printf("En moyenne on a une performance de %5.3f GBytes/s\n", val_moy4/NB_FOIS);


return 0;
}
