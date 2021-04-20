#include <stdio.h>
#include <x86intrin.h>

#include "../include/mnblas.h"
#include "../include/complexe.h"

#include "flop.h"

#define NB_FOIS 250
#define NB_FOIS2 250


#define MATSIZE 250000 // pour des matrices de taille 500 par 500

typedef double mdouble[MATSIZE];
typedef complexe_double_t mComplexedouble[MATSIZE];

mdouble mat1, mat2, mat3;
mComplexedouble matc1, matc2, matc3;

void MAT_init(mdouble m, double x)
{
  register unsigned int i;

  for (i = 0; i < MATSIZE; i++){
    m[i] = x;
  }

  return;
}

void MAT_init_comp(mComplexedouble m, double x, double y)
{
  register unsigned int i;

  for (i = 0; i < MATSIZE; i++){
    m[i].real = x;
    m[i].imaginary = y;
  }


  return;
}

int main(int argc, char **argv)
{
    // unsigned long long start, end;
    // register int i;
/*
    int len = 6;

    const float matriceA[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // MxK = 2x3
    const float matrice1[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // KxN = 3x2
    const float matriceB[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; // KxN = 3x2
    float* matriceC = calloc(4, sizeof(float));         // MxN = 2x2
    mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans,
                  2, 2, 3, 1.0, matriceA, 3 , matriceB, 2,
                  1.0, matriceC, 2);
        */
       /*
    for(int i = 0; i < len; i++){
        printf("matrice %i: %f\n",i,matriceA[i]);
    }
    //transposition
    mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans,
                  2, 2, 3, 1.0, matriceA, 3, matrice1, 2,
                  1.0, matriceC, 2);

     for(int i = 0; i < len; i++){
        printf("matrice %i: %f\n",i,matriceC[i]);
    }

    const float* matriceAA = calloc(100, sizeof(float)); // MxK = 2x3
    const float* matriceBB = calloc(100, sizeof(float)); // MxK = 2x3
    float* matriceC2 = calloc(100, sizeof(float)); // KxN = 3x2

    for(int i = 0; i < 100; i++){
        printf("matrice %i: %f\n",i,matriceC2[i]);
    }


    mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans,
                  10, 10, 10, 1.0, matriceAA, 10, matriceBB, 10,
                  1.0, matriceC2, 10);

    for(int i = 0; i < 100; i++){
        printf("matrice %i: %f\n",i,matriceC2[i]);
    }
    */
    //============================================= Ecriture des vrais tests ================================
    double somme_double = 0.f;

    unsigned long long start, end;
    //float res[MATSIZE];
    int i;

    init_flop();

    complexe_double_t alpha = {1.0, 2.0};
    void* new_alpha = (void* ) &alpha;
    complexe_double_t beta = {2.0, 1.0};
    void* new_beta = (void* ) &beta;

    unsigned int nb_ope = 125000000 * 2 + 250000;
    for (i = 0; i < NB_FOIS; i++)
    {
    //on initialise nos 3 vecteurs, sinon risque de dépassement dans C vu que le resultat y est stocké
    MAT_init(mat1, 1.0);
    MAT_init(mat2, 1.0);
    MAT_init(mat3, 15.0);

    start = _rdtsc();

    mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans,
                  500, 500, 500, 1.0, mat1, 500, mat2, 500,
                  1.0, mat3, 2);

    end = _rdtsc();

    printf("mncblas_dgemm %d : nombre de cycles: %Ld \n", i, end - start);
    somme_double += calcul_flop_return("mncblas_dgemm", nb_ope , end - start); //nombre d'opération = m*n*p = 500^3 = 125000000
    }
    printf("\n");


    double somme_double_comp = 0.f;
    unsigned int nb_ope2 = 40*125000000+250000;
    for (i = 0; i < NB_FOIS2; i++)
    {
    //on initialise nos 3 vecteurs, sinon risque de dépassement dans C vu que le resultat y est stocké
    MAT_init_comp(matc1, 1.0, 5.0);
    MAT_init_comp(matc2, 5.0, 1.0);
    MAT_init_comp(matc3, 15.0, 3.0);



    start = _rdtsc();

    mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans,
                  500, 500, 500, new_alpha, matc1, 500, matc2, 500,
                  new_beta, matc3, 2);

    end = _rdtsc();

    printf("mncblas_zgemm %d : nombre de cycles: %Ld \n", i, end - start);
    somme_double_comp += calcul_flop_return("mncblas_zgemm", nb_ope2 , end - start); //nombre d'opération = m*n*p = 500^3 = 125000000
    }
    printf("\n");

    printf("\n=========== Résultats BLAS3 double précision ===========\n");
    printf("En moyenne on a une performance de %5.3f GFlop/s\n", somme_double/NB_FOIS);
      printf("\n=========== Résultats BLAS3 double précision complexe ===========\n");
    printf("En moyenne on a une performance de %5.3f GFlop/s\n", somme_double_comp/NB_FOIS2);
    return 0;
}
