#include <stdio.h>
#include <x86intrin.h>

#include "../include/mnblas.h"
#include "../include/complexe2.h"

#include "flop.h"

#define NB_FOIS 12000
#define VECSIZE 700

int main(int argc, char **argv)
{
    unsigned long long start, end;
    register int i;

    init_flop();
    float matrice[VECSIZE * VECSIZE];
    float vecteur1[VECSIZE];
    float *vecteur2 = malloc(sizeof(float) * VECSIZE);
    float moyenne_FLOP = 0;
    float moyenne_FLOP2 = 0;

    for (i = 0; i < NB_FOIS; i++)
    {
        for (int j = 0; j < VECSIZE; j++)
        {

            vecteur1[j] = 5.0;
            vecteur2[j] = 7.0;
        }

        for (int k = 0; k < VECSIZE * VECSIZE; k++)
        {
            matrice[k] = 3.0;
        }

        start = _rdtsc();
        mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, 6, matrice, VECSIZE * VECSIZE, vecteur1, 1, 4, vecteur2, 1);
        end = _rdtsc();

        //printf("%f\n", vecteur2[0]);

        printf("mncblas_sgemv %d : nombre de cycles: %Ld \n", i, end - start);
        moyenne_FLOP += calcul_flop_return("sdot ", 3 * VECSIZE * VECSIZE + 2 * VECSIZE, end - start);

        printf("\n");
    }
    moyenne_FLOP /= (float)NB_FOIS;
    complexe_float_t mat[VECSIZE * VECSIZE];
    complexe_float_t v1[VECSIZE];
    complexe_float_t v2[VECSIZE];
    complexe_float_t alp, bet;
    /*
    En rowMajor :
    5 + i2 5 + i2 5 + i2             4+ 2i         6+3i               234+372i        6+3i             258 +i399
    5 + i2 5 + i2 5 + i2  * 6 + i *  4+ 2i     +   6+3i * 5 +i2   =   234+372i   +    6+3i  * 5+i2 =   258 +i399
    5 + i2 5 + i2 5 + i2             4+ 2i         6+3i               234+372i        6+3i             258 +i399

    Attendu 450 + i135 en MNCblasConjTrans et M et N= 3
    */

    alp.real = 6.;
    alp.imaginary = 1.;

    bet.real = 5.;
    bet.imaginary = 2.;

    for (int j = 0; j < VECSIZE * VECSIZE; j++)
    {
        mat[j].real = 5.;
        mat[j].imaginary = 2.0;
    }

    for (int j = 0; j < VECSIZE; j++)
    {
        v1[j].real = 4.;
        v1[j].imaginary = 2.;
    }
    for (i = 0; i < NB_FOIS; i++)
    {
        for (int j = 0; j < VECSIZE; j++)
        {

            v2[j].real = 6.;
            v2[j].imaginary = 3.;
        }

        start = _rdtsc();
        mncblas_cgemv(MNCblasColMajor, MNCblasConjTrans, VECSIZE, VECSIZE, &alp, mat, VECSIZE * VECSIZE, v1, 1, &bet, v2, 1);
        end = _rdtsc();

        //printf("%f + i%f\n", v2[0].real, v2[0].imaginary);

        printf("mncblas_cgemv %d :  nombre de cycles: %Ld \n", i, end - start);
        moyenne_FLOP2 += calcul_flop_return("sdot ", 14 * VECSIZE * VECSIZE + VECSIZE * VECSIZE + 8 * VECSIZE, end - start);

        printf("\n");
    }
    moyenne_FLOP2 /= (float)NB_FOIS;
    free(vecteur2);
        printf("\n");

    printf("\n=========== Résultats BLAS2 simple précision ===========\n");
    printf("En moyenne on a une performance de %5.3f GFlop/s\n", moyenne_FLOP);
      printf("\n=========== Résultats BLAS2 simple précision complexe ===========\n");
    printf("En moyenne on a une performance de %5.3f GFlop/s\n", moyenne_FLOP2);

    /*  const float matrice[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    const float vecteur1[3] = {1.0, 2.0, 3.0};
    float vecteur2[3] = {7.0, 8.0, 9.0};*/
    /*
    En rowMajor :
    1 2 3         1         7           84        7         112
    4 5 6  * 6 *  2     +   8 * 4   =   192   +   8 * 4 =   224
    7 8 9         3         9           300       9         336
    */
    /*mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, 3, 3, 6, matrice, 9, vecteur1, 1, 4, vecteur2, 1);
    for (int i = 0; i < 3; i++)
    {
        printf("%f\n", vecteur2[i]);
    }
    printf("\n");*/
    /*
    En colMajor :
    1 4 7         1         7           180       7         208
    2 5 8  * 6 *  2     +   8 * 4   =   216   +   8 * 4 =   248
    3 6 9         3         9           252       9         288
    */

    /*  vecteur2[0] = 7.0;
    vecteur2[1] = 8.0;
    vecteur2[2] = 9.0;

    mncblas_sgemv(MNCblasColMajor, MNCblasNoTrans, 3, 3, 6, matrice, 9, vecteur1, 1, 4, vecteur2, 1);
    for (int i = 0; i < 3; i++)
    {
        printf("%f\n", vecteur2[i]);
    }

    complexe_float_t mat[9];
    complexe_float_t v1[3];
    complexe_float_t v2[3];
    complexe_float_t alp, bet;
    alp.real = 6;
    alp.imaginary = 1;

    bet.real = 5;
    bet.imaginary = 2;
    for (int i = 0; i < 9; i++)
    {
        mat[i].real = 5;
        mat[i].imaginary = 2;
    }
    for (int i = 0; i < 3; i++)
    {
        v1[i].real = 4;
        v1[i].imaginary = 2;

        v2[i].real = 6;
        v2[i].imaginary = 3;
    }*/

    /*
    En rowMajor :
    5 + i2 5 + i2 5 + i2             4+ 2i         6+3i               234+372i        6+3i             258 +i399
    5 + i2 5 + i2 5 + i2  * 6 + i *  4+ 2i     +   6+3i * 5 +i2   =   234+372i   +    6+3i  * 5+i2 =   258 +i399
    5 + i2 5 + i2 5 + i2             4+ 2i         6+3i               234+372i        6+3i             258 +i399

    Attendu 450 + i135 en MNCblasConjTrans
    */

    /* mncblas_cgemv(MNCblasColMajor, MNCblasConjTrans, 3, 3, &alp, mat, 9, v1, 1, &bet, v2, 1);
    for (int i = 0; i < 3; i++)
    {
        printf("%f + i%f\n", v2[i].real, v2[i].imaginary);
    }
    return 0;*/
}
