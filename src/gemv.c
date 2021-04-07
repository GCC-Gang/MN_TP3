#include "../include/mnblas.h"
#include "../include/complexe2.h"
#include <stdlib.h>
#include <stdio.h>

/*
 layout : indique si trié par colonne MNCblasColMajor ou par ligne MNCblasRowMajor ( voir compte rendu )

 TransA :
    Specifies the operation:
    if trans = CblasNoTrans, then y:= alpha*A*x+ beta*y;
    if trans = CblasTrans, then y:= alpha*A'*x+ beta*y ;
    if trans = CblasConjTrans , then y:= alpha*conjg(A')*x + beta*y.
f
m : nombre de lignes
n : nombre de colonnes

alpha : scalaire
beta : scalaire

a :  matrice de taille lda*k
    si CblasColMajor : k=n
    si CblasRowMajor : k=m

lda :  taille du tableau en une dimension
    si CblasColMajor : lda vaut au minimum entre 1 et m
    si CblasRowMajor : lda vaut au minimum entre 1 et n

x : vecteur
incX :incrément de X

y : vecteur
incY : incrément de Y
*/

/*  3 x N/incX x M/incY + 2 x M/incY FLOP */
void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha,
                   const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY)
{

    register int indice_matrice = 0;

    register int colonne = 0;
    register int ligne = 0;
    /*
    Exemple de matrice (avec  m = 2 et n = 3)
    1 2 3
    4 5 6
    et sa transposée :
    1 4
    2 5
    3 6

    Stockée dans le tableau ainsi en RowMajor ainsi :
     A = [ 1, 2, 3, 4, 5, 6 ]
     et en colMajor ainsi :
     A =  [ 1, 4, 2, 5, 3, 6 ]
     sa transposée A' = [1, 2, 3, 4, 5, 6]
     En lisant A dans l'ordre naturel des indices c'est donc comme si on lisait sa transposée ( en colMajor).


    */
    float y_temp[M];
    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {
        for (ligne = 0; ligne < M; ligne += incY)
        {
            y_temp[ligne] = 0.;
            for (colonne = 0; colonne < N; colonne += incX)
            {
                y_temp[ligne] = y_temp[ligne] + A[indice_matrice] * alpha * X[colonne]; // 3 x N x M FLOPS
                indice_matrice++;
            }
            Y[ligne] = y_temp[ligne] + beta * Y[ligne]; // 2 x M FLOP
        }
    }

    else if ((layout == MNCblasColMajor && TransA == MNCblasNoTrans) || (layout == MNCblasRowMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {

        for (ligne = 0; ligne < M; ligne += incY)
        {
            y_temp[ligne] = 0.;
            indice_matrice = ligne;
            for (colonne = 0; colonne < N; colonne += incX)
            {
                y_temp[ligne] = y_temp[ligne] + A[indice_matrice] * alpha * X[colonne];

                indice_matrice = indice_matrice + M;
            }
            Y[ligne] = y_temp[ligne] + beta * Y[ligne];
        }
    }
}

/* 3 x N/incX x M/incY + 2 x M/incY FLOP */
void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda,
                   const double *X, const int incX, const double beta, double *Y, const int incY)
{
    register int indice_matrice = 0;
    register int colonne = 0;
    register int ligne = 0;
    double y_temp[M];
    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {
        for (ligne = 0; ligne < M; ligne += incY)
        {
            y_temp[ligne] = 0.;
            for (colonne = 0; colonne < N; colonne += incX)
            {
                y_temp[ligne] = y_temp[ligne] + A[indice_matrice] * alpha * X[colonne];
                indice_matrice++;
            }
            Y[ligne] = y_temp[ligne] + beta * Y[ligne];
        }
    }
    else if ((layout == MNCblasColMajor && TransA == MNCblasNoTrans) || (layout == MNCblasRowMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {

        for (ligne = 0; ligne < M; ligne += incY)
        {
            y_temp[ligne] = 0.;
            indice_matrice = ligne;
            for (colonne = 0; colonne < N; colonne += incX)
            {
                y_temp[ligne] = y_temp[ligne] + A[indice_matrice] * alpha * X[colonne];

                indice_matrice = indice_matrice + M;
            }
            Y[ligne] = y_temp[ligne] + beta * Y[ligne];
        }
    }
}

/* 14 * N/incX x M/incY + (N/incX x M/incY)&MNCblasConjTrans + 8 x M/incY  FLOP */
void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta, void *Y, const int incY)
{
    register int indice_matrice = 0;
    register int colonne = 0;
    register int ligne = 0;
    complexe_float_t y_temp[M];
    complexe_float_t y_temp_colonne[M];
    complexe_float_t *A_complexe = malloc(sizeof(complexe_float_t) * lda);

    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {
        for (ligne = 0; ligne < M; ligne += incY)
        {
            y_temp[ligne].real = 0.;
            y_temp[ligne].imaginary = 0.;
            y_temp_colonne[ligne].real = 0.;
            y_temp_colonne[ligne].imaginary = 0.;

            for (colonne = 0; colonne < N; colonne += incX)
            {
                A_complexe[indice_matrice].real = ((complexe_float_t *)A)[indice_matrice].real;
                if (TransA == MNCblasConjTrans)
                {
                    A_complexe[indice_matrice].imaginary = -((complexe_float_t *)A)[indice_matrice].imaginary; // N/incX x M/incY
                }
                else
                {
                    A_complexe[indice_matrice].imaginary = ((complexe_float_t *)A)[indice_matrice].imaginary;
                }

                y_temp_colonne[ligne] = mult_complexe_float(A_complexe[indice_matrice], *(complexe_float_t *)alpha); // 6 * N/incX x M/incY
                // printf("- %f + i%f\n", ((complexe_float_t *)alpha)->real, ((complexe_float_t *)alpha)->imaginary);
                y_temp_colonne[ligne] = mult_complexe_float(y_temp_colonne[ligne], ((complexe_float_t *)X)[colonne]); // 6 * N/incX x M/incY
                y_temp[ligne] = add_complexe_float(y_temp[ligne], y_temp_colonne[ligne]);                             //2 * N/incX x M/incY
                indice_matrice++;
            }
            ((complexe_float_t *)Y)[ligne] = mult_complexe_float(((complexe_float_t *)Y)[ligne], *(complexe_float_t *)beta); // 6 x M/incY
            ((complexe_float_t *)Y)[ligne] = add_complexe_float(y_temp[ligne], ((complexe_float_t *)Y)[ligne]);              // 2 x M/incY
        }
    }
    else if ((layout == MNCblasColMajor && TransA == MNCblasNoTrans) || (layout == MNCblasRowMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {

        for (ligne = 0; ligne < M; ligne += incY)
        {
            indice_matrice = ligne;
            y_temp[ligne].real = 0.;
            y_temp[ligne].imaginary = 0.;
            y_temp_colonne[ligne].real = 0.;
            y_temp_colonne[ligne].imaginary = 0.;
            for (colonne = 0; colonne < N; colonne += incX)
            {
                A_complexe[indice_matrice].real = ((complexe_float_t *)A)[indice_matrice].real;
                if (TransA == MNCblasConjTrans)
                {
                    A_complexe[indice_matrice].imaginary = -((complexe_float_t *)A)[indice_matrice].imaginary; // N/incX x M/incY
                }
                else
                {
                    A_complexe[indice_matrice].imaginary = ((complexe_float_t *)A)[indice_matrice].imaginary;
                }
                y_temp_colonne[ligne] = mult_complexe_float(((complexe_float_t *)A)[indice_matrice], *(complexe_float_t *)alpha);
                y_temp_colonne[ligne] = mult_complexe_float(y_temp_colonne[ligne], ((complexe_float_t *)X)[colonne]);
                y_temp[ligne] = add_complexe_float(y_temp[ligne], y_temp_colonne[ligne]);

                indice_matrice = indice_matrice + M;
            }
            ((complexe_float_t *)Y)[ligne] = mult_complexe_float(((complexe_float_t *)Y)[ligne], *(complexe_float_t *)beta);
            ((complexe_float_t *)Y)[ligne] = add_complexe_float(y_temp[ligne], ((complexe_float_t *)Y)[ligne]);
        }
    }

    free(A_complexe);
}

/* 14 * N/incX x M/incY + (N/incX x M/incY)&MNCblasConjTrans + 8 x M/incY  FLOP */
void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta, void *Y, const int incY)
{
    register int indice_matrice = 0;
    register int colonne = 0;
    register int ligne = 0;
    complexe_double_t y_temp[M];
    complexe_double_t y_temp_colonne[M];
    complexe_double_t *A_complexe = malloc(sizeof(complexe_double_t) * lda);
    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {
        for (ligne = 0; ligne < M; ligne += incY)
        {
            y_temp[ligne].real = 0.;
            y_temp[ligne].imaginary = 0.;
            y_temp_colonne[ligne].real = 0.;
            y_temp_colonne[ligne].imaginary = 0.;
            for (colonne = 0; colonne < N; colonne += incX)
            {
                A_complexe[indice_matrice].real = ((complexe_double_t *)A)[indice_matrice].real;
                if (TransA == MNCblasConjTrans)
                {
                    A_complexe[indice_matrice].imaginary = -((complexe_double_t *)A)[indice_matrice].imaginary; // N/incX x M/incY
                }
                else
                {
                    A_complexe[indice_matrice].imaginary = ((complexe_double_t *)A)[indice_matrice].imaginary;
                }
                y_temp_colonne[ligne] = mult_complexe_double(A_complexe[indice_matrice], *(complexe_double_t *)alpha);
                y_temp_colonne[ligne] = mult_complexe_double(y_temp_colonne[ligne], ((complexe_double_t *)X)[colonne]);
                y_temp[ligne] = add_complexe_double(y_temp[ligne], y_temp_colonne[ligne]);

                indice_matrice++;
            }
            ((complexe_double_t *)Y)[ligne] = mult_complexe_double(((complexe_double_t *)Y)[ligne], *(complexe_double_t *)beta);
            ((complexe_double_t *)Y)[ligne] = add_complexe_double(y_temp[ligne], ((complexe_double_t *)Y)[ligne]);
        }
    }
    else if ((layout == MNCblasColMajor && TransA == MNCblasNoTrans) || (layout == MNCblasRowMajor && (TransA == MNCblasTrans || TransA == MNCblasConjTrans)))
    {

        for (ligne = 0; ligne < M; ligne += incY)
        {
            y_temp[ligne].real = 0.;
            y_temp[ligne].imaginary = 0.;
            y_temp_colonne[ligne].real = 0.;
            y_temp_colonne[ligne].imaginary = 0.;
            indice_matrice = ligne;
            for (colonne = 0; colonne < N; colonne += incX)
            {
                A_complexe[indice_matrice].real = ((complexe_double_t *)A)[indice_matrice].real;
                if (TransA == MNCblasConjTrans)
                {
                    A_complexe[indice_matrice].imaginary = -((complexe_double_t *)A)[indice_matrice].imaginary; // N/incX x M/incY
                }
                else
                {
                    A_complexe[indice_matrice].imaginary = ((complexe_double_t *)A)[indice_matrice].imaginary;
                }
                y_temp_colonne[ligne] = mult_complexe_double(((complexe_double_t *)A)[indice_matrice], *(complexe_double_t *)alpha);
                y_temp_colonne[ligne] = mult_complexe_double(y_temp_colonne[ligne], ((complexe_double_t *)X)[colonne]);
                y_temp[ligne] = add_complexe_double(y_temp[ligne], y_temp_colonne[ligne]);

                indice_matrice = indice_matrice + M;
            }
            ((complexe_double_t *)Y)[ligne] = mult_complexe_double(((complexe_double_t *)Y)[ligne], *(complexe_double_t *)beta);
            ((complexe_double_t *)Y)[ligne] = add_complexe_double(y_temp[ligne], ((complexe_double_t *)Y)[ligne]);
        }
    }

    free(A_complexe);
}
