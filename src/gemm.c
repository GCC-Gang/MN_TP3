#include "../include/mnblas.h"
#include "../include/complexe2.h"
#include <stdlib.h>



/*
— float, simple précision, s
— double, double précision, d
— complexe simple précision, sc
— complexe double précision, dz
*/

/*
void cblas_sgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const float alpha, const float *a, const MKL_INT lda, const float *b, const MKL_INT ldb, const float beta, float *c, const MKL_INT ldc);

C := alpha*op(A)*op(B) + beta*C
where:
op(X) is one of op(X) = X, or op(X) = XT, or op(X) = XH,
alpha and beta are scalars,
A, B and C are matrices:
op(A) is an m-by-k matrix,
op(B) is a k-by-n matrix,
C is an m-by-n matrix.

*/

//What is CBLAS_transpose ? --> I add MN

/*====================== Algo général =================
 A * B
 A = m * k
 B = k * n
X = ligne * colones

                    | 3 | 2 | 1 |
                    | 3 | 2 | 1 |
                    | 3 | 2 | 1 |

    | 1 | 2 | 3 |
    | 1 | 2 | 3 |
    | 1 | 2 | 3 |

            ligne,colone = m,k * k,n

Opérartions: (0,0)  = (0,0) * (0,0) + (0,1) * (1,0) + (0,2) * (2,0)
             (0,1) = ()



Pour chaque ligne


layout : indique si trié par colonne MNCblasColMajor ou par ligne MNCblasRowMajor ( voir compte rendu )

 TransA :
    Specifies the form of op(A) used in the matrix multiplication:
        if transa=CblasNoTrans, then op(A) = A;
        if transa=CblasTrans, then op(A) = AT;
        if transa=CblasConjTrans, then op(A) = AH.
TransB
    Specifies the form of op(B) used in the matrix multiplication:
        if transb=CblasNoTrans, then op(B) = B;
        if transb=CblasTrans, then op(B) = BT;
        if transb=CblasConjTrans, then op(B) = BH.
m : nombre de lignes
n : nombre de colonnes

a :  matrice de taille lda*k
    si CblasColMajor : k=n
    si CblasRowMajor : k=m

lda :  taille du tableau en une dimension
    si CblasColMajor : lda vaut au minimum entre 1 et m
    si CblasRowMajor : lda vaut au minimum entre 1 et n


C := alpha*op(A)*op(B) + beta*C

 MH = mattrice adjointe
 Exemple:
    |4+i    5|   =>     |4-i 2-2i|
    |2+2i   i|   =>     |5     -i|

*/

 float* transposition_float(const float* X, const int ligne,const int colonne,const int ldx){
    float* table2 = malloc(sizeof(float)*colonne*ligne);
    //float table2[colonne*ligne];
    #pragma omp parallel for
        for(register int i = 0; i < ligne; ++i){
            for(register int j =0; j < colonne; ++j){
                table2[i*ldx + j] = X[j*ldx + i];
            }
        } //on a maintenant la transposé de X dans table2
    return table2;
 }

 double* transposition_double(const double* X, const int ligne,const int colonne,const int ldx){
    double* table2 = malloc(sizeof(float)*colonne*ligne);
    //double table2[colonne*ligne];
        #pragma omp parallel for
        for(register int i = 0; i < ligne; ++i){
            for(register int j =0; j < colonne; ++j){
                table2[i*ldx + j] = X[j*ldx + i];
            }
        } //on a maintenant la transposé de X dans table2
    return table2;
 }

 complexe_float_t* transposition_float_comp(const void* X, const int ligne,const int colonne,const int ldx){
    complexe_float_t* table2 = malloc(sizeof(float)*colonne*ligne);
    //complexe_float_t table2[colonne*ligne];
    complexe_float_t* new_X = (complexe_float_t*) X;
        #pragma omp parallel for
        for(register int i = 0; i < ligne; ++i){
                for(register int j =0; j < colonne; ++j){
                table2[i*ldx + j].real = new_X[j*ldx + i].real;
                table2[i*ldx + j].imaginary = new_X[j*ldx + i].imaginary;

            }
        } //on a maintenant la transposé de X dans table2
    return table2;
 }

 complexe_double_t* transposition_double_comp(const void* X, const int ligne,const int colonne,const int ldx){
    complexe_double_t* table2 = malloc(sizeof(float)*colonne*ligne);
    //complexe_double_t table2[colonne*ligne];
    complexe_float_t* new_X = (complexe_float_t*) X;
        #pragma omp parallel for
        for(register int i = 0; i < ligne; ++i){
            for(register int j =0; j < colonne; ++j){
                table2[i*ldx + j].real = new_X[j*ldx + i].real;
                table2[i*ldx + j].real = new_X[j*ldx + i].real;

            }
        } //on a maintenant la transposé de X dans table2
    return table2;
 }

 complexe_float_t* adjointement_float(const void* X,const int ligne,const int colonne,const int ldx){
    complexe_float_t* table2 = malloc(sizeof(float)*colonne*ligne);
    //complexe_float_t table2[ligne*colonne];
    complexe_float_t* new_X = (complexe_float_t*)  X;
        #pragma omp parallel for
        for(register int i = 0; i < ligne; ++i){
            for(register int j =0; j < colonne; ++j){
                table2[i*ldx+j].real = new_X[j*ldx + i].real;
                table2[i*ldx+j].real = -new_X[j*ldx + i].imaginary;

                //table2[i*ldx + j] = X[j*ldx + i];
            }
        } //on a maintenant la transposé conjugué de X dans table2
    return table2;
 }

 complexe_double_t* adjointement_double(const void* X,const int ligne,const int colonne,const int ldx){
    complexe_double_t* table2 = malloc(sizeof(float)*colonne*ligne);
    //complexe_double_t table2[ligne*colonne];
    complexe_double_t* new_X = (complexe_double_t*) X;
        #pragma omp parallel for
        for(register int i = 0; i < ligne; ++i){
            for(register int j =0; j < colonne; ++j){
                table2[i*ldx+j].real = new_X[j*ldx + i].real;
                table2[i*ldx+j].imaginary = -new_X[j*ldx + i].imaginary;

                //table2[i*ldx + j] = X[j*ldx + i];
            }
        } //on a maintenant la transposé conjugué de X dans table2
    return table2;
 }


float* cpymat_f(const float* X, const int ligne, const int colonne, const int ldx){
    float* table2 = malloc(sizeof(float)*colonne*ligne);
    #pragma omp parallel for
    for(register int i = 0; i < ligne; ++i){
        for(register int j =0; j < colonne; ++j){

            table2[i*ldx+j] = X[i*ldx+j];
        }
    }
    return table2;
}

double* cpymat_d(const double* X, const int ligne,const int colonne, const int ldx){
    double* table2 = malloc(sizeof(double)*colonne*ligne);
    #pragma omp parallel for
    for(register int i = 0; i < ligne; ++i){
        for(register int j =0; j < colonne; ++j){
            table2[i*ldx+j] = X[i*ldx+j];
        }
    }
    return table2;
}

complexe_float_t* cpymat_c(const void* X, const int ligne,const int colonne, const int ldx){
    complexe_float_t* table2 = malloc(sizeof(complexe_float_t)*colonne*ligne);
    const complexe_float_t* new_X = (complexe_float_t*) X;
    #pragma omp parallel for
    for(register int i = 0; i < ligne; ++i){
        for(register int j =0; j < colonne; ++j){
            table2[i*ldx+j].real = new_X[i*ldx+j].real;
            table2[i*ldx+j].imaginary = new_X[i*ldx+j].imaginary;
        }
    }
    return table2;
}

complexe_double_t* cpymat_z(const void* X, const int ligne,const int colonne, const int ldx){
    complexe_double_t* table2 = malloc(sizeof(complexe_double_t)*colonne*ligne);
    const complexe_double_t* new_X = (complexe_double_t*) X;
    #pragma omp parallel for

    for(register int i = 0; i < ligne; ++i){
        for(register int j =0; j < colonne; ++j){
            table2[i*ldx+j].real = new_X[i*ldx+j].real;
            table2[i*ldx+j].imaginary = new_X[i*ldx+j].imaginary;
        }
    }
    return table2;
}


/*
if(TransA == MNCblasConjTrans){
            new_A = adjointement(A, M, K, lda);
        }else{
            new_A = A;
        }
        if(TransB == MNCblasConjTrans){
            new_A = adjointement(B, K, N, ldb);
        }else{
            new_B = B;
        }
*/

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc){

        register unsigned int f = 0;
        register unsigned int i = 0;
        register unsigned int j = 0;
        register unsigned int k = 0;
        register float somme;
        register unsigned int lg_C = M*N; // TODO : value is already in ldc
        float table_tmp[lg_C];


        register float *new_A;
        register float *new_B;

        if (TransA == MNCblasNoTrans) {
            new_A = cpymat_f(A, M, K, lda);
        } else {
            new_A = transposition_float(A, M, K, lda);
        }

        if (TransB == MNCblasNoTrans) {
            new_B = cpymat_f(B, K, N, ldb);
        } else {
            new_B = transposition_float(B, K, N, ldb);
        }
        //une fois les operations sur les matrices faites, on multiplie

       if (layout == MNCblasRowMajor){
        //cas ou trié par ligne

            //alpha*A*B + beta*C
            //on multiplie A par son scalaire

            //on multiplie A par B et on ajoute beta*C
                #pragma omp parallel for private(somme,j, k)

                for(i = 0; i < M; ++i){
                    for(j = 0; j < N; ++j){
                        //C[i] *=beta;
                        somme = 0;
           
                        for(k = 0; k < K; k++){
                            //lda represente taille ligne
                            //on multiplie A et B
                            somme += new_A[i*lda + k] * new_B[k * ldb + j];
                        }
                        table_tmp[i*ldc + j] = somme * alpha;
                    }
                }
            }
        if(layout == MNCblasColMajor){
            #pragma omp parallel for private(somme,j, k)

            for(i = 0; i < M; ++i){
                    for(j = 0; j < N; ++j){
                        //C[i] *=beta;
                        somme = 0;
                     
                        for(k = 0; k < K; k++){
                            //lda represente taille ligne
                            //on multiplie A et B
                            somme += new_A[k*lda + i] * new_B[j * ldb + k];
                        }
                        table_tmp[j*ldc + i] = somme * alpha;
                    }
                }
        }

        //il ne reste plus qu'a ajouter C correctement correpsond au + beta*C
        #pragma omp parallel for
        for(f = 0; f < lg_C; ++f){
            C[f] *= beta;
            C[f] += table_tmp[f];
        }

        free(new_A);
        free(new_B);
    }


void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){

        register unsigned int f = 0;
        register unsigned int i = 0;
        register unsigned int j = 0;
        register unsigned int k = 0;
        register double somme;
        register unsigned int lg_C = M*N;
        double table_tmp[lg_C];


        register double *new_A;
        register double *new_B;

        if(TransA == MNCblasNoTrans){
            new_A = cpymat_d(A, M, K, lda);
        }else{
            new_A = transposition_double(A, M, K, lda);
        }

        if(TransB == MNCblasNoTrans){
            new_B = cpymat_d(B, K, N, ldb);
        }else{
            new_B = transposition_double(B, K, N, ldb);
        }

        if (layout == MNCblasRowMajor){
        //cas ou trié par ligne

            //alpha*A*B + beta*C
            //on multiplie A par son scalaire

            //on multiplie A par B et on ajoute beta*C
                #pragma omp parallel for private(somme,j, k)
                for(i = 0; i < M; ++i){
                    for(j = 0; j < N; ++j){
                        //C[i] *=beta;
                        somme = 0;
                        for(k = 0; k < K; k++){
                            //lda represente taille ligne
                            //on multiplie A et B
                            somme += new_A[i*lda + k] * new_B[k * ldb + j];
                        }
                        table_tmp[i*ldc + j] = somme * alpha;
                    }
                }
            }
        if(layout == MNCblasColMajor){
            #pragma omp parallel for private(somme,j, k)

            for(i = 0; i < M; ++i){
                    for(j = 0; j < N; ++j){
                        //C[i] *=beta;
                        somme = 0;
                        for(k = 0; k < K; k++){
                            //lda represente taille ligne
                            //on multiplie A et B
                            somme += new_A[k*lda + i] * new_B[j * ldb + k];
                        }
                        table_tmp[j*ldc + i] = somme * alpha;
                    }
                }
        }

        //il ne reste plus qu'a ajouter C correctement correpsond au + beta*C
        #pragma omp parallel for
        for(f = 0; f < lg_C; ++f){
            C[f] *= beta;
            C[f] += table_tmp[f];
        }
        free(new_A);
        free(new_B);

}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){

        register unsigned int f = 0;
        register unsigned int i = 0;
        register unsigned int j = 0;
        register unsigned int k = 0;
        register complexe_float_t somme, tmp;
        register unsigned int lg_C = M*N;
        complexe_float_t table_tmp[lg_C];


        register complexe_float_t *new_A;
        register complexe_float_t *new_B;

        if(TransA == MNCblasTrans){
            new_A = transposition_float_comp(A, M, K, lda);
        }else if(TransA == MNCblasConjTrans){
            new_A = adjointement_float(A, M, K, lda);
        }else{
            new_A = cpymat_c(A, M, K, lda);
        }

        if(TransB == MNCblasTrans || TransB == MNCblasConjTrans){
            new_B = transposition_float_comp(B, K, N, ldb);
        }else if(TransB == MNCblasConjTrans){
            new_B = adjointement_float(B, K, N, ldb);
        }else{
            new_B = cpymat_c(B, K, N, ldb);
        }


        if (layout == MNCblasRowMajor){
        //cas ou trié par colone

            //alpha*A*B + beta*C
            //on multiplie A par son scalaire

            //on multiplie A par B et on ajoute beta*C
            #pragma omp parallel for private(tmp, somme, k)
            for(i = 0; i < M; ++i){
                for(j = 0; j < N; ++j){
                    //C[i] *=beta;
                    somme.real = 0;
                    somme.imaginary = 0;
                    for(k = 0; k < K; k++){
                        //lda represente taille ligne
                        //on multiplie A et B
                        tmp = mult_complexe_float(((complexe_float_t *) new_A)[i*lda + k], ((complexe_float_t *) new_B)[k * ldb + j]);
                        somme.real = tmp.real;
                        somme.imaginary = tmp.imaginary;
                    }

                    table_tmp[i*ldc + j] = mult_complexe_float(somme, *(complexe_float_t*) alpha);
                    }
                }
            }
        if(layout == MNCblasColMajor){
            #pragma omp parallel for private(tmp, somme, k)
            for(i = 0; i < M; ++i){
                for(j = 0; j < N; ++j){
                    //C[i] *=beta;
                    somme.real = 0;
                    somme.imaginary = 0;
                    for(k = 0; k < K; k++){
                        //lda represente taille ligne
                        //on multiplie A et B
                        tmp = mult_complexe_float(((complexe_float_t *) new_A)[k*lda + i], ((complexe_float_t *) new_B)[j * ldb + k]);
                        somme.real = tmp.real;
                        somme.imaginary = tmp.imaginary;
                    }

                    table_tmp[j*ldc + i] = mult_complexe_float(somme, *(complexe_float_t*) alpha);
                    }
                }

        }

        //il ne reste plus qu'a ajouter C correctement correpsond au + beta*C
        #pragma omp parallel for
        for(f = 0; f < lg_C; ++f){
            tmp.real = 0;
            tmp.imaginary = 0;
            tmp = mult_complexe_float(((complexe_float_t *) C)[f], *(complexe_float_t*) beta);
            ((complexe_float_t *) C)[f].real = table_tmp[f].real + tmp.real;
            ((complexe_float_t *) C)[f].real = table_tmp[f].real + tmp.real;

        }

        free(new_A);
        free(new_B);

}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){

        register unsigned int f = 0;
        register unsigned int i = 0;
        register unsigned int j = 0;
        register unsigned int k = 0;
        register complexe_float_t somme, tmp;
        register unsigned int lg_C = M*N;
        complexe_float_t table_tmp[lg_C];

        register complexe_double_t *new_A;
        register complexe_double_t *new_B;
        if(TransA == MNCblasTrans){
            new_A = transposition_double_comp(A, M, K, lda);
        }else if(TransA == MNCblasConjTrans){
            new_A = adjointement_double(A, M, K, lda);
        }else{
            new_A = cpymat_z(A, M, K, lda);
        }
        if(TransB == MNCblasTrans || TransB == MNCblasConjTrans){
            new_B = transposition_double_comp(B, K, N, ldb);
        }else if(TransB == MNCblasConjTrans){
            new_B = adjointement_double(B, K, N, ldb);
        }else{
            new_B = cpymat_z(B, K, N, ldb);
        }

        if (layout == MNCblasRowMajor){
        //cas ou trié par colone

            //alpha*A*B + beta*C
            //on multiplie A par son scalaire

            //on multiplie A par B et on ajoute beta*C
            #pragma omp parallel for private(tmp, somme, k)
            for(i = 0; i < M; ++i){
                for(j = 0; j < N; ++j){
                    //C[i] *=beta;
                    somme.real = 0;
                    somme.imaginary = 0;
                    for(k = 0; k < K; k++){
                        //lda represente taille ligne
                        //on multiplie A et B
                        tmp = mult_complexe_float(((complexe_float_t *) new_A)[i*lda + k], ((complexe_float_t *) new_B)[k * ldb + j]);
                        somme.real = tmp.real;
                        somme.imaginary = tmp.imaginary;
                    }

                    table_tmp[i*ldc + j] = mult_complexe_float(somme, *(complexe_float_t*) alpha);
                    }
                }
            }
        if(layout == MNCblasColMajor){
            #pragma omp parallel for private(tmp, somme, k)
            for(i = 0; i < M; ++i){
                for(j = 0; j < N; ++j){
                    //C[i] *=beta;
                    somme.real = 0;
                    somme.imaginary = 0;
                    for(k = 0; k < K; k++){
                        //lda represente taille ligne
                        //on multiplie A et B
                        tmp = mult_complexe_float(((complexe_float_t *) new_A)[k*lda + i], ((complexe_float_t *) new_B)[j * ldb + k]);
                        somme.real = tmp.real;
                        somme.imaginary = tmp.imaginary;
                    }

                    table_tmp[j*ldc + i] = mult_complexe_float(somme, *(complexe_float_t*) alpha);
                    }
                }

        }

        //il ne reste plus qu'a ajouter C correctement correpsond au + beta*C
        #pragma omp parallel for
        for(f = 0; f < lg_C; ++f){
            tmp.real = 0;
            tmp.imaginary = 0;
            tmp = mult_complexe_float(((complexe_float_t *) C)[f], *(complexe_float_t*) beta);
            ((complexe_float_t *) C)[f].real = table_tmp[f].real + tmp.real;
            ((complexe_float_t *) C)[f].real = table_tmp[f].real + tmp.real;

        }

        free(new_A);
        free(new_B);

}
