#include "../include/mnblas.h"
#include "../include/complexe.h"


/*
— float, simple précision, s
— double, double précision, d
— complexe simple précision, sc
— complexe double précision, dz
*/


void mncblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY)
{ //swap simple précision float
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY)
{ //Swap double précision double
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;
  // c'est le point virgule qui définit les différentes instructions ainsi:
  //i et j sont déjà déclarés
  // tant que i < N ET j < N
  // on incremente i de incX ou j de incY
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY){
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
  }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX, void *Y, const int incY){
//Swap complexe simple précision float
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  //complexe_float_t* save = 0x0;
  complexe_float_t save;

  //il faut trasnformer en complexe;
  //les parenthèses servent à indiquer que l'on transtype X et non la valeur du tableau de X
  complexe_float_t* new_X = X;
  complexe_float_t* new_Y = Y;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY){
    save.real = new_Y[j].real;
    save.imaginary = new_Y[j].imaginary;

    new_Y[j].real =  new_X[i].real ;
    new_Y[j].imaginary =  new_X[i].imaginary ;

    new_X[i].real = save.real ;
    new_X[i].imaginary = save.imaginary ;
  }


  return ;

}

void mncblas_zswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{// Swap complexe double precision double
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t save;
  //il faut trasnformer en complexe;
  complexe_float_t* new_X = X;
  complexe_float_t* new_Y = Y;
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY){
    save.real = new_Y[j].real;
    save.imaginary = new_Y[j].imaginary;

    new_Y[j].real =  new_X[i].real ;
    new_Y[j].imaginary =  new_X[i].imaginary ;

    new_X[i].real = save.real ;
    new_X[i].imaginary = save.imaginary ;
  }

  return ;
}
