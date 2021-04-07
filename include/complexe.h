#ifndef COMPLEXE_H
#define COMPLEXE_H



//#include <xmmintrin.h>

/*
#include <xmmintrin.h> //– SSE:
#include <emmintrin.h> //– SSE:
#include <pmmintrin.h> //– SSE3: 
#include <tmmintrin.h> //– SSSE3:
#include <smmintrin.h> //– SSE4.1:
#include <nmmintrin.h> // – SSSE4.2:
*/
//dot, axpy, copy
// __m128
// __m128d 

typedef struct {
  float real ;
  float imaginary ;
} complexe_float_t ;

typedef struct {
  double real ;
  double imaginary ;
} complexe_double_t ;

complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2) ;

complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2) ;

complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2) ;

complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2) ;

complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2) ;

complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2) ;


#endif
