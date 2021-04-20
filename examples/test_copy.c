#include <stdio.h>
#include <x86intrin.h>

#include "../include/mnblas.h"
#include "../include/complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    100000

typedef __m128 vfloat[VECSIZE];
typedef complexe_float_t vcomplexe[VECSIZE];

vfloat vec1, vec2;
vcomplexe vec3, vec4;

void vector_init(vfloat V, float x) {
    unsigned int i;

    float tab[4] __attribute__ ((aligned(16))) = {x, x, x, x};

    for (i = 0; i < VECSIZE; i++)
        V[i] = _mm_load_ps(tab);

    return;
}

void vectorC_init(vcomplexe V, float a, float b) {
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++) {
        V[i].real = a;
        V[i].imaginary = b;
    }

    return;
}

//void vector_print(vfloat V) {
//    register unsigned int i;
//
//    for (i = 0; i < VECSIZE; i++)
//        printf("%f ", V[i]);
//    printf("\n");
//
//    return;
//}

void vectorC_print(vcomplexe V) {
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("%f %f, ", V[i].real, V[i].imaginary);
    printf("\n");

    return;
}

int main(int argc, char **argv) {
    unsigned long long start, end;
    int i;
    float time_sum = 0.f;

    init_flop();

    for (i = 0; i < NB_FOIS; i++) {
        vector_init(vec1, 1.0);

        start = _rdtsc();
        mncblas_scopy(VECSIZE, vec1, 1, vec2, 1);
        end = _rdtsc();

        // printf("mncblas_scopy %d : nombre de cycles: %Ld \n", i, end - start);
        time_sum += calcul_byte_swap("scopy ", 4 * VECSIZE * sizeof(float), end - start);
    }

    printf("Average copy time = %f\n", time_sum / NB_FOIS);

//   printf("=========================================\n");
//
//   for (i = 0; i < NB_FOIS; i++) {
//       vectorC_init(vec3, 1.0, 2.0);
//
//       start = _rdtsc () ;
//        mncblas_ccopy (VECSIZE, vec3, 1, vec4, 1) ;
//       end = _rdtsc () ;
//
//       printf ("mncblas_ccopy %d : nombre de cycles: %Ld \n", i, end-start) ;
//       calcul_byte ("ccopy ", VECSIZE * sizeof(float) * 2, end-start) ;
//   }

}
