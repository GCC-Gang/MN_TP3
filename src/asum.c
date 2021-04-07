#include "../include/mnblas.h"
#include "../include/complexe2.h"

/*
* 1 FLOP
*/
float mnblas_sasum(const int N, const float* X, const int incX) {
    register int i = 0;
    register float res = 0.f;

    for (; i < N; i += incX) {
        register float tmp = X[i];
        res += (tmp >= 0 ? tmp : -tmp);
    }

    return res;
}

double mnblas_dasum(const int N, const double* X, const int incX) {
    register int i = 0;
    register double res = 0.f;

    for (; i < N; i += incX) {
        register double tmp = X[i];
        res += (tmp >= 0 ? tmp : -tmp);
    }

    return res;
}

/*
* 2 FLOP
*/
float mnblas_scasum(const int N, const void* X, const int incX) {
    register int i = 0;
    register float res = 0.f;

    for (; i < N; i += incX) {
        register float tmp = ((complexe_float_t*) X)[i].real;
        res += (tmp >= 0 ? tmp : -tmp);

        tmp = ((complexe_float_t*) X)[i].imaginary;
        res += (tmp >= 0 ? tmp : -tmp);
    }

    return res;
}

double mnblas_dzasum(const int N, const void* X, const int incX) {
    register int i = 0;
    register double res = 0.f;

    for (; i < N; i += incX) {
        register double tmp = ((complexe_double_t*) X)[i].real;
        res += (tmp >= 0 ? tmp : -tmp);

        tmp = ((complexe_double_t*) X)[i].imaginary;
        res += (tmp >= 0 ? tmp : -tmp);
    }

    return res;
}
