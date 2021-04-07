#include "../include/mnblas.h"
#include "../include/complexe2.h"


/*
* no FLOP in this context
*/
CBLAS_INDEX mnblas_isamax(const int N, const float* X, const int incX) {
    register int i = 0;
    register CBLAS_INDEX index = 0;
    register float max = 0.f;

    for (; i < N; i += incX) {
        register float tmp = X[i];
        if (tmp < 0) {
            tmp = -tmp;
        }

        if (tmp > max) {
            max = tmp;
            index = i;
        }
    }

    return index;
}

CBLAS_INDEX mnblas_idamax(const int N, const double *X, const int incX) {
    register int i = 0;
    register CBLAS_INDEX index = 0;
    register double max = 0.f;

    for (; i < N; i += incX) {
        register double tmp = X[i];
        if (tmp < 0) {
            tmp = -tmp;
        }

        if (tmp > max) {
            max = tmp;
            index = i;
        }
    }

    return index;
}

/*
* no FLOP in this context
*/
CBLAS_INDEX mnblas_icamax(const int N, const void* X, const int incX) {
    register int i = 0;
    register CBLAS_INDEX index = 0;
    register float max = 0.f;

    for (; i < N; i += incX) {
        register float tmp1 = ((complexe_float_t*) X)[i].real;
        if (tmp1 < 0) {
            tmp1 = -tmp1;
        }

        register float tmp2 = ((complexe_float_t*) X)[i].imaginary;
        if (tmp2 < 0) {
            tmp2 = -tmp2;
        }

        register float tmp3 = tmp1 + tmp2;

        if (tmp3 > max) {
            max = tmp3;
            index = i;
        }
    }

    return index;
}

CBLAS_INDEX mnblas_izamax(const int N, const void* X, const int incX) {
    register int i = 0;
    register CBLAS_INDEX index = 0;
    register double max = 0.f;

    for (; i < N; i += incX) {
        register double tmp1 = ((complexe_float_t*) X)[i].real;
        if (tmp1 < 0) {
            tmp1 = -tmp1;
        }

        register double tmp2 = ((complexe_float_t*) X)[i].imaginary;
        if (tmp2 < 0) {
            tmp2 = -tmp2;
        }

        register double tmp3 = tmp1 + tmp2;

        if (tmp3 > max) {
            max = tmp3;
            index = i;
        }
    }

    return index;
}
