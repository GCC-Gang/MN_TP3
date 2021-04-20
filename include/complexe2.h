#include <stdio.h>

#include <xmmintrin.h>

typedef struct {
    __m128 real;
    __m128 imaginary;
} complexe_float_t;

typedef struct {
    double real;
    double imaginary;
} complexe_double_t;

inline complexe_float_t add_complexe_float(const complexe_float_t c1, const complexe_float_t c2) {
    complexe_float_t r;

    r.real = _mm_add_ps(c1.real, c2.real);
    r.imaginary = _mm_add_ps(c1.imaginary, c2.imaginary);

    return r;
}

inline complexe_double_t add_complexe_double(const complexe_double_t c1, const complexe_double_t c2) {
    complexe_double_t r;

    r.real = c1.real + c2.real;
    r.imaginary = c1.imaginary + c2.imaginary;

    return r;
}

inline complexe_float_t mult_complexe_float(const complexe_float_t c1, const complexe_float_t c2) {
    complexe_float_t r;

//  r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
//  r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;

    r.real = _mm_sub_ps(_mm_mul_ps(c1.real, c2.real), _mm_mul_ps(c1.imaginary, c2.imaginary));
    r.imaginary = _mm_add_ps(_mm_mul_ps(c1.real, c2.imaginary), _mm_mul_ps(c1.imaginary, c2.real));

    return r;
}

inline complexe_double_t mult_complexe_double(const complexe_double_t c1, const complexe_double_t c2) {
    complexe_double_t r;

    r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
    r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;

    return r;
}

inline complexe_float_t div_complexe_float(const complexe_float_t c1, const complexe_float_t c2) {
//    complexe_float_t div;
//    complexe_float_t c2_conjuge = c2;
//
//    c2_conjuge.imaginary = -c2.imaginary;
//
//    float denominateur = c2.real * c2_conjuge.real - c2.imaginary * c2_conjuge.imaginary;
//
//    complexe_float_t nominateur = mult_complexe_float(c1, c2_conjuge);
//
//    div.real = nominateur.real / denominateur;
//    div.imaginary = nominateur.imaginary / denominateur;

    complexe_float_t div;
    complexe_float_t c2_conjuge = c2;

    float tab[4] __attribute__ ((aligned(16))) = {-1, -1, -1, -1};
    __m128 temp = _mm_load_ps(tab);
    c2_conjuge.imaginary = -_mm_mul_ps(c2.imaginary, temp);

    __m128 denominateur = _mm_sub_ps(_mm_mul_ps(c2.real, c2_conjuge.real), _mm_mul_ps(c2.imaginary, c2_conjuge.imaginary));

    complexe_float_t nominateur = mult_complexe_float(c1, c2_conjuge);

    div.real = _mm_div_ps(nominateur.real, denominateur);
    div.imaginary = _mm_div_ps(nominateur.imaginary, denominateur);

    return div;
}

inline complexe_double_t div_complexe_double(const complexe_double_t c1, const complexe_double_t c2) {
    complexe_double_t div;
    complexe_double_t c2_conjuge = c2;
    c2_conjuge.imaginary = -c2.imaginary;
    float denominateur = c2.real * c2_conjuge.real - c2.imaginary * c2_conjuge.imaginary;
    complexe_double_t nominateur = mult_complexe_double(c1, c2_conjuge);
    div.real = nominateur.real / denominateur;
    div.imaginary = nominateur.imaginary / denominateur;

    return div;
}
