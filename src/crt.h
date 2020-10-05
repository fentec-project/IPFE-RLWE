#ifndef CRT_H
#define CRT_H

#include <gmp.h>

void crt_convert (const uint32_t a[SIFE_N], uint32_t a_crt[SIFE_NMODULI][SIFE_N]);

void crt_convert_gmp (const mpz_t a[SIFE_N], uint32_t a_crt[SIFE_NMODULI][SIFE_N]);

void crt_convert_generic (const uint32_t a[SIFE_L], uint32_t a_crt[SIFE_NMODULI][SIFE_L], const int len);

void crt_reverse_gmp(mpz_t a[SIFE_N], const uint32_t a_crt[SIFE_NMODULI][SIFE_N]);

#endif
