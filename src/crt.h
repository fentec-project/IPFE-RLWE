#ifndef CRT_H
#define CRT_H

#include <gmp.h>

void crt_convert (const uint64_t a[MIFE_N], uint32_t a_crt[MIFE_NMODULI][MIFE_N]);

void crt_convert_generic (const uint32_t a[MIFE_L], uint32_t a_crt[MIFE_NMODULI][MIFE_L], const int len);

void crt_reverse_gmp(mpz_t a[MIFE_N], const uint32_t a_crt[MIFE_NMODULI][MIFE_N]);

#endif
