#ifndef RLWE_MIFE_H
#define RLWE_MIFE_H

#include "params.h"

void rlwe_mife_setup (uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N]);

void rlwe_mife_encrypt (uint32_t m[MIFE_L], const uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N]);

void rlwe_mife_keygen (const uint32_t y[MIFE_L], const uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N], uint32_t sk_y[MIFE_NMODULI][MIFE_N]);

void rlwe_mife_decrypt (const uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N], const uint32_t y[MIFE_L], const uint32_t sk_y[MIFE_NMODULI][MIFE_N], uint64_t dy[MIFE_N]);

void rlwe_mife_encrypt_vec (uint32_t m[MIFE_L][MIFE_N], const uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N]);

void rlwe_mife_keygen_vec (const uint32_t y[MIFE_L][MIFE_N], const uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N], uint32_t sk_y[MIFE_NMODULI][MIFE_N]);

void rlwe_mife_decrypt_vec (const uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N], const uint32_t y[MIFE_L][MIFE_N], const uint32_t sk_y[MIFE_NMODULI][MIFE_N], uint64_t dy[MIFE_N]);

void round_extract (uint64_t a[MIFE_N]);

#endif