#ifndef POLY_MUL_H
#define POLY_MUL_H


void add_mod (const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod);

void sub_mod (const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod);

void poly_mul_mod (const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod);

void poly_mul_mac_mod (const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod);

#endif