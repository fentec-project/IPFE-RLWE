#ifndef ARITH_RNS_H
#define ARITH_RNS_H


uint32_t add_mod_ntt (uint32_t a, uint32_t b, uint32_t sel);
uint32_t sub_mod_ntt (uint32_t a, uint32_t b, uint32_t sel);
uint32_t mul_mod_ntt (uint32_t a, uint32_t b, uint32_t sel);

// Careful: sel is the index for the array of moduli, not the modulo
uint32_t mod_prime (uint64_t m, uint32_t sel);

void poly_add_mod (const uint32_t a[SIFE_N], const uint32_t b[SIFE_N], uint32_t c[SIFE_N], const uint32_t sel);

void poly_sub_mod (const uint32_t a[SIFE_N], const uint32_t b[SIFE_N], uint32_t c[SIFE_N], const uint32_t sel);

#endif