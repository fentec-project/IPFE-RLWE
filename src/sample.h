#ifndef SAMPLE_H
#define SAMPLE_H

#include "params.h"
#include <gmp.h>

void sample_zeros (uint64_t a[SIFE_N]);

void sample_ones (uint64_t a[SIFE_N]);

//void sample_uniform (uint64_t a[SIFE_N]);
void sample_uniform(mpz_t a[SIFE_N]);

void sample_polya(unsigned char *seed, uint32_t poly_a[SIFE_NMODULI][SIFE_N]);


void sample_m(uint32_t a[SIFE_N]);

void sample_x(uint32_t a[SIFE_L]);

void sample_y(uint32_t a[SIFE_L]);

#endif
