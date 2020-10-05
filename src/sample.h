#ifndef SAMPLE_H
#define SAMPLE_H

#include "params.h"
#include <gmp.h>

void sample_zeros (uint64_t a[MIFE_N]);

void sample_ones (uint64_t a[MIFE_N]);

//void sample_uniform (uint64_t a[MIFE_N]);
void sample_uniform(mpz_t a[MIFE_N]);

void sample_polya(unsigned char *seed, uint32_t poly_a[MIFE_NMODULI][MIFE_N]);

void sample_sigma1 (uint64_t a[MIFE_N]);
void sample_sigma2 (uint64_t a[MIFE_N]);
void sample_sigma3 (uint64_t a[MIFE_N]);

void sample_message (uint64_t a[MIFE_L]);

void sample_x(uint32_t a[MIFE_L]);
void sample_y(uint32_t a[MIFE_L]);

#endif
