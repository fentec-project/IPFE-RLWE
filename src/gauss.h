#ifndef GAUSS_H
#define GAUSS_H

#include "aes256ctr.h"
#include "params.h"

/*
void gaussian_sampler_S1(unsigned char *seed, int64_t sample[SIFE_NMODULI][SIFE_N], uint32_t slen);

void gaussian_sampler_S2(unsigned char *seed, int64_t sample[SIFE_NMODULI][SIFE_N], uint32_t slen);

void gaussian_sampler_S3(unsigned char *seed, int64_t sample[SIFE_NMODULI][SIFE_N], uint32_t slen);
*/

void gaussian_sampler_S1(aes256ctr_ctx *state, uint32_t sample[SIFE_NMODULI][SIFE_N], uint32_t slen);

void gaussian_sampler_S2(aes256ctr_ctx *state, uint32_t sample[SIFE_NMODULI][SIFE_N], uint32_t slen);

void gaussian_sampler_S3(aes256ctr_ctx *state, uint32_t sample[SIFE_NMODULI][SIFE_N], uint32_t slen);

#endif
