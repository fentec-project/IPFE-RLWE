#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "modred.h"
#include "params.h"
#include "sample.h"
#include <stdlib.h>

//#define UINT64_RAND_MAX ((RAND_MAX << 32)|RAND_MAX)

void
sample_ones
(uint64_t a[MIFE_N])
{
	int i;
	for (i = 0; i < MIFE_N; ++i) {
		a[i] = 1;
	}
}

void
sample_uniform
(uint64_t a[MIFE_N])
{
	int i;
	uint64_t num;
	for (i = 0; i < MIFE_N; ++i) {
		do {
			num = rand();
			num = (num << 32) | (rand());
		//} while (num >= (UINT64_RAND_MAX - (UINT64_RAND_MAX % MIFE_Q)));
		//a[i] = num % MIFE_Q;
		} while (num >= MIFE_Q);
		a[i] = num;
	}
}

void
sample_sigma1
(uint64_t a[MIFE_N])
{
	sample_uniform(a);
}

void
sample_sigma2
(uint64_t a[MIFE_N])
{
	sample_uniform(a);
}

void
sample_sigma3
(uint64_t a[MIFE_N])
{
	sample_uniform(a);
}
