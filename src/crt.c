#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "modred.h"
#include "crt.h"

void
crt_convert
(const uint64_t a[MIFE_N], uint32_t a_crt[MIFE_NMODULI][MIFE_N])
{
	int i, j;
	for (i = 0; i < MIFE_N; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			//a_crt[j][i] = a[i] % MIFE_MOD_Q_I[j];
			a_crt[j][i] = mod_red(a[i], MIFE_MOD_Q_I[j]);
		}
	}
}

void
crt_convert_generic
(const uint32_t a[MIFE_L], uint32_t a_crt[MIFE_NMODULI][MIFE_L], const int len)
{
	int i, j;
	for (i = 0; i < len; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			//a_crt[j][i] = a[i] % MIFE_MOD_Q_I[j];
			a_crt[j][i] = mod_red(a[i], MIFE_MOD_Q_I[j]);
		}
	}
}

static
int64_t
int_mul_mod
(uint64_t a, uint64_t b, uint64_t m)
{
	//I haven't debugged the algorithm, I took it from here and assume it works
	//https://math.stackexchange.com/questions/33923/modular-multiplication-with-machine-word-limitations

	// If the problem is here, we could also try another algorithm as in
	// https://stackoverflow.com/questions/2566010/fastest-way-to-calculate-a-128-bit-integer-modulo-a-64-bit-integer?rq=1

	if (a >= m) a %= m;
	if (b >= m) b %= m;

	long double x = a;

	uint64_t c = x * b / m;
	int64_t r = (int64_t)(a * b - c * m) % (int64_t)m;
	return r < 0 ? r + m : r;
}

void
crt_reverse__old
(uint64_t a[MIFE_N], const uint32_t a_crt[MIFE_NMODULI][MIFE_N])
{
	int i, j;
	uint64_t mac;
	for (i = 0; i < MIFE_N; ++i) {
		a[i] = 0;
	}
	for (i = 0; i < MIFE_N; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			//a[i] = (a[i] + a_crt[j][i] * MIFE_CRT_REVERSE_I[j]) % MIFE_Q;
			a[i] = (a[i] + (uint64_t)int_mul_mod(a_crt[j][i], MIFE_CRT_REVERSE_I[j], MIFE_Q)) % MIFE_Q;

			// XXX - problem is mod_red takes uint32_t modulus

			//mac = (uint64_t)int_mul_mod(a_crt[j][i], MIFE_CRT_REVERSE_I[j], MIFE_Q);
			//mac = mac + a[i];
			//a[i] = mod_red(mac, MIFE_Q);
		}
	}
}


//Ad-hoc CRT reverse transformation
/*Algo from http://cacr.uwaterloo.ca/hac/about/chap14.pdf, page 612 unrolled for two primes*/
uint64_t
inverse_CRT
(uint32_t v1, uint32_t v2) //only 2 primes since we are testing for 64 bits
{
	uint64_t u, x;

	u = v1;
	x = (uint64_t)u;
	u = v2 - x;
	u = u * MIFE_C2;
	u = mod_red(u, MIFE_Q2);
	x = x + u*MIFE_Q1;

	return x;
}

//TODO: make inverse_CRT generic

void
crt_reverse
(uint64_t a[MIFE_N], const uint32_t a_crt[MIFE_NMODULI][MIFE_N])
{
	int i, j;

	for (i = 0; i < MIFE_N; ++i) {
		a[i] = 0;
	}
	for (i = 0; i < MIFE_N; ++i) {
		// XXX - This function is not generic anymore:
		// It assumes modulus is broken into two
		a[i] = inverse_CRT(a_crt[0][i], a_crt[1][i]);
	}
}
