#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "rlwe_mife.h"
#include "crt.h"
#include "poly_mul.h"
#include "sample.h"

void
rlwe_mife_setup
(uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N])
{
	int i, j;
	uint64_t a[MIFE_N], s[MIFE_N], e[MIFE_N];
	uint32_t e_crt[MIFE_NMODULI][MIFE_N];

	// Sample a from U(0,q-1)
	sample_uniform(a);
	crt_convert(a, mpk[MIFE_L]);

	// Sample s_i and e_i with i = 1...l from D_sigma1
	// pk_i = a * s_i + e_i
	for (i = 0; i < MIFE_L; ++i) {
		sample_sigma1(s);
		sample_sigma1(e);
		crt_convert(s, msk[i]);
		crt_convert(e, e_crt);
		for (j = 0; j < MIFE_NMODULI; ++j) {
			poly_mul_mod(mpk[MIFE_L][j], msk[i][j], mpk[i][j], MIFE_MOD_Q_I[j]);
			add_mod(mpk[i][j], e_crt[j], mpk[i][j], MIFE_MOD_Q_I[j]);
		}
	}
}

void
rlwe_mife_encrypt
(uint32_t m[MIFE_L][MIFE_NMODULI][MIFE_N], const uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N])
{
	int i, j, k;
	uint64_t r[MIFE_N], f[MIFE_N];
	uint32_t r_crt[MIFE_NMODULI][MIFE_N], f_crt[MIFE_NMODULI][MIFE_N];

/**	//Message (m1...ml) with coeffs in (-B,B)
	//XXX - I AM ASSUMING IT'S PASSED ALREADY IN CRT DOMAIN, BUT NOT SCALED
	//XXX - IS IT SCALED LIKE THIS??
	for (i = 0; i < MIFE_L; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			for (k = 0; k < MIFE_N; ++k) {
				m[i][j][k] = (MIFE_SCALE_M * m[i][j][k]) % MIFE_MOD_Q_I[j];
			}
		}
	}/**/

	// Sample r, f_0 from D_sigma2
	sample_sigma2(r);
	sample_sigma2(f);
	// c_0 = a * r + f_0
	crt_convert(r, r_crt);
	crt_convert(f, f_crt);
	for (i = 0; i < MIFE_NMODULI; ++i) {
		poly_mul_mod(mpk[MIFE_L][i], r_crt[i], c[MIFE_L][i], MIFE_MOD_Q_I[i]);
		add_mod(c[MIFE_L][i], f_crt[i], c[MIFE_L][i], MIFE_MOD_Q_I[i]);
	}

	// Sample f_i with i = 1...l from D_sigma3
	// c_i = pk_i * r + f_i + (floor(q/p)m_i)1_R
	for (i = 0; i < MIFE_L; ++i) {
		sample_sigma3(f);
		crt_convert(f, f_crt);
		for (j = 0; j < MIFE_NMODULI; ++j) {
			poly_mul_mod(mpk[i][j], r_crt[j], c[i][j], MIFE_MOD_Q_I[j]);
			add_mod(c[i][j], f_crt[j], c[i][j], MIFE_MOD_Q_I[j]);
			add_mod(c[i][j], m[i][j], c[i][j], MIFE_MOD_Q_I[j]);
		}
	}
}

void
rlwe_mife_keygen
(const uint32_t y[MIFE_L][MIFE_NMODULI][MIFE_N], const uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N], uint32_t sk_y[MIFE_NMODULI][MIFE_N])
{
	int i, j;
	//For a key that decrypts m·y
	// s_y = sum_i=1...l( y_i * s_i )
	for (i = 0; i < MIFE_NMODULI; ++i) {
		for (j = 0; j < MIFE_N; ++j) {
			sk_y[i][j] = 0;
		}
	}
	for (i = 0; i < MIFE_L; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			poly_mul_mac_mod(y[i][j], msk[i][j], sk_y[j], MIFE_MOD_Q_I[j]);
		}
	}
}

void
rlwe_mife_decrypt
(const uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N], const uint32_t y[MIFE_L][MIFE_NMODULI][MIFE_N], const uint32_t sk_y[MIFE_NMODULI][MIFE_N], uint32_t d_y[MIFE_NMODULI][MIFE_N])
{
	int i, j;
	uint32_t c0sy[MIFE_NMODULI][MIFE_N];
	// d = ( sum_i=1...l ( y_i * c_i ) ) - c_0 * s_y [to get m·y(q/p)1_R]
	for (i = 0; i < MIFE_NMODULI; ++i) {
		for (j = 0; j < MIFE_N; ++j) {
			d_y[i][j] = 0;
		}
	}
	for (i = 0; i < MIFE_L; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			poly_mul_mac_mod(y[i][j], c[i][j], d_y[j], MIFE_MOD_Q_I[j]);
		}
	}
	for (i = 0; i < MIFE_NMODULI; ++i) {
		poly_mul_mod(c[MIFE_L][i], sk_y[i], c0sy[i], MIFE_MOD_Q_I[i]);
		sub_mod(d_y[i], c0sy[i], d_y[i], MIFE_MOD_Q_I[i]);
	}
}
