#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "rlwe_mife.h"
#include "crt.h"
#include "poly_mul.h"
#include "sample.h"

#include <gmp.h>

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
		//sample_zeros(e);
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
(uint32_t m[MIFE_L], const uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N])
{
	int i, j, k;
	uint64_t r[MIFE_N], f[MIFE_N];
	uint32_t r_crt[MIFE_NMODULI][MIFE_N], f_crt[MIFE_NMODULI][MIFE_N];

	uint32_t m_crt[MIFE_NMODULI][MIFE_L];

	// CRT and scaled message
	crt_convert_generic(m, m_crt, MIFE_L);

	uint64_t mxm;
	for (i = 0; i < MIFE_L; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			mxm = (uint64_t)m_crt[j][i] * MIFE_SCALE_M_MOD_Q_I[j];
			m_crt[j][i] = mod_red(mxm, MIFE_MOD_Q_I[j]);
		}
	}

	// Sample r, f_0 from D_sigma2
	sample_sigma2(r);
	sample_sigma2(f);
	//sample_zeros(f);
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
		//sample_zeros(f);
		crt_convert(f, f_crt);
		for (j = 0; j < MIFE_NMODULI; ++j) {
			poly_mul_mod(mpk[i][j], r_crt[j], c[i][j], MIFE_MOD_Q_I[j]);
			add_mod(c[i][j], f_crt[j], c[i][j], MIFE_MOD_Q_I[j]);
			//add_mod(c[i][j], m[i][j], c[i][j], MIFE_MOD_Q_I[j]);
			for (int k = 0; k < MIFE_N; ++k) {
				c[i][j][k] = mod_red((c[i][j][k] + m_crt[j][i]), MIFE_MOD_Q_I[j]);
			}
		}
	}
}

void
rlwe_mife_keygen
(const uint32_t y[MIFE_L], const uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N], uint32_t sk_y[MIFE_NMODULI][MIFE_N])
{
	int i, j, k;
	uint64_t mac;

	uint32_t y_crt[MIFE_NMODULI][MIFE_L];

	crt_convert_generic(y, y_crt, MIFE_L);

	for (i = 0; i < MIFE_NMODULI; ++i) {
		for (j = 0; j < MIFE_N; ++j) {
			sk_y[i][j] = 0;
		}
	}

	for (i = 0; i < MIFE_L; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			//poly_mul_mac_mod(y[i][j], msk[i][j], sk_y[j], MIFE_MOD_Q_I[j]);
			for (k = 0; k < MIFE_N; ++k) {
				mac = (uint64_t)y_crt[j][i]*msk[i][j][k];
				mac = mac + sk_y[j][k];
				sk_y[j][k] = mod_red(mac, MIFE_MOD_Q_I[j]);
			}
		}
	}
}

void
rlwe_mife_decrypt
(const uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N], const uint32_t y[MIFE_L], const uint32_t sk_y[MIFE_NMODULI][MIFE_N], uint64_t dy[MIFE_N])
{
	int i, j, k;
	uint64_t mac;

	uint32_t c0sy[MIFE_NMODULI][MIFE_N];
	uint32_t d_y[MIFE_NMODULI][MIFE_N];

	uint32_t y_crt[MIFE_NMODULI][MIFE_L];

	crt_convert_generic(y, y_crt, MIFE_L);

	for (i = 0; i < MIFE_NMODULI; ++i) {
		for (j = 0; j < MIFE_N; ++j) {
			d_y[i][j] = 0;
		}
	}

	for (i = 0; i < MIFE_L; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			//poly_mul_mac_mod(y[i][j], c[i][j], d_y[j], MIFE_MOD_Q_I[j]);
			for (k = 0; k < MIFE_N; ++k) {
				mac = (uint64_t)y_crt[j][i]*c[i][j][k];
				mac = mac + d_y[j][k];
				d_y[j][k] = mod_red(mac, MIFE_MOD_Q_I[j]);
			}
		}
	}

	for (i = 0; i < MIFE_NMODULI; ++i) {
		poly_mul_mod(c[MIFE_L][i], sk_y[i], c0sy[i], MIFE_MOD_Q_I[i]);
		sub_mod(d_y[i], c0sy[i], d_y[i], MIFE_MOD_Q_I[i]);
	}

	crt_reverse(dy, d_y);
}

void
round_extract
(uint64_t a[MIFE_N])
{
	int i;

	mpz_t quotient, rem, a_i;
	mpz_init(quotient);
	mpz_init(rem);
	mpz_init(a_i);

	for (i = 0; i < MIFE_N; ++i) {
		mpz_set_ui(a_i, a[i]);
		mpz_mul_ui(a_i, a_i, MIFE_P);
		mpz_fdiv_qr_ui(quotient, rem, a_i, MIFE_Q);
		if( mpz_cmp_ui(rem, (MIFE_Q >> 1)) > 0 ) {
			mpz_add_ui(quotient, quotient, 1);
		}
		a[i] = mpz_get_ui(quotient);
	}
	mpz_clear(quotient);
	mpz_clear(rem);
	mpz_clear(a_i);

	//for (i = 0; i < MIFE_N; ++i) {
	//	a[i] = (MIFE_P*a[i])/MIFE_Q;
	//}
}
