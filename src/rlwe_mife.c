#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include "params.h"
#include "rlwe_mife.h"
#include "crt.h"
#include "sample.h"
#include "randombytes.h"
#include "ntt.h"
#include "arith_rns.h"
#include "gauss.h"
#include "aes256ctr.h"


void rlwe_mife_setup(uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N]) 
{
	int i, j, k;
	uint32_t e_crt[MIFE_NMODULI][MIFE_N];
	uint32_t msk_ntt[MIFE_N];

	aes256ctr_ctx state_secret, state_error;

	unsigned char seed[32];

	randombytes(seed, 32);

	sample_polya(seed, mpk[MIFE_L]);

	randombytes(seed, 32);
	aes256ctr_init(&state_error, seed, 0);

	// Store a in NTT domain
	for (i = 0; i < MIFE_NMODULI; ++i) {
		CT_forward(mpk[MIFE_L][i], i);
	}

	// Sample s_i and e_i with i = 1...l from D_sigma1
	// pk_i = a * s_i + e_i
	for (i = 0; i < MIFE_L; ++i) {
		gaussian_sampler_S1(&state_secret, msk[i], MIFE_N);
		gaussian_sampler_S1(&state_error, e_crt, MIFE_N);
		for (j = 0; j < MIFE_NMODULI; ++j) {
			// Store pk_i in NTT domain but not s_i
			for (k = 0; k < MIFE_N; ++k) {
				msk_ntt[k] = msk[i][j][k];
			}
			CT_forward(msk_ntt, j);
			CT_forward(e_crt[j], j);
			point_mul(mpk[MIFE_L][j], msk_ntt, mpk[i][j], j);
			//poly_mul_ntt(mpk[MIFE_L][j], msk[i][j], mpk[i][j], j);
			poly_add_mod(mpk[i][j], e_crt[j], mpk[i][j], j);
		}
	}

}

void rlwe_mife_encrypt(uint32_t m[MIFE_L], uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N]) //add const keywords
{
	int i, j, k;
	uint32_t r_crt[MIFE_NMODULI][MIFE_N], f_crt[MIFE_NMODULI][MIFE_N];

	uint32_t m_crt[MIFE_NMODULI][MIFE_L];

	aes256ctr_ctx state_s2, state_s3;
	unsigned char seed[32];

	// CRT and scaled message
	crt_convert_generic(m, m_crt, MIFE_L); // needs to be changed. messagges are small no need for reduction

	uint64_t mxm;
	for (i = 0; i < MIFE_L; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			mxm = (uint64_t)m_crt[j][i] * MIFE_SCALE_M_MOD_Q_I[j];
			m_crt[j][i] = mod_prime(mxm, j);
		}
	}

	randombytes(seed, 32);
	aes256ctr_init(&state_s2, seed, 0);

	// Sample r, f_0 from D_sigma2

	gaussian_sampler_S2(&state_s2, r_crt, MIFE_N);
	gaussian_sampler_S2(&state_s2, f_crt, MIFE_N);

	// r in NTT domain
	for (i = 0; i < MIFE_NMODULI; ++i) {
		CT_forward(r_crt[i], i);
	}

	for (i = 0; i < MIFE_NMODULI; ++i) {
		//poly_mul_ntt(mpk[MIFE_L][i], r_crt[i], c[MIFE_L][i], i);
		point_mul(mpk[MIFE_L][i], r_crt[i], c[MIFE_L][i], i);
		GS_reverse(c[MIFE_L][i], i);
		poly_add_mod(c[MIFE_L][i], f_crt[i], c[MIFE_L][i], i);
	}

	// Sample f_i with i = 1...l from D_sigma3
	// c_i = pk_i * r + f_i + (floor(q/p)m_i)1_R
	for (i = 0; i < MIFE_L; ++i) {
		gaussian_sampler_S3(&state_s3, f_crt, MIFE_N);
		
		for (j = 0; j < MIFE_NMODULI; ++j) {
			//poly_mul_ntt(mpk[i][j], r_crt[j], c[i][j], j);
			point_mul(mpk[i][j], r_crt[j], c[i][j], j);
			GS_reverse(c[i][j], j);
			poly_add_mod(c[i][j], f_crt[j], c[i][j], j);
			for (k = 0; k < MIFE_N; ++k) {
				c[i][j][k] = mod_prime((c[i][j][k] + m_crt[j][i]), j);
			}
		}
	}
}

void rlwe_mife_keygen(const uint32_t y[MIFE_L], const uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N], uint32_t sk_y[MIFE_NMODULI][MIFE_N])
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
			for (k = 0; k < MIFE_N; ++k) {
				mac = (uint64_t)y_crt[j][i]*msk[i][j][k];
				mac = mac + sk_y[j][k];
				sk_y[j][k] = mod_prime(mac, j);
			}
		}
	}
}


void rlwe_mife_decrypt_gmp(uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N], const uint32_t y[MIFE_L], uint32_t sk_y[MIFE_NMODULI][MIFE_N], mpz_t dy[MIFE_N])
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
			for (k = 0; k < MIFE_N; ++k) {
				mac = (uint64_t)y_crt[j][i]*c[i][j][k];
				mac = mac + d_y[j][k];
				d_y[j][k] = mod_prime(mac, j);
			}
		}
	}

	for (i = 0; i < MIFE_NMODULI; ++i) {
		poly_mul_ntt(c[MIFE_L][i], sk_y[i], c0sy[i], i);
		poly_sub_mod(d_y[i], c0sy[i], d_y[i], i);
	}

	//crt_reverse(dy, d_y);
	crt_reverse_gmp(dy, d_y);
}

void rlwe_mife_encrypt_vec(uint32_t m[MIFE_L][MIFE_N], uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N], uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N]) //add const keywords
{
	int i, j, k;
	uint32_t r_crt[MIFE_NMODULI][MIFE_N], f_crt[MIFE_NMODULI][MIFE_N];

	uint32_t m_crt[MIFE_L][MIFE_NMODULI][MIFE_N];

	aes256ctr_ctx state_s2, state_s3;
	unsigned char seed[32];

	// CRT and scaled message
	uint64_t mxm;
	for (i = 0; i < MIFE_L; ++i) {
		crt_convert(m[i], m_crt[i]);
		for (j = 0; j < MIFE_NMODULI; ++j) {
			for (k = 0; k < MIFE_N; ++k) {
				mxm = (uint64_t)m_crt[i][j][k] * MIFE_SCALE_M_MOD_Q_I[j];
				m_crt[i][j][k] = mod_prime(mxm, j);
			}
		}
	}

	randombytes(seed, 32);
	aes256ctr_init(&state_s2, seed, 0);

	gaussian_sampler_S2(&state_s2, r_crt, MIFE_N);
	gaussian_sampler_S2(&state_s2, f_crt, MIFE_N);

	// r in NTT domain
	for (i = 0; i < MIFE_NMODULI; ++i) {
		CT_forward(r_crt[i], i);
	}

	for (i = 0; i < MIFE_NMODULI; ++i) {
		//poly_mul_ntt(mpk[MIFE_L][i], r_crt[i], c[MIFE_L][i], i);
		point_mul(mpk[MIFE_L][i], r_crt[i], c[MIFE_L][i], i);
		GS_reverse(c[MIFE_L][i], i);
		poly_add_mod(c[MIFE_L][i], f_crt[i], c[MIFE_L][i], i);
	}

	// Sample f_i with i = 1...l from D_sigma3
	// c_i = pk_i * r + f_i + (floor(q/p)m_i)1_R
	for (i = 0; i < MIFE_L; ++i) {
		gaussian_sampler_S3(&state_s3, f_crt, MIFE_N);
		
		for (j = 0; j < MIFE_NMODULI; ++j) {
			//poly_mul_ntt(mpk[i][j], r_crt[j], c[i][j], j);
			point_mul(mpk[i][j], r_crt[j], c[i][j], j);
			GS_reverse(c[i][j], j);
			poly_add_mod(c[i][j], f_crt[j], c[i][j], j);
			poly_add_mod(c[i][j], m_crt[i][j], c[i][j], j);
		}
	}
}


void rlwe_mife_decrypt_gmp_vec(uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N], const uint32_t y[MIFE_L], uint32_t sk_y[MIFE_NMODULI][MIFE_N], mpz_t dy[MIFE_N])
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
			for (k = 0; k < MIFE_N; ++k) {
				mac = (uint64_t)y_crt[j][i]*c[i][j][k];
				mac = mac + d_y[j][k];
				d_y[j][k] = mod_prime(mac, j);
			}
		}
	}

	for (i = 0; i < MIFE_NMODULI; ++i) {
		poly_mul_ntt(c[MIFE_L][i], sk_y[i], c0sy[i], i);
		poly_sub_mod(d_y[i], c0sy[i], d_y[i], i);
	}

	//crt_reverse(dy, d_y);
	crt_reverse_gmp(dy, d_y);
}



void round_extract_gmp(mpz_t a[MIFE_N])
{
	int i;

	mpz_t quotient, rem, a_i, MIFE_Q_gmp, MIFE_Q_gmp_by2, MIFE_P_gmp;
	mpz_init(quotient);
	mpz_init(rem);
	mpz_init(a_i);
	mpz_init(MIFE_Q_gmp);
	mpz_init(MIFE_P_gmp);
	mpz_init(MIFE_Q_gmp_by2);

	if(mpz_set_str(MIFE_Q_gmp, MIFE_Q_str, 10)!=0){
		printf("--ERROR unable to set Q to gmp--\n");
		return;
	}

	if(mpz_set_str(MIFE_P_gmp, MIFE_P_str, 10)!=0){
		printf("--ERROR unable to set P to gmp--\n");
		return;
	}

	mpz_fdiv_q_ui(MIFE_Q_gmp_by2, MIFE_Q_gmp, 2);
	//gmp_printf("d[0]: %Zd\n", a[0]);
	for (i = 0; i < MIFE_N; ++i) {
		//mpz_set_ui(a_i, a[i]);
		mpz_set(a_i, a[i]);
		//mpz_mul_ui(a_i, a_i, MIFE_P);
		mpz_mul(a_i, a_i, MIFE_P_gmp);
		//mpz_fdiv_qr_ui(quotient, rem, a_i, MIFE_Q);
		mpz_fdiv_qr(quotient, rem, a_i, MIFE_Q_gmp);
		//if( mpz_cmp_ui(rem, (MIFE_Q_gmp >> 1)) > 0 ) {
		if( mpz_cmp(rem, MIFE_Q_gmp_by2 ) > 0 ) {
			mpz_add_ui(quotient, quotient, 1);
		}
		//a[i] = mpz_get_ui(quotient);
		mpz_set(a[i], quotient);
	}
	mpz_clear(quotient);
	mpz_clear(rem);
	mpz_clear(a_i);
	mpz_clear(MIFE_Q_gmp);
	mpz_clear(MIFE_P_gmp);
	mpz_clear(MIFE_Q_gmp_by2);
}
