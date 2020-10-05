#include <stdint.h>
#include "params.h"
#include "arith_rns.h"
#include "crt.h"


void
crt_convert
(const uint64_t a[MIFE_N], uint32_t a_crt[MIFE_NMODULI][MIFE_N])
{
	int i, j;
	for (i = 0; i < MIFE_N; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			a_crt[j][i] = mod_prime(a[i], j);
		}
	}
}

void
crt_convert_gmp
(const mpz_t a[MIFE_N], uint32_t a_crt[MIFE_NMODULI][MIFE_N])
{
	int i, j;
	for (i = 0; i < MIFE_N; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			a_crt[j][i]= mpz_fdiv_ui(a[i], MIFE_MOD_Q_I[j]);

		}
	}
}


void
crt_convert_generic
(const uint32_t a[MIFE_L], uint32_t a_crt[MIFE_NMODULI][MIFE_L], const int len)//fix. len not needed
{
	int i, j;
	for (i = 0; i < len; ++i) {
		for (j = 0; j < MIFE_NMODULI; ++j) {
			a_crt[j][i] = mod_prime(a[i], j);
		}
	}
}

void
crt_reverse_gmp
(mpz_t a[MIFE_N], const uint32_t split_a[MIFE_NMODULI][MIFE_N])
{
	uint64_t i,j,k;
	mpz_t gmp_u, gmp_c, gmp_x;

	mpz_init(gmp_u);
	mpz_init(gmp_c);
	mpz_init(gmp_x);

	for (k = 0; k < MIFE_N; k++) {
		mpz_set_ui(gmp_u, split_a[0][k]);
		mpz_set(gmp_x, gmp_u);
		for (i = 1; i <= MIFE_NMODULI-1; i++) {
			mpz_set_ui(gmp_c, split_a[i][k]);	//c=vi
			mpz_sub(gmp_u, gmp_c, gmp_x);	//vi-x
			mpz_mul_ui(gmp_u, gmp_u, MIFE_CRT_CONSTS[i]);//(vi-x)Ci
			mpz_mod_ui(gmp_u, gmp_u, MIFE_MOD_Q_I[i]);//(vi-x)Ci mod mi
			mpz_set_str(gmp_c, "1", 10);
			for(j = 0; j <= i-1; j++) {
				//mpz_mul(gmp_c, gmp_c, gmp_q_ar[j]);
				mpz_mul_ui(gmp_c, gmp_c, MIFE_MOD_Q_I[j]);
			} 
			mpz_mul( gmp_u, gmp_u, gmp_c );
			mpz_add(gmp_x, gmp_x, gmp_u);				
		}
		mpz_set(a[k],gmp_x);
	}

	mpz_clear(gmp_u);
	mpz_clear(gmp_c);
	mpz_clear(gmp_x);
}
