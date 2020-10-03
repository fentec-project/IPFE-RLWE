#include <stdint.h>
#include <stdio.h>
#include <x86intrin.h>
#include "params.h"
#include "arith_rns.h"
#include "randombytes.h"
#include "aes256ctr.h"
#include "gauss.h"


#define CDT_ENTRY_SIZE 16
#define CDT_LOW_MASK 0x7fffffffffffffff
#define CDT_LENGTH 9 /* [0..tau*sigma]=[0..9] */

#define BERNOULLI_ENTRY_SIZE 9 /* 72bit randomness */

/*ceil(DESIRED_SIGMA/sigma_0) sigma_0=sqrt( 1/(2*log(2) ) )*/
#define BINARY_SAMPLER_K_S1 266		//sigma_1
#define BINARY_SAMPLER_K_S2 304214978LL	//sigma_2
#define BINARY_SAMPLER_K_S3 608429953LL	//sigma_3

/* -1/k^2 */
#define BINARY_SAMPLER_K_2_INV_S1 (-1.0/(BINARY_SAMPLER_K_S1 * BINARY_SAMPLER_K_S1))
#define BINARY_SAMPLER_K_2_INV_S2 (-1.0/(BINARY_SAMPLER_K_S2 * BINARY_SAMPLER_K_S2))
#define BINARY_SAMPLER_K_2_INV_S3 (-1.0/(BINARY_SAMPLER_K_S3 * BINARY_SAMPLER_K_S3))


#define EXP_MANTISSA_PRECISION 52
#define EXP_MANTISSA_MASK ((1LL << EXP_MANTISSA_PRECISION) - 1)
#define R_MANTISSA_PRECISION (EXP_MANTISSA_PRECISION + 1)
#define R_MANTISSA_MASK ((1LL << R_MANTISSA_PRECISION) - 1)
#define R_EXPONENT_L (8 * BERNOULLI_ENTRY_SIZE - R_MANTISSA_PRECISION)

#define DOUBLE_ONE (1023LL << 52)

#define UNIFORM_SIZE 4
#define UNIFORM_REJ 20
#define BARRETT_BITSHIFT (UNIFORM_SIZE * 8)

#define BARRETT_FACTOR_S1 ((1LL << BARRETT_BITSHIFT) / BINARY_SAMPLER_K_S1)
#define BARRETT_FACTOR_S2 ((1LL << BARRETT_BITSHIFT) / BINARY_SAMPLER_K_S2)
#define BARRETT_FACTOR_S3 ((1LL << BARRETT_BITSHIFT) / BINARY_SAMPLER_K_S3)

#define UNIFORM_Q_S1 (BINARY_SAMPLER_K_S1 * BARRETT_FACTOR_S1)
#define UNIFORM_Q_S2 (BINARY_SAMPLER_K_S2 * BARRETT_FACTOR_S2)
#define UNIFORM_Q_S3 (BINARY_SAMPLER_K_S3 * BARRETT_FACTOR_S3)


#define BASE_TABLE_SIZE (4 * CDT_ENTRY_SIZE)
#define BERNOULLI_TABLE_SIZE (4 * BERNOULLI_ENTRY_SIZE)

/* CDT table */
static const __m256i V_CDT[][2] = {{{2200310400551559144, 2200310400551559144, 2200310400551559144, 2200310400551559144}, {3327841033070651387, 3327841033070651387, 3327841033070651387, 3327841033070651387}},
{{7912151619254726620, 7912151619254726620, 7912151619254726620, 7912151619254726620}, {380075531178589176, 380075531178589176, 380075531178589176, 380075531178589176}},
{{5167367257772081627, 5167367257772081627, 5167367257772081627, 5167367257772081627}, {11604843442081400, 11604843442081400, 11604843442081400, 11604843442081400}},
{{5081592746475748971, 5081592746475748971, 5081592746475748971, 5081592746475748971}, {90134450315532, 90134450315532, 90134450315532, 90134450315532}},
{{6522074513864805092, 6522074513864805092, 6522074513864805092, 6522074513864805092}, {175786317361, 175786317361, 175786317361, 175786317361}},
{{2579734681240182346, 2579734681240182346, 2579734681240182346, 2579734681240182346}, {85801740, 85801740, 85801740, 85801740}},
{{8175784047440310133, 8175784047440310133, 8175784047440310133, 8175784047440310133}, {10472, 10472, 10472, 10472}},
{{2947787991558061753, 2947787991558061753, 2947787991558061753, 2947787991558061753}, {0, 0, 0, 0}},
{{22489665999543, 22489665999543, 22489665999543, 22489665999543}, {0, 0, 0, 0}}};

static const __m256i V_CDT_LOW_MASK = {CDT_LOW_MASK, CDT_LOW_MASK, CDT_LOW_MASK, CDT_LOW_MASK};

static const __m256i V_K_K_K_K_S1 = {BINARY_SAMPLER_K_S1, BINARY_SAMPLER_K_S1, BINARY_SAMPLER_K_S1, BINARY_SAMPLER_K_S1};
static const __m256i V_K_K_K_K_S2 = {BINARY_SAMPLER_K_S2, BINARY_SAMPLER_K_S2, BINARY_SAMPLER_K_S2, BINARY_SAMPLER_K_S2};
static const __m256i V_K_K_K_K_S3 = {BINARY_SAMPLER_K_S3, BINARY_SAMPLER_K_S3, BINARY_SAMPLER_K_S3, BINARY_SAMPLER_K_S3};


/* coefficients of the exp evaluation polynomial */
static const __m256i EXP_COFF[] = {{0x3e833b70ffa2c5d4, 0x3e833b70ffa2c5d4, 0x3e833b70ffa2c5d4, 0x3e833b70ffa2c5d4},
								   {0x3eb4a480fda7e6e1, 0x3eb4a480fda7e6e1, 0x3eb4a480fda7e6e1, 0x3eb4a480fda7e6e1},
								   {0x3ef01b254493363f, 0x3ef01b254493363f, 0x3ef01b254493363f, 0x3ef01b254493363f},
								   {0x3f242e0e0aa273cc, 0x3f242e0e0aa273cc, 0x3f242e0e0aa273cc, 0x3f242e0e0aa273cc},
								   {0x3f55d8a2334ed31b, 0x3f55d8a2334ed31b, 0x3f55d8a2334ed31b, 0x3f55d8a2334ed31b},
								   {0x3f83b2aa56db0f1a, 0x3f83b2aa56db0f1a, 0x3f83b2aa56db0f1a, 0x3f83b2aa56db0f1a},
								   {0x3fac6b08e11fc57e, 0x3fac6b08e11fc57e, 0x3fac6b08e11fc57e, 0x3fac6b08e11fc57e},
								   {0x3fcebfbdff556072, 0x3fcebfbdff556072, 0x3fcebfbdff556072, 0x3fcebfbdff556072},
								   {0x3fe62e42fefa7fe6, 0x3fe62e42fefa7fe6, 0x3fe62e42fefa7fe6, 0x3fe62e42fefa7fe6},
								   {0x3ff0000000000000, 0x3ff0000000000000, 0x3ff0000000000000, 0x3ff0000000000000}};
								   
static const __m256d V_INT64_DOUBLE = {0x0010000000000000, 0x0010000000000000, 0x0010000000000000, 0x0010000000000000};
static const __m256d V_DOUBLE_INT64 = {0x0018000000000000, 0x0018000000000000, 0x0018000000000000, 0x0018000000000000};

static const __m256i V_EXP_MANTISSA_MASK = {EXP_MANTISSA_MASK, EXP_MANTISSA_MASK, EXP_MANTISSA_MASK, EXP_MANTISSA_MASK};
static const __m256i V_RES_MANTISSA = {1LL << EXP_MANTISSA_PRECISION, 1LL << EXP_MANTISSA_PRECISION, 1LL << EXP_MANTISSA_PRECISION, 1LL << EXP_MANTISSA_PRECISION};
static const __m256i V_RES_EXPONENT = {R_EXPONENT_L - 1023 + 1, R_EXPONENT_L - 1023 + 1, R_EXPONENT_L - 1023 + 1, R_EXPONENT_L - 1023 + 1};
static const __m256i V_R_MANTISSA_MASK = {R_MANTISSA_MASK, R_MANTISSA_MASK, R_MANTISSA_MASK, R_MANTISSA_MASK};
static const __m256i V_1 = {1, 1, 1, 1};
static const __m256i V_DOUBLE_ONE = {DOUBLE_ONE, DOUBLE_ONE, DOUBLE_ONE, DOUBLE_ONE};

static const __m256d V_K_2_INV_S1 = {BINARY_SAMPLER_K_2_INV_S1, BINARY_SAMPLER_K_2_INV_S1, BINARY_SAMPLER_K_2_INV_S1, BINARY_SAMPLER_K_2_INV_S1};
static const __m256d V_K_2_INV_S2 = {BINARY_SAMPLER_K_2_INV_S2, BINARY_SAMPLER_K_2_INV_S2, BINARY_SAMPLER_K_2_INV_S2, BINARY_SAMPLER_K_2_INV_S2};
static const __m256d V_K_2_INV_S3 = {BINARY_SAMPLER_K_2_INV_S3, BINARY_SAMPLER_K_2_INV_S3, BINARY_SAMPLER_K_2_INV_S3, BINARY_SAMPLER_K_2_INV_S3};

/* constant time CDT sampler */
static inline __m256i cdt_sampler(unsigned char *r)
{
	__m256i x = _mm256_setzero_si256();
	__m256i r1, r2;
	__m256i r1_lt_cdt0, r2_lt_cdt1;
	__m256i r2_eq_cdt1;
	__m256i b;
	
	uint32_t i;
	
	r1 = _mm256_loadu_si256((__m256i *)r);
	r2 = _mm256_loadu_si256((__m256i *)(r + 32));
	
	r1 = _mm256_and_si256(r1, V_CDT_LOW_MASK);
	r2 = _mm256_and_si256(r2, V_CDT_LOW_MASK);

	for (i = 0; i < CDT_LENGTH; i++)
	{
		r1_lt_cdt0 = _mm256_sub_epi64(r1, V_CDT[i][0]);

		r2_lt_cdt1 = _mm256_sub_epi64(r2, V_CDT[i][1]);
		r2_eq_cdt1 = _mm256_cmpeq_epi64(r2, V_CDT[i][1]);

		b = _mm256_and_si256(r1_lt_cdt0, r2_eq_cdt1);
		b = _mm256_or_si256(b, r2_lt_cdt1);
		b = _mm256_srli_epi64(b, 63);

		x = _mm256_add_epi64(x, b);
	}

	return x;
}

/* constant time Bernoulli sampler
 * we directly compute exp(-x/(2*sigma_0^2)), 
 * since sigma_0=sqrt(1/2ln2), exp(-x/(2*sigma_0^2))=2^(-x/k^2), 
 * we use a polynomial to directly evaluate 2^(-x/k^2) */

static inline void bernoulli_sampler_S1(uint64_t *b, __m256i x, unsigned char *r)
{	
	__m256d vx, vx_1, vx_2, vsum;
	__m256i vt, k, vres, vres_mantissa, vres_exponent, vr_mantissa, vr_exponent, vr_exponent2, vres_eq_1, vr_lt_vres_mantissa, vr_lt_vres_exponent;

	/* 2^x=2^(floor(x)+a)=2^(floor(x))*2^a, where a is in [0,1]
	 * we only evaluate 2^a by using a polynomial */
	x = _mm256_or_si256(x, _mm256_castpd_si256(V_INT64_DOUBLE));
	vx = _mm256_sub_pd(_mm256_castsi256_pd(x), V_INT64_DOUBLE);
	vx = _mm256_mul_pd(vx, V_K_2_INV_S1);
	
	vx_1 = _mm256_floor_pd(vx);
	vx_2 = _mm256_add_pd(vx_1, V_DOUBLE_INT64);
	vt = _mm256_sub_epi64(_mm256_castpd_si256(vx_2), _mm256_castpd_si256(V_DOUBLE_INT64));	
	vt = _mm256_slli_epi64(vt, 52);
	
	/* evaluate 2^a */
	vx_2 = _mm256_sub_pd(vx, vx_1);
	vsum = _mm256_fmadd_pd(_mm256_castsi256_pd(EXP_COFF[0]), vx_2, _mm256_castsi256_pd(EXP_COFF[1]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[2]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[3]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[4]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[5]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[6]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[7]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[8]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[9]));
	
	/* combine to compute 2^x */
	vres = _mm256_add_epi64(vt, _mm256_castpd_si256(vsum));

	/* compute the Bernoulli value */
	vres_mantissa = _mm256_and_si256(vres, V_EXP_MANTISSA_MASK);
	vres_mantissa = _mm256_or_si256(vres_mantissa, V_RES_MANTISSA);
	
	vres_exponent = _mm256_srli_epi64(vres, EXP_MANTISSA_PRECISION);
	vres_exponent = _mm256_add_epi64(vres_exponent, V_RES_EXPONENT);
	vres_exponent = _mm256_sllv_epi64(V_1, vres_exponent);
	
	vr_mantissa = _mm256_loadu_si256((__m256i *)r);
	vr_exponent = _mm256_srli_epi64(vr_mantissa, R_MANTISSA_PRECISION);
	vr_mantissa = _mm256_and_si256(vr_mantissa, V_R_MANTISSA_MASK);
	vr_exponent2 = _mm256_set_epi64x(r[35], r[34], r[33], r[32]);
	vr_exponent2 = _mm256_slli_epi64(vr_exponent2, 64 - R_MANTISSA_PRECISION);
	vr_exponent = _mm256_or_si256(vr_exponent, vr_exponent2);

	/* (res == 1.0) || ((r_mantissa < res_mantissa) && (r_exponent < (1 << res_exponent))) */
	vres_eq_1 = _mm256_cmpeq_epi64(vres, V_DOUBLE_ONE);
	vr_lt_vres_mantissa = _mm256_sub_epi64(vr_mantissa, vres_mantissa);	
	vr_lt_vres_exponent = _mm256_sub_epi64(vr_exponent, vres_exponent);
	
	k = _mm256_and_si256(vr_lt_vres_mantissa, vr_lt_vres_exponent);
	k = _mm256_or_si256(k, vres_eq_1);

	_mm256_store_si256((__m256i *)(b), k);
}


static inline void bernoulli_sampler_S2(uint64_t *b, __m256i x, unsigned char *r)
{	
	__m256d vx, vx_1, vx_2, vsum;
	__m256i vt, k, vres, vres_mantissa, vres_exponent, vr_mantissa, vr_exponent, vr_exponent2, vres_eq_1, vr_lt_vres_mantissa, vr_lt_vres_exponent;

	/* 2^x=2^(floor(x)+a)=2^(floor(x))*2^a, where a is in [0,1]
	 * we only evaluate 2^a by using a polynomial */
	x = _mm256_or_si256(x, _mm256_castpd_si256(V_INT64_DOUBLE));
	vx = _mm256_sub_pd(_mm256_castsi256_pd(x), V_INT64_DOUBLE);
	vx = _mm256_mul_pd(vx, V_K_2_INV_S2);
	
	vx_1 = _mm256_floor_pd(vx);
	vx_2 = _mm256_add_pd(vx_1, V_DOUBLE_INT64);
	vt = _mm256_sub_epi64(_mm256_castpd_si256(vx_2), _mm256_castpd_si256(V_DOUBLE_INT64));	
	vt = _mm256_slli_epi64(vt, 52);
	
	/* evaluate 2^a */
	vx_2 = _mm256_sub_pd(vx, vx_1);
	vsum = _mm256_fmadd_pd(_mm256_castsi256_pd(EXP_COFF[0]), vx_2, _mm256_castsi256_pd(EXP_COFF[1]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[2]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[3]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[4]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[5]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[6]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[7]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[8]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[9]));
	
	/* combine to compute 2^x */
	vres = _mm256_add_epi64(vt, _mm256_castpd_si256(vsum));

	/* compute the Bernoulli value */
	vres_mantissa = _mm256_and_si256(vres, V_EXP_MANTISSA_MASK);
	vres_mantissa = _mm256_or_si256(vres_mantissa, V_RES_MANTISSA);
	
	vres_exponent = _mm256_srli_epi64(vres, EXP_MANTISSA_PRECISION);
	vres_exponent = _mm256_add_epi64(vres_exponent, V_RES_EXPONENT);
	vres_exponent = _mm256_sllv_epi64(V_1, vres_exponent);
	
	vr_mantissa = _mm256_loadu_si256((__m256i *)r);
	vr_exponent = _mm256_srli_epi64(vr_mantissa, R_MANTISSA_PRECISION);
	vr_mantissa = _mm256_and_si256(vr_mantissa, V_R_MANTISSA_MASK);
	vr_exponent2 = _mm256_set_epi64x(r[35], r[34], r[33], r[32]);
	vr_exponent2 = _mm256_slli_epi64(vr_exponent2, 64 - R_MANTISSA_PRECISION);
	vr_exponent = _mm256_or_si256(vr_exponent, vr_exponent2);

	/* (res == 1.0) || ((r_mantissa < res_mantissa) && (r_exponent < (1 << res_exponent))) */
	vres_eq_1 = _mm256_cmpeq_epi64(vres, V_DOUBLE_ONE);
	vr_lt_vres_mantissa = _mm256_sub_epi64(vr_mantissa, vres_mantissa);	
	vr_lt_vres_exponent = _mm256_sub_epi64(vr_exponent, vres_exponent);
	
	k = _mm256_and_si256(vr_lt_vres_mantissa, vr_lt_vres_exponent);
	k = _mm256_or_si256(k, vres_eq_1);

	_mm256_store_si256((__m256i *)(b), k);
}

static inline void bernoulli_sampler_S3(uint64_t *b, __m256i x, unsigned char *r)
{	
	__m256d vx, vx_1, vx_2, vsum;
	__m256i vt, k, vres, vres_mantissa, vres_exponent, vr_mantissa, vr_exponent, vr_exponent2, vres_eq_1, vr_lt_vres_mantissa, vr_lt_vres_exponent;

	/* 2^x=2^(floor(x)+a)=2^(floor(x))*2^a, where a is in [0,1]
	 * we only evaluate 2^a by using a polynomial */
	x = _mm256_or_si256(x, _mm256_castpd_si256(V_INT64_DOUBLE));
	vx = _mm256_sub_pd(_mm256_castsi256_pd(x), V_INT64_DOUBLE);
	vx = _mm256_mul_pd(vx, V_K_2_INV_S3);
	
	vx_1 = _mm256_floor_pd(vx);
	vx_2 = _mm256_add_pd(vx_1, V_DOUBLE_INT64);
	vt = _mm256_sub_epi64(_mm256_castpd_si256(vx_2), _mm256_castpd_si256(V_DOUBLE_INT64));	
	vt = _mm256_slli_epi64(vt, 52);
	
	/* evaluate 2^a */
	vx_2 = _mm256_sub_pd(vx, vx_1);
	vsum = _mm256_fmadd_pd(_mm256_castsi256_pd(EXP_COFF[0]), vx_2, _mm256_castsi256_pd(EXP_COFF[1]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[2]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[3]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[4]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[5]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[6]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[7]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[8]));
	vsum = _mm256_fmadd_pd(vsum, vx_2, _mm256_castsi256_pd(EXP_COFF[9]));
	
	/* combine to compute 2^x */
	vres = _mm256_add_epi64(vt, _mm256_castpd_si256(vsum));

	/* compute the Bernoulli value */
	vres_mantissa = _mm256_and_si256(vres, V_EXP_MANTISSA_MASK);
	vres_mantissa = _mm256_or_si256(vres_mantissa, V_RES_MANTISSA);
	
	vres_exponent = _mm256_srli_epi64(vres, EXP_MANTISSA_PRECISION);
	vres_exponent = _mm256_add_epi64(vres_exponent, V_RES_EXPONENT);
	vres_exponent = _mm256_sllv_epi64(V_1, vres_exponent);
	
	vr_mantissa = _mm256_loadu_si256((__m256i *)r);
	vr_exponent = _mm256_srli_epi64(vr_mantissa, R_MANTISSA_PRECISION);
	vr_mantissa = _mm256_and_si256(vr_mantissa, V_R_MANTISSA_MASK);
	vr_exponent2 = _mm256_set_epi64x(r[35], r[34], r[33], r[32]);
	vr_exponent2 = _mm256_slli_epi64(vr_exponent2, 64 - R_MANTISSA_PRECISION);
	vr_exponent = _mm256_or_si256(vr_exponent, vr_exponent2);

	/* (res == 1.0) || ((r_mantissa < res_mantissa) && (r_exponent < (1 << res_exponent))) */
	vres_eq_1 = _mm256_cmpeq_epi64(vres, V_DOUBLE_ONE);
	vr_lt_vres_mantissa = _mm256_sub_epi64(vr_mantissa, vres_mantissa);	
	vr_lt_vres_exponent = _mm256_sub_epi64(vr_exponent, vres_exponent);
	
	k = _mm256_and_si256(vr_lt_vres_mantissa, vr_lt_vres_exponent);
	k = _mm256_or_si256(k, vres_eq_1);

	_mm256_store_si256((__m256i *)(b), k);
}


/* make sure that Pr(rerun the PRG)<=2^(-64) */
static inline void uniform_sampler_S1(unsigned char *r, __m256i *y1, __m256i *y2)
{
	uint64_t sample[8] __attribute__ ((aligned (32)));
	uint32_t i = 0, j = 0;
	uint64_t x;
	
	while (j < 8)
	{
		do
		{	/* we ignore the low probability of rerunning the PRG */
			x = *((uint32_t *)(r + UNIFORM_SIZE * (i++)));
		} while (1 ^ ((x - UNIFORM_Q_S1) >> 63));

		x = x - ((((x * BARRETT_FACTOR_S1) >> BARRETT_BITSHIFT) + 1) * BINARY_SAMPLER_K_S1);
		x = x + (x >> 63) * BINARY_SAMPLER_K_S1;
		
		sample[j++] = x;
		//printf("%lu\t", x);
	}
	
	*y1 = _mm256_load_si256((__m256i *)(sample));
	*y2 = _mm256_load_si256((__m256i *)(sample + 4));
}

static inline void uniform_sampler_S2(unsigned char *r, __m256i *y1, __m256i *y2)
{
	uint64_t sample[8] __attribute__ ((aligned (32)));
	uint32_t i = 0, j = 0;
	uint64_t x;
	
	while (j < 8)
	{
		do
		{	/* we ignore the low probability of rerunning the PRG */
			x = *((uint32_t *)(r + UNIFORM_SIZE * (i++)));
		} while (1 ^ ((x - UNIFORM_Q_S2) >> 63));

		x = x - ((((x * BARRETT_FACTOR_S2) >> BARRETT_BITSHIFT) + 1) * BINARY_SAMPLER_K_S2);
		x = x + (x >> 63) * BINARY_SAMPLER_K_S2;
		
		sample[j++] = x;
		//printf("%lu\t", x);
	}
	
	*y1 = _mm256_load_si256((__m256i *)(sample));
	*y2 = _mm256_load_si256((__m256i *)(sample + 4));
}

static inline void uniform_sampler_S3(unsigned char *r, __m256i *y1, __m256i *y2)
{
	uint64_t sample[8] __attribute__ ((aligned (32)));
	uint32_t i = 0, j = 0;
	uint64_t x;
	
	while (j < 8)
	{
		do
		{	/* we ignore the low probability of rerunning the PRG */
			x = *((uint32_t *)(r + UNIFORM_SIZE * (i++)));
		} while (1 ^ ((x - UNIFORM_Q_S3) >> 63));

		x = x - ((((x * BARRETT_FACTOR_S3) >> BARRETT_BITSHIFT) + 1) * BINARY_SAMPLER_K_S3);
		x = x + (x >> 63) * BINARY_SAMPLER_K_S3;
		
		sample[j++] = x;
		//printf("%lu\t", x);
	}
	
	*y1 = _mm256_load_si256((__m256i *)(sample));
	*y2 = _mm256_load_si256((__m256i *)(sample + 4));
}


//void gaussian_sampler_S1(unsigned char *seed, int64_t sample[MIFE_NMODULI][MIFE_N], uint32_t slen)
void gaussian_sampler_S1(aes256ctr_ctx *state, uint32_t sample[MIFE_NMODULI][MIFE_N], uint32_t slen)
{
	__m256i v_x, v_y1, v_y2, v_z, v_b_in;
	uint64_t z[8] __attribute__ ((aligned (32)));
	uint64_t b[8] __attribute__ ((aligned (32)));
	
	//unsigned char r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE + 1] __attribute__ ((aligned (32)));
	unsigned char r[384];
	unsigned char *r1;
	
	uint32_t i = 8, j = 0;
	uint64_t k;
	//aes256ctr_ctx state;
	//aes256ctr_init(&state, seed, 0);
	
	const uint32_t AES_ROUNDS=3;

	uint64_t mod;
	uint32_t mod1, mod2, mod3;

	while (j < slen)
	{
		do
		{
			if (i == 8)
			{
				/* x<--D_sigma_0, y<--U([0,k-1]), z=kx+y */
				//fastrandombytes(r, 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE + 1);
				//aes256ctr_squeezeblocks(r, AES_ROUNDS, &state);
				aes256ctr_squeezeblocks(r, AES_ROUNDS, state);
				
				uniform_sampler_S1(r + 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE), &v_y1, &v_y2);
				
				r1 = r;
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, V_K_K_K_K_S1);
				v_z = _mm256_add_epi64(v_x, v_y1);
				_mm256_store_si256((__m256i *)(z), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y1);
				bernoulli_sampler_S1(b, v_b_in, r1 + BASE_TABLE_SIZE);
				
				r1 = r + (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE);
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, V_K_K_K_K_S1);
				v_z = _mm256_add_epi64(v_x, v_y2);
				_mm256_store_si256((__m256i *)(z + 4), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y2);
				bernoulli_sampler_S1(b + 4, v_b_in, r1 + BASE_TABLE_SIZE);

				i = 0;
			}
			
			k = (r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE] >> i) & 0x1;
			i++;			
		} while (1 ^ ((b[i - 1] & ((z[i - 1] | -z[i - 1]) | (k | -k))) >> 63)); /* rejection condition: b=0 or ((b=1) && (z=0) && (k=0)) */

		mod=z[i-1];

		/*
		if(mod>10*225){
			printf("\n----LARGE----\n");
		}
		*/

		mod1=mod_prime(mod, 0);
		mod2=mod_prime(mod, 1);
		mod3=mod_prime(mod, 2);		
		
		sample[0][j]=(1-k)*mod1+k*mod_prime(MIFE_MOD_Q_I[0]-mod1, 0);
		sample[1][j]=(1-k)*mod2+k*mod_prime(MIFE_MOD_Q_I[1]-mod2, 1);
		sample[2][j]=(1-k)*mod3+k*mod_prime(MIFE_MOD_Q_I[2]-mod3, 2);
		j++;

	}

}

//void gaussian_sampler_S2(unsigned char *seed, uint32_t sample[MIFE_NMODULI][MIFE_N], uint32_t slen)
void gaussian_sampler_S2(aes256ctr_ctx *state, uint32_t sample[MIFE_NMODULI][MIFE_N], uint32_t slen)
{
	__m256i v_x, v_y1, v_y2, v_z, v_b_in;
	uint64_t z[8] __attribute__ ((aligned (32)));
	uint64_t b[8] __attribute__ ((aligned (32)));
	
	//unsigned char r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE + 1] __attribute__ ((aligned (32)));
	unsigned char r[384];
	unsigned char *r1;
	
	uint32_t i = 8, j = 0;
	uint64_t k;
	//aes256ctr_ctx state;
	//aes256ctr_init(&state, seed, 0);
	
	const uint32_t AES_ROUNDS=3;

	uint64_t mod;
	uint32_t mod1, mod2, mod3;

	while (j < slen)
	{
		do
		{
			if (i == 8)
			{
				/* x<--D_sigma_0, y<--U([0,k-1]), z=kx+y */
				//fastrandombytes(r, 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE + 1);
				//aes256ctr_squeezeblocks(r, AES_ROUNDS, &state);
				aes256ctr_squeezeblocks(r, AES_ROUNDS, state);
				
				uniform_sampler_S2(r + 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE), &v_y1, &v_y2);
				
				r1 = r;
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, V_K_K_K_K_S2);
				v_z = _mm256_add_epi64(v_x, v_y1);
				_mm256_store_si256((__m256i *)(z), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y1);
				bernoulli_sampler_S2(b, v_b_in, r1 + BASE_TABLE_SIZE);
				
				r1 = r + (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE);
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, V_K_K_K_K_S2);
				v_z = _mm256_add_epi64(v_x, v_y2);
				_mm256_store_si256((__m256i *)(z + 4), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y2);
				bernoulli_sampler_S2(b + 4, v_b_in, r1 + BASE_TABLE_SIZE);

				i = 0;
			}
			
			k = (r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE] >> i) & 0x1;
			i++;			
		} while (1 ^ ((b[i - 1] & ((z[i - 1] | -z[i - 1]) | (k | -k))) >> 63)); /* rejection condition: b=0 or ((b=1) && (z=0) && (k=0)) */
		
		mod=z[i-1];
		/*
		if(mod>10*258376412UL){
			printf("\n----LARGE----\n");
		}
		*/

		mod1=mod_prime(mod, 0);
		mod2=mod_prime(mod, 1);
		mod3=mod_prime(mod, 2);		

		sample[0][j]=(1-k)*mod1+k*mod_prime(MIFE_MOD_Q_I[0]-mod1, 0);
		sample[1][j]=(1-k)*mod2+k*mod_prime(MIFE_MOD_Q_I[1]-mod2, 1);
		sample[2][j]=(1-k)*mod3+k*mod_prime(MIFE_MOD_Q_I[2]-mod3, 2);
		j++;

	}
}

//void gaussian_sampler_S3(unsigned char *seed, uint32_t sample[MIFE_NMODULI][MIFE_N], uint32_t slen)
void gaussian_sampler_S3(aes256ctr_ctx *state, uint32_t sample[MIFE_NMODULI][MIFE_N], uint32_t slen)
{
	__m256i v_x, v_y1, v_y2, v_z, v_b_in;
	uint64_t z[8] __attribute__ ((aligned (32)));
	uint64_t b[8] __attribute__ ((aligned (32)));
	
	//unsigned char r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE + 1] __attribute__ ((aligned (32)));
	unsigned char r[384];
	unsigned char *r1;
	
	uint32_t i = 8, j = 0;
	uint64_t k;
	//aes256ctr_ctx state;
	//aes256ctr_init(&state, seed, 0);
	
	const uint32_t AES_ROUNDS=3;

	uint64_t mod;
	uint32_t mod1, mod2, mod3;

	while (j < slen)
	{
		do
		{
			if (i == 8)
			{
				/* x<--D_sigma_0, y<--U([0,k-1]), z=kx+y */
				//fastrandombytes(r, 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE + 1);
				//aes256ctr_squeezeblocks(r, AES_ROUNDS, &state);
				aes256ctr_squeezeblocks(r, AES_ROUNDS, state);
				
				uniform_sampler_S3(r + 2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE), &v_y1, &v_y2);
				
				r1 = r;
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, V_K_K_K_K_S3);
				v_z = _mm256_add_epi64(v_x, v_y1);
				_mm256_store_si256((__m256i *)(z), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y1);
				bernoulli_sampler_S3(b, v_b_in, r1 + BASE_TABLE_SIZE);
				
				r1 = r + (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE);
				v_x = cdt_sampler(r1);
				v_x = _mm256_mul_epu32(v_x, V_K_K_K_K_S3);
				v_z = _mm256_add_epi64(v_x, v_y2);
				_mm256_store_si256((__m256i *)(z + 4), v_z);
				/* b<--Bernoulli(exp(-y(y+2kx)/2sigma_0^2)) */
				v_b_in = _mm256_add_epi64(v_z, v_x);
				v_b_in = _mm256_mul_epu32(v_b_in, v_y2);
				bernoulli_sampler_S3(b + 4, v_b_in, r1 + BASE_TABLE_SIZE);

				i = 0;
			}
			
			k = (r[2 * (BASE_TABLE_SIZE + BERNOULLI_TABLE_SIZE) + UNIFORM_REJ * UNIFORM_SIZE] >> i) & 0x1;
			i++;			
		} while (1 ^ ((b[i - 1] & ((z[i - 1] | -z[i - 1]) | (k | -k))) >> 63)); /* rejection condition: b=0 or ((b=1) && (z=0) && (k=0)) */
		
		mod=z[i-1];
		/*
		if(z[i-1]>10UL*516752822UL){
			printf("\n----LARGE----\n");
		}
		*/

		mod1=mod_prime(mod, 0);
		mod2=mod_prime(mod, 1);
		mod3=mod_prime(mod, 2);		

		sample[0][j]=(1-k)*mod1+k*mod_prime(MIFE_MOD_Q_I[0]-mod1, 0);
		sample[1][j]=(1-k)*mod2+k*mod_prime(MIFE_MOD_Q_I[1]-mod2, 1);
		sample[2][j]=(1-k)*mod3+k*mod_prime(MIFE_MOD_Q_I[2]-mod3, 2);
		j++;
	}
}

