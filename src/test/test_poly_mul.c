#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <gmp.h>

#include "../modred.h"

#include "../rng.h"
#include "../rlwe_mife.h"

#include "../crt.h"
#include "../poly_mul.h"
#include "../sample.h"

#define N_TESTS 10

long long cpucycles(void)
{
	unsigned long long result;
	asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
		: "=a" (result) ::  "%rdx");
	return result;
}

void
fprintBstr(char *S, unsigned char *A, unsigned long long L)
{
	unsigned long long  i;

	printf("%s", S);

	for ( i=0; i<L; i++ )
		printf("%02X", A[i]);

	if ( L == 0 )
		printf("00");

	printf("\n");
}

static
void
fprintPoly_gp(char * s, uint64_t * a, int n, uint64_t mod)
{
	int i;

	printf("%s", s);
	
	printf("{");
	
	//("Mod(%ld,%ld)*x^%d + ",a[i],mod,i)
	for (i = 0; i < n; ++i) {
		printf("Mod(%ld,%ld)*x^%d + ",a[i],mod,i);
	}
	
	printf("}");
	
	printf("\n");
}

static
void
fprintPoly_small_gp(char * s, uint32_t * a, int n, uint64_t mod)
{
	int i;

	printf("%s", s);
	
	printf("{");
	
	//("Mod(%ld,%ld)*x^%d + ",a[i],mod,i)
	for (i = 0; i < n; ++i) {
		printf("Mod(%ld,%ld)*x^%d + ",a[i],mod,i);
	}
	
	printf("}");
	
	printf("\n");
}

static
void
fprintPoly(char * s, uint64_t * a, int n)
{
	int i;

	printf("%s", s);
	
	printf("{");
	
	for (i = 0; i < n; ++i) {
		printf("%ld, ", a[i]);
	}
	
	printf("}");
	
	printf("\n");
}

static
void
fprintPoly_small(char * s, uint32_t * a, int n)
{
	int i;

	printf("%s", s);
	
	printf("{");
	
	for (i = 0; i < n; ++i) {
		printf("%ld, ", a[i]);
	}
	
	printf("}");
	
	printf("\n");
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

static
int64_t
int_mul_mod_gmp
(uint64_t a, uint64_t b, uint64_t m)
{
	//TODO
	return a*b % m;
}

void
poly_mul_mod_ref
(const uint64_t a[MIFE_N], const uint64_t b[MIFE_N], uint64_t c[MIFE_N], const uint64_t mod)
{
	int i, j;
	uint64_t res[2*MIFE_N-1];
	for (i = 0; i < MIFE_N; ++i) {
		c[i] = 0;
	}
	for (i = 0; i < MIFE_N; ++i) {
		for (j = 0; j < MIFE_N; ++j) {
			res[i+j] = (res[i+j] + ((uint64_t)int_mul_mod( a[i], b[j], mod )) ) % mod;
		}
	}
	for (i = MIFE_N; i < 2*MIFE_N-1; ++i) {
		c[i - MIFE_N] = (c[i - MIFE_N] + res[i - MIFE_N] - res[i]) % mod;
	}
	c[MIFE_N-1] = (c[MIFE_N-1] + res[MIFE_N-1]) % mod;
}

void
poly_mul_mod_gmp
(const uint64_t a[MIFE_N], const uint64_t b[MIFE_N], uint64_t c[MIFE_N], const uint64_t mod)
{
	int i, j;
	uint64_t res[2*MIFE_N-1];
	for (i = 0; i < MIFE_N; ++i) {
		c[i] = 0;
	}
	for (i = 0; i < MIFE_N; ++i) {
		for (j = 0; j < MIFE_N; ++j) {
			res[i+j] = (res[i+j] + ((uint64_t)int_mul_mod_gmp( a[i], b[j], mod )) ) % mod;
		}
	}
	for (i = MIFE_N; i < 2*MIFE_N-1; ++i) {
		c[i - MIFE_N] = (c[i - MIFE_N] + res[i - MIFE_N] - res[i]) % mod;
	}
	c[MIFE_N-1] = (c[MIFE_N-1] + res[MIFE_N-1]) % mod;
}


int test_rlwe_mife_poly_mul()
{
	// Declarate variables
	uint64_t a[MIFE_N], b[MIFE_N], c_ref[MIFE_N], c[MIFE_N], a_crt_inv[MIFE_N];
	uint32_t a_crt[MIFE_NMODULI][MIFE_N], b_crt[MIFE_NMODULI][MIFE_N], c_crt[MIFE_NMODULI][MIFE_N];
	uint32_t c_crt_ntt[MIFE_NMODULI][MIFE_N];
	uint64_t c_ntt[MIFE_N];

	unsigned char entropy_input[48];

	uint64_t i, j, k;
	uint64_t CLOCK1, CLOCK2;
	uint64_t CLOCK_gen, CLOCK_enc, CLOCK_ref, CLOCK_mul, CLOCK_dec;

	CLOCK1 = 0;
	CLOCK2 = 0;
	CLOCK_gen = CLOCK_ref = CLOCK_mul = CLOCK_enc = CLOCK_dec = 0;

	time_t t;
	// Intializes random number generator
	srand((unsigned) time(&t));

	for (i=0; i<48; i++){
		//entropy_input[i] = rand()%256;
		entropy_input[i] = i;
	}
	randombytes_init(entropy_input, NULL, 256);

	// Print parameters
	printf("\n\n\n");
	printf("MIFE_Q1=%d\n", MIFE_Q1);
	printf("MIFE_Q2=%d\n", MIFE_Q2);
	printf("MIFE_Q=%lu\n", MIFE_Q);
	printf("\n");

/*
	printf("\n### TEST GMP MULTIPLICATION ###\n");
	// initilze the state object for the random generator functions
	gmp_randstate_t rstate;
	// initialize state for a Mersenne Twister algorithm. This algorithm is fast and has good randomness properties.
	gmp_randinit_mt(rstate);
	// create the generator seed for the random engine to reference 
	gmp_randseed_ui(rstate, rand());
	mpz_t a1,a2, mul, mod;
	for(i = 0; i < N_TESTS; i++) {
		mpz_init2(a1,64);//initialize a 100 bit number
		mpz_init2(a2,64);//initialize a 100 bit number
		mpz_init2(mod,64);//initialize a 100 bit number
		mpz_init2(mul,128);//initialize a 100 bit number
		mpz_urandomb(a1, rstate, 64); //generate a 100 bit random number
		mpz_urandomb(a2, rstate, 64); //generate a 100 bit random number
		mpz_urandomb(mod, rstate, 64); //generate a 100 bit random number

		mpz_mul(mul, a1, a2);
		mpz_mod(mul, mul, mod);

		gmp_printf("%Zd * %Zd mod %Zd = %Zd\n",a1,a2,mod, mul);
		mpz_clear(a1);
		mpz_clear(a2);
		mpz_clear(mod);
		mpz_clear(mul);
	}

	printf("\n### TEST CRT MULTIPLICATION OF INTEGERS ###\n");

	for (i = 0; i < N_TESTS; ++i) {

		printf("Test #%ld - ", i);
	}

	printf("\nTest of CRT for integer arithmetic DONE!\n\n");
*/
#ifndef TEST_ADDITION
	printf("\n### TEST POLYNOMIAL ADDITION USING CRT ###\n");

	for (i = 0; i < N_TESTS; ++i) {

		printf("Test #%ld - ", i);

		// Random operands
		sample_uniform(a);
		sample_uniform(b);

		// CRT
		crt_convert(a, a_crt);
		crt_convert(b, b_crt);

		// Functional verification of convert-reverse
		crt_reverse(a_crt_inv, a_crt);
		for (j = 0; j < MIFE_N; ++j) {
			if (a[j] != a_crt_inv[j]) {
				printf("Error in the CRT conversion at index %ld\n", j);
				//break;
			}
		}

		// Addition without CRT
		//add_mod(a, b, c_ref, MIFE_Q);
		mpz_t c_i, mod;
		mpz_init(c_i);
		mpz_init(mod);
		mpz_set_ui(mod, MIFE_Q);
		for (j = 0; j < MIFE_N; ++j) {
			c_ref[j] = (a[j] + b[j]);
			mpz_set_ui(c_i, c_ref[j]);
			mpz_mod(c_i, c_i, mod);
			c_ref[j] = mpz_get_ui(c_i);
		}
		mpz_clear(c_i);
		mpz_clear(mod);

		// Addition with CRT
		for (j = 0; j < MIFE_NMODULI; j++) {
			add_mod(a_crt[j], b_crt[j], c_crt[j], MIFE_MOD_Q_I[j]);
		}
		crt_reverse(c, c_crt);
/*
		fprintPoly("\n\na = ",a,1);
		fprintPoly("b = ",b,1);
		printf("MIFE_Q1 = %lu\n", MIFE_Q1);
		fprintPoly("a_crt = ",a_crt,1);
		fprintPoly("b_crt = ",b_crt,1);
		fprintPoly("c_crt = ",c_crt,1);
		printf("MIFE_Q2 = %lu\n", MIFE_Q2);
		fprintPoly("a_crt = ",a_crt[1],1);
		fprintPoly("b_crt = ",b_crt[1],1);
		fprintPoly("c_crt = ",c_crt[1],1);
		printf("MIFE_Q = %lu\n", MIFE_Q);
		fprintPoly("c_ref = ",c_ref,1);
		fprintPoly("c = ",c,1);
		printf("\n\n");
*/
		// Functional verification:
		for (j = 0; j < MIFE_N; ++j) {
			if (c[j] != c_ref[j]) {
				printf("Error at iteration %ld\n", j);
				break;
			}
		}
		printf("\n");
	}

	printf("\nTest of polynomial addition DONE!\n\n");
#endif


	printf("\n### TEST POLYNOMIAL MULTIPLICATION ###\n");

	for(i = 0; i < N_TESTS; i++) {
		printf("i : %lu\n",i);
		for (j = 0; j < MIFE_N; ++j) {
			c_ref[j] = 0;
			c[j] = 0;
			c_crt[0][j] = 0;
			c_crt[1][j] = 0;
		}

		//Generation of the operand polynomials a and b
		CLOCK1=cpucycles();
		sample_uniform(a);
		sample_uniform(b);
		//sample_ones(a);
		//sample_ones(b);
		CLOCK2=cpucycles();
		CLOCK_gen += (CLOCK2-CLOCK1);

		//fprintPoly("a = ", a, 10);
		//fprintPoly("b = ", b, 10);
//		fprintPoly_gp("a = ", a, 4, MIFE_Q);
//		fprintPoly_gp("b = ", b, 4, MIFE_Q);

		//Conversion of a anb d to CRT domain
		CLOCK1=cpucycles();
		crt_convert(a, a_crt);
		crt_convert(b, b_crt);
		CLOCK2=cpucycles();
		CLOCK_enc += (CLOCK2-CLOCK1);

		//fprintPoly_small("a1 = ", a_crt[0], 10);
		//fprintPoly_small("a2 = ", a_crt[1], 10);
		//fprintPoly_small("b1 = ", b_crt[0], 10);
		//fprintPoly_small("b2 = ", b_crt[1], 10);
//		fprintPoly_small_gp("a1 = ", a_crt[0], 4, MIFE_Q1);
//		fprintPoly_small_gp("a2 = ", a_crt[1], 4, MIFE_Q2);
//		fprintPoly_small_gp("b1 = ", b_crt[0], 4, MIFE_Q1);
//		fprintPoly_small_gp("b2 = ", b_crt[1], 4, MIFE_Q2);

		// Functional verification of convert-reverse
		crt_reverse(a_crt_inv, a_crt);
		for (j = 0; j < MIFE_N; ++j) {
			if (a[j] != a_crt_inv[j]) {
				printf("Error in the CRT conversion at index %ld\n", j);
				//break;
			}
		}

		//Reference multiplication
		CLOCK1=cpucycles();
		//poly_mul_mod_ref(a, b, c_ref, MIFE_Q);
		uint64_t res[2*MIFE_N-1] = {0};
		uint64_t index, mac;
		mpz_t c_index, q, aa, bb;
		mpz_init(c_index);
		mpz_init(q);
		mpz_init(aa);
		mpz_init(bb);
		mpz_set_ui(q, MIFE_Q);
/*		for (j = 0; j < MIFE_N; ++j) {
			for (k = 0; k < MIFE_N; ++k) {
				if (j+k >= MIFE_N) {
					index = j + k - MIFE_N;
				} else {
					index = j + k;
				}
				//printf("a[j] = %ld\n", a[j]);
				//printf("b[k] = %ld\n", b[k]);
				//printf("acc = %ld\n", c_ref[index]);
				mpz_set_ui(aa, a[j]);
				mpz_set_ui(bb, b[k]);
				mpz_set_ui(c_index, c_ref[index]);
//				gmp_printf("%Zd * %Zd + %Zd = \n",aa, bb, c_index);
				mpz_mul(aa, aa, bb);
				mpz_add(aa, aa, c_index);
				mpz_mod(aa, aa, q);
				c_ref[index] = mpz_get_ui(aa);
				//printf("c[index] = %ld\n\n", c_ref[index]);
			}
		}*/
		//Multiplication
		for (j = 0; j < MIFE_N; ++j) {
			for (k = 0; k < MIFE_N; ++k) {
				mpz_set_ui(aa, a[j]);
				mpz_set_ui(bb, b[k]);
				mpz_set_ui(c_index, res[j+k]);
				mpz_mul(aa, aa, bb);
				mpz_add(aa, aa, c_index);
				mpz_mod(aa, aa, q);
				res[j+k] = mpz_get_ui(aa);
			}
		}
		//fprintPoly("res = ", res, 10);
		//Negative wrap
		for (j = MIFE_N; j < 2*MIFE_N-1; ++j) {
			mpz_set_ui(aa, res[j - MIFE_N]);
			mpz_set_ui(bb, res[j]);
			mpz_sub(aa, aa, bb);
			mpz_set_ui(bb, c_ref[j - MIFE_N]);
			mpz_add(aa, aa, bb);
			mpz_mod(aa, aa, q);
			c_ref[j - MIFE_N] = mpz_get_ui(aa);
		}
		mac = (c_ref[MIFE_N-1] + res[MIFE_N-1]);
		mpz_set_ui(aa, c_ref[MIFE_N-1]);
		mpz_set_ui(bb, res[MIFE_N-1]);
		mpz_add(aa, aa, bb);
		mpz_mod(aa, aa, q);
		c_ref[MIFE_N-1] =  mpz_get_ui(aa);
		mpz_clear(c_index);
		mpz_clear(q);
		mpz_clear(aa);
		mpz_clear(bb);
		CLOCK2=cpucycles();
		CLOCK_ref += (CLOCK2-CLOCK1);
		//fprintPoly("c_ref = ", c_ref, 10);


		//CRT domain multiplication
		CLOCK1=cpucycles();
		for (j = 0; j < MIFE_NMODULI; j++) {
			//fprintPoly_small("a = ", a_crt[j],10);
			//fprintPoly_small("b = ", b_crt[j],10);
			//fprintPoly_small("c_pre = ", c_crt[j],10);
			poly_mul_mod(a_crt[j], b_crt[j], c_crt[j], MIFE_MOD_Q_I[j]);
			poly_mul_mod_ntt(a_crt[j], b_crt[j], c_crt_ntt[j], MIFE_MOD_Q_I[j]);
			//fprintPoly_small("c_post = ", c_crt[j],10);
//			fprintPoly_small_gp("c_crt = ", c_crt[j], 4, MIFE_MOD_Q_I[j]);
		}
		CLOCK2=cpucycles();
		CLOCK_mul += (CLOCK2-CLOCK1);

		//Conversion of a and b from CRT domain
		//fprintPoly_small("c1 = ", c_crt[0],10);
		//fprintPoly_small("c2 = ", c_crt[1],10);
		CLOCK1=cpucycles();
		crt_reverse(c, c_crt);
		crt_reverse(c_ntt, c_crt_ntt);
		CLOCK2=cpucycles();
		CLOCK_dec += (CLOCK2-CLOCK1);
		//fprintPoly("c = ", c, 10);
//		fprintPoly_gp("c = ", c, 4, MIFE_Q);
//		fprintPoly_gp("c_ref = ", c_ref, 4, MIFE_Q);

		// Functional verification:
		for (j = 0; j < MIFE_N; ++j) {
			if (c[j] != c_ref[j]) {
				printf("c_ref = %ld --- c = %ld\n", c_ref[j], c[j]);
				printf("Error at iteration %ld\n", j);
				break;
			}
			if (c_ntt[j] != c_ref[j]) {
				//printf("c_ref = %ld --- c_ntt = %ld\n", c_ref[j], c_ntt[j]);
				//printf("Error at iteration %ld\n", j);
				break;
			}
		}
		printf("\n");
	}

	printf("Repeat is : %d\n",N_TESTS);
	printf("Average times generation \t %lu \n", CLOCK_gen/N_TESTS);
	printf("Average times crt_convert \t %lu \n",CLOCK_enc/N_TESTS);
	printf("Average times crt_reverse \t %lu \n",CLOCK_dec/N_TESTS);
	printf("Average times ref: \t %lu \n",CLOCK_ref/N_TESTS);
	printf("Average times mul: \t %lu \n",CLOCK_mul/N_TESTS);

	printf("\nTest of polynomial multiplication DONE!\n\n");

	return 0;
}

int main()
{
	test_rlwe_mife_poly_mul();
	return 0;
}
