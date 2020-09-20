#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "../rng.h"
#include "../rlwe_mife.h"

#define N_TESTS 3

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


int test_rlwe_mife()
{
	// Declarate variables
	uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N];
	uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N];
	uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N];
	uint32_t sk_y[MIFE_NMODULI][MIFE_N];
	uint32_t d_y[MIFE_NMODULI][MIFE_N];

	uint32_t m[MIFE_L];
	uint32_t y[MIFE_L];

	uint32_t m_crt[MIFE_NMODULI][MIFE_L];
	uint32_t y_crt[MIFE_NMODULI][MIFE_L];

	uint64_t dy[MIFE_N];

	unsigned char entropy_input[48];

	uint64_t i, j, k;
	uint64_t CLOCK1, CLOCK2;
	uint64_t CLOCK_su, CLOCK_enc, CLOCK_kp, CLOCK_dec;

	CLOCK1 = 0;
	CLOCK2 = 0;
	CLOCK_su = CLOCK_kp = CLOCK_enc = CLOCK_dec = 0;

	time_t t;
	// Intializes random number generator
	srand((unsigned) time(&t));

	for (i=0; i<48; i++){
		//entropy_input[i] = rand()%256;
		entropy_input[i] = i;
	}
	randombytes_init(entropy_input, NULL, 256);

	// Print parameters
	printf("MIFE_Q1=%d\n", MIFE_Q1);
	printf("MIFE_Q2=%d\n", MIFE_Q2);
	printf("MIFE_Q=%lu\n", MIFE_Q);
	printf("\n");

	// Sample message and y
	//m[0] = MIFE_B-1;
	//y[1] = MIFE_B-1;
	for (i = 0; i < MIFE_L; ++i) {
		m[i] = 1;
		y[i] = 1;
	}

	for(i = 0; i < N_TESTS; i++) {
		printf("i : %lu\n",i);

		//Generation of master secret key sk and master public key pk pair
		CLOCK1=cpucycles();
		rlwe_mife_setup(mpk, msk);
		CLOCK2=cpucycles();	
		CLOCK_su += (CLOCK2-CLOCK1);
		printf("Keysetup done \n");
	
		//Encryption of the message m
		CLOCK1=cpucycles();
		rlwe_mife_encrypt(m, mpk, c);
		CLOCK2=cpucycles();	
		CLOCK_enc += (CLOCK2-CLOCK1);
		printf("Encryption done \n");

		//Generation of the key for decrypting m·y
		CLOCK1=cpucycles();
		rlwe_mife_keygen(y, msk, sk_y);
		CLOCK2=cpucycles();	
		CLOCK_kp += (CLOCK2-CLOCK1);
		printf("Keygen done \n");

		//Decryption of m·y
		CLOCK1=cpucycles();
		rlwe_mife_decrypt(c, y, sk_y, d_y);
		CLOCK2=cpucycles();	
		CLOCK_dec += (CLOCK2-CLOCK1);
		printf("Decrypt done \n");

		// Functional verification:
		crt_reverse(dy, d_y);

		uint64_t xy;
		xy = 0;
		for (j = 0; j < MIFE_L; ++j) {
			xy += (uint64_t)m[j]*y[j];
		}
		//xy = xy*MIFE_SCALE_M;

		printf("xy = %ld and dy = %ld\n", xy, (MIFE_P*dy[0]/MIFE_Q) );

		// TODO: Double check the result in CRT domain
		crt_convert_generic(m, m_crt, MIFE_L);
		crt_convert_generic(y, y_crt, MIFE_L);
		uint64_t mxm;
		mxm = 0;
		for (k = 0; k < MIFE_L; ++k) {
			for (j = 0; j < MIFE_NMODULI; ++j) {
				//mxm = (uint64_t)m_crt[j][k] * y_crt;
				//mxm = (uint64_t)m_crt[j][k] * MIFE_SCALE_M_MOD_Q_I[j];
				//m_crt[j][k] = mod_red(mxm, MIFE_MOD_Q_I[j]);
			}
		}
		printf("TEST %lu DONE!\n", i);
	}

	printf("Repeat is : %ld\n",N_TESTS);
	printf("Average times setup \t \t %lu \n", CLOCK_su/N_TESTS);
	printf("Average times enc: \t \t %lu \n",CLOCK_enc/N_TESTS);
	printf("Average times key_pair: \t %lu \n",CLOCK_kp/N_TESTS);
	printf("Average times dec: \t \t %lu \n",CLOCK_dec/N_TESTS);

	return 0;
}

int main()
{

	test_rlwe_mife();
	return 0;
}
