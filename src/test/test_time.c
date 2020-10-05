#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gmp.h>

#include "../sample.h"
#include "../rlwe_mife.h"

#define N_TESTS 5

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


int test_rlwe_mife_vec_vec()			/*Only vector-vector multiplication*/
{
	// Declarate variables
	uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N];
	uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N];
	uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N];
	uint32_t sk_y[MIFE_NMODULI][MIFE_N];
	//uint32_t d_y[MIFE_NMODULI][MIFE_N];

	uint32_t m[MIFE_L];
	uint32_t y[MIFE_L];

	//uint32_t m_crt[MIFE_NMODULI][MIFE_L];
	//uint32_t y_crt[MIFE_NMODULI][MIFE_L];

	//uint64_t dy[MIFE_N];
	mpz_t dy[MIFE_N];
	mpz_t noise_gmp, Q_gmp, scale_M_gmp;


	//unsigned char entropy_input[48];

	uint64_t i, j, k;
	uint64_t CLOCK1, CLOCK2;
	uint64_t CLOCK_su, CLOCK_enc, CLOCK_kp, CLOCK_dec, CLOCK_extract;

	for(i=0;i<MIFE_N;i++){
		mpz_init(dy[i]);
	}
	mpz_init(noise_gmp);
	mpz_init(Q_gmp);
	mpz_init(scale_M_gmp);

	CLOCK1 = 0;
	CLOCK2 = 0;
	CLOCK_su = CLOCK_kp = CLOCK_enc = CLOCK_dec = 0;
	CLOCK_extract = 0;

	time_t t;
	// Intializes random number generator
	srand((unsigned) time(&t));

	/*
	for (i=0; i<48; i++){
		//entropy_input[i] = rand()%256;
		entropy_input[i] = i;
	}
	*/
	//randombytes_init(entropy_input, NULL, 256);
	
	//const char Q_string[] ="259809622039819";

	if(mpz_set_str(Q_gmp, MIFE_Q_str, 10)!=0){

		printf("--ERROR unable to set Q to gmp--\n");
		return 0;
	}
	//mpz_set_ui(noise_gmp, MIFE_P);
	if(mpz_set_str(noise_gmp, MIFE_P_str, 10)!=0){

		printf("--ERROR unable to set P to gmp--\n");
		return 0;
	}
	if(mpz_set_str(scale_M_gmp, MIFE_SCALE_M_str, 10)!=0){

		printf("--ERROR unable to set scaling factor M to gmp--\n");
		return 0;
	}

	mpz_mul_ui(noise_gmp, noise_gmp, 2);
	mpz_fdiv_q(noise_gmp, Q_gmp, noise_gmp);
	// Print parameters
	//printf("MIFE_Q1=%d\n", MIFE_Q1);
	//printf("MIFE_Q2=%d\n", MIFE_Q2);
	//printf("MIFE_Q=%lu\n", MIFE_Q);
	//printf("Noise tolerance=%lu\n", MIFE_Q/MIFE_P);

	for(i=0;i<MIFE_NMODULI;i++){
		printf("Q[i] : %u\n", MIFE_MOD_Q_I[i]);	
	}	

	gmp_printf("MIFE_Q=%Zd\n", Q_gmp);
	gmp_printf("Noise tolerance=%Zd\n", noise_gmp);
	printf("\n");

	for(i = 0; i < N_TESTS; i++) {
		printf("i : %lu\n",i);

		// Sample message and y
		sample_x(m);
		sample_y(y);

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
		rlwe_mife_decrypt_gmp(c, y, sk_y, dy);
		CLOCK2=cpucycles();
		CLOCK_dec += (CLOCK2-CLOCK1);
		printf("Decrypt done \n");

		// Functional verification
		k = 0;
		for (j = 0; j < MIFE_L; ++j) {
			k += (uint64_t)m[j]*y[j];
		}

		mpz_set_ui(noise_gmp, k);
		//mpz_mul_ui(noise_gmp, noise_gmp, MIFE_SCALE_M);
		mpz_mul(noise_gmp, noise_gmp, scale_M_gmp);
		mpz_sub(noise_gmp, dy[0], noise_gmp);

		//Extraction of the result (cancel scaling)
		CLOCK1=cpucycles();
		//round_extract(dy);
		round_extract_gmp(dy);		
		CLOCK2=cpucycles();
		CLOCK_extract += (CLOCK2-CLOCK1);
		printf("Extraction done \n");


		gmp_printf("xy = %ld and dy = %Zd\n", k, dy[0]);
		
		gmp_printf("dy: %Zd\n", dy[0]);
		gmp_printf("Noise is : %Zd\n",noise_gmp);

		if(mpz_cmp_ui(dy[0], k)!=0){
			printf("--ERROR---\n");
			break;
		}

		printf("TEST %lu DONE!\n\n", i);
	}

	printf("Repeat is : %d\n",N_TESTS);
	printf("Average times setup: \t \t %lu \n", CLOCK_su/N_TESTS);
	printf("Average times enc: \t \t %lu \n",CLOCK_enc/N_TESTS);
	printf("Average times key_pair: \t %lu \n",CLOCK_kp/N_TESTS);
	printf("Average times dec: \t \t %lu \n",CLOCK_dec/N_TESTS);
	printf("Average times extract: \t \t %lu \n",CLOCK_extract/N_TESTS);

	for(i=0;i<MIFE_N;i++){
		mpz_clear(dy[i]);
	}
	mpz_clear(noise_gmp);
	mpz_clear(Q_gmp);
	mpz_clear(scale_M_gmp);


	return 0;
}

int test_rlwe_mife_mat_vec()			/*Only matrix-vector multiplication*/
{
	// Declarate variables
	uint32_t mpk[MIFE_L+1][MIFE_NMODULI][MIFE_N];
	uint32_t msk[MIFE_L][MIFE_NMODULI][MIFE_N];
	uint32_t c[MIFE_L+1][MIFE_NMODULI][MIFE_N];
	uint32_t sk_y[MIFE_NMODULI][MIFE_N];
	//uint32_t d_y[MIFE_NMODULI][MIFE_N];

	uint32_t m[MIFE_L][MIFE_N];
	uint32_t y[MIFE_L];
	uint32_t res_vec[MIFE_N];

	//uint32_t m_crt[MIFE_NMODULI][MIFE_L];
	//uint32_t y_crt[MIFE_NMODULI][MIFE_L];

	//uint64_t dy[MIFE_N];
	mpz_t dy[MIFE_N];
	mpz_t noise_gmp, Q_gmp, scale_M_gmp;


	//unsigned char entropy_input[48];

	uint64_t i, j, k;
	uint64_t CLOCK1, CLOCK2;
	uint64_t CLOCK_su, CLOCK_enc, CLOCK_kp, CLOCK_dec, CLOCK_extract;

	for(i=0;i<MIFE_N;i++){
		mpz_init(dy[i]);
	}
	mpz_init(noise_gmp);
	mpz_init(Q_gmp);
	mpz_init(scale_M_gmp);

	CLOCK1 = 0;
	CLOCK2 = 0;
	CLOCK_su = CLOCK_kp = CLOCK_enc = CLOCK_dec = 0;
	CLOCK_extract = 0;

	time_t t;
	// Intializes random number generator
	srand((unsigned) time(&t));

	if(mpz_set_str(Q_gmp, MIFE_Q_str, 10)!=0){

		printf("--ERROR unable to set Q to gmp--\n");
		return 0;
	}
	if(mpz_set_str(noise_gmp, MIFE_P_str, 10)!=0){

		printf("--ERROR unable to set P to gmp--\n");
		return 0;
	}
	if(mpz_set_str(scale_M_gmp, MIFE_SCALE_M_str, 10)!=0){

		printf("--ERROR unable to set scaling factor M to gmp--\n");
		return 0;
	}

	mpz_mul_ui(noise_gmp, noise_gmp, 2);
	mpz_fdiv_q(noise_gmp, Q_gmp, noise_gmp);

	for(i=0;i<MIFE_NMODULI;i++){
		printf("Q[i] : %u\n", MIFE_MOD_Q_I[i]);	
	}	

	gmp_printf("MIFE_Q=%Zd\n", Q_gmp);
	gmp_printf("Noise tolerance=%Zd\n", noise_gmp);
	printf("\n");

	for(i = 0; i < N_TESTS; i++) {
		printf("i : %lu\n",i);

		// Sample message and y
		sample_y(y);
		for(k=0;k<MIFE_L;k++){
			sample_m(m[k]);
		}

		//Generation of master secret key sk and master public key pk pair
		CLOCK1=cpucycles();
		rlwe_mife_setup(mpk, msk);
		CLOCK2=cpucycles();	
		CLOCK_su += (CLOCK2-CLOCK1);
		printf("Keysetup done \n");
	
		//Encryption of the message m
		CLOCK1=cpucycles();
		rlwe_mife_encrypt_vec(m, mpk, c);
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
		rlwe_mife_decrypt_gmp_vec(c, y, sk_y, dy);
		CLOCK2=cpucycles();
		CLOCK_dec += (CLOCK2-CLOCK1);
		printf("Decrypt done \n");

		// Functional verification
		for(k=0;k<MIFE_N;k++){
			res_vec[k]=0;	
			for (j = 0; j < MIFE_L; ++j) {
				res_vec[k]= res_vec[k]+m[j][k]*y[j];
			}
		}

		mpz_set_ui(noise_gmp, res_vec[0]);
		//mpz_mul_ui(noise_gmp, noise_gmp, MIFE_SCALE_M);
		mpz_mul(noise_gmp, noise_gmp, scale_M_gmp);
		mpz_sub(noise_gmp, dy[0], noise_gmp);

		//Extraction of the result (cancel scaling)
		CLOCK1=cpucycles();
		//round_extract(dy);
		round_extract_gmp(dy);		
		CLOCK2=cpucycles();
		CLOCK_extract += (CLOCK2-CLOCK1);
		printf("Extraction done \n");

		/*
		for(k=0;k<MIFE_N;k++){
			gmp_printf("xy[%lu] = %ld and dy[%lu] = %Zd\n", k, res_vec[k], k, dy[k]);
		}
		*/
		
		//gmp_printf("dy: %Zd\n", dy[0]);
		//gmp_printf("Noise is : %Zd\n",noise_gmp);

		for(k=0;k<MIFE_N;k++){
			if(mpz_cmp_ui(dy[k], res_vec[k])!=0){
				printf("--ERROR---\n");
				return 0;
			}
		}
		
		printf("TEST %lu DONE!\n\n", i);
	}

	printf("Repeat is : %d\n",N_TESTS);
	printf("Average times setup: \t \t %lu \n", CLOCK_su/N_TESTS);
	printf("Average times enc: \t \t %lu \n",CLOCK_enc/N_TESTS);
	printf("Average times key_pair: \t %lu \n",CLOCK_kp/N_TESTS);
	printf("Average times dec: \t \t %lu \n",CLOCK_dec/N_TESTS);
	printf("Average times extract: \t \t %lu \n",CLOCK_extract/N_TESTS);

	for(i=0;i<MIFE_N;i++){
		mpz_clear(dy[i]);
	}
	mpz_clear(noise_gmp);
	mpz_clear(Q_gmp);
	mpz_clear(scale_M_gmp);


	return 0;
}


int main()
{

	//test_rlwe_mife_vec_vec();
	test_rlwe_mife_mat_vec();
	return 0;
}
