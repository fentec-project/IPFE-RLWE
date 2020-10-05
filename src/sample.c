#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "sample.h"
#include "aes256ctr.h"

//#define UINT64_RAND_MAX ((RAND_MAX << 32)|RAND_MAX)

#define BUF_SIZE 45056 //buf_size=88*N/8
#define AES_ROUNDS 352//AES rounds=88*N/(8*128)

void
sample_zeros
(uint64_t a[SIFE_N])
{
	int i;
	for (i = 0; i < SIFE_N; ++i) {
		a[i] = 0;
	}
}

void
sample_ones
(uint64_t a[SIFE_N])
{
	int i;
	for (i = 0; i < SIFE_N; ++i) {
		a[i] = 1;
	}
}

void
sample_polya
(unsigned char *seed, uint32_t poly_a[SIFE_NMODULI][SIFE_N])
{
	uint64_t i, counter=0, coeff_t[2],j;

	aes256ctr_ctx state;
	aes256ctr_init(&state, seed, 0);
	unsigned char buf[BUF_SIZE];	// buffer to hold 4096 coeffs. generate 88 bits at a time then 

	unsigned char small_buf[128];

	mpz_t Q_gmp, coeff_gmp, r_gmp;

	mpz_init(Q_gmp);
	mpz_init(coeff_gmp);
	mpz_init(r_gmp);

	if(mpz_set_str(Q_gmp, SIFE_Q_str,10)!=0){
		printf("--ERROR unable to set Q to gmp--\n");
		return;
	}

	aes256ctr_squeezeblocks(buf, AES_ROUNDS, &state); //as AES generates 128 bytes at a time

	for(i=0;i<SIFE_N;i++){

		coeff_t[0]=0;
		for(j=0;j<8;j++){
			coeff_t[0]=coeff_t[0] | ((uint64_t) buf[11*i+j]<<(8*j)  );
		}
		coeff_t[1]=buf[11*i+8] | buf[11*i+9]<<8 | (buf[11*i+10] & 0x01)<<16; //81 bits

		mpz_set_ui(coeff_gmp,coeff_t[1]);
		mpz_mul_2exp(coeff_gmp, coeff_gmp, 64);
		mpz_add_ui(coeff_gmp, coeff_gmp, coeff_t[0]);

		if(mpz_cmp(coeff_gmp, Q_gmp)<=0){
			
			poly_a[0][counter]=mpz_mod_ui(r_gmp, coeff_gmp, SIFE_MOD_Q_I[0]);	//set the reduced values mod q1, q2, q3
			poly_a[1][counter]=mpz_mod_ui(r_gmp, coeff_gmp, SIFE_MOD_Q_I[1]);
			poly_a[2][counter]=mpz_mod_ui(r_gmp, coeff_gmp, SIFE_MOD_Q_I[2]);

			counter++;	//increase the counter

		}
	}

	while(counter<SIFE_N){
		aes256ctr_squeezeblocks(small_buf, 1, &state); 

		for(i=0;i<11;i++){
			coeff_t[0]=0;
			for(j=0;j<8;j++){
				coeff_t[0]=coeff_t[0] | ((uint64_t) small_buf[11*i+j]<<(8*j)  );
			}
			coeff_t[1]=small_buf[11*i+8] | small_buf[11*i+9]<<8 | (small_buf[11*i+10] & 0x01)<<16; //81 bits

			mpz_set_ui(coeff_gmp,coeff_t[1]);
			mpz_mul_2exp(coeff_gmp, coeff_gmp, 64);
			mpz_add_ui(coeff_gmp, coeff_gmp, coeff_t[0]);

			if(mpz_cmp(coeff_gmp, Q_gmp)<=0){

				poly_a[0][counter]=mpz_mod_ui(r_gmp, coeff_gmp, SIFE_MOD_Q_I[0]);	//set the reduced values mod q1, q2, q3
				poly_a[1][counter]=mpz_mod_ui(r_gmp, coeff_gmp, SIFE_MOD_Q_I[1]);
				poly_a[2][counter]=mpz_mod_ui(r_gmp, coeff_gmp, SIFE_MOD_Q_I[2]);

				counter++;	//increase the counter
				if(counter==SIFE_N)
					break;
			}
		}
	}

	mpz_clear(Q_gmp);
	mpz_clear(coeff_gmp);
	mpz_clear(r_gmp);
}

void
sample_uniform
(mpz_t a[SIFE_N])
{
	int i;
	uint64_t num;
	//const char Q_string[] ="259809622039819";
	mpz_t Q_gmp, num_gmp;
	
	mpz_init(Q_gmp);
	mpz_init(num_gmp);

	//mpz_set_str(Q_gmp, Q_string, 10);
	if(mpz_set_str(Q_gmp, SIFE_Q_str, 10)!=0){

		printf("--ERROR unable to set Q to gmp--\n");
		return;
	}

	for (i = 0; i < SIFE_N; ++i) {
		
		num = rand();
		mpz_set_ui(num_gmp, num);

		mpz_mul_2exp(num_gmp, num_gmp, 32);
		num=rand();
		mpz_add_ui(num_gmp, num_gmp, num);	

		mpz_mul_2exp(num_gmp, num_gmp, 32);
		num=rand();
		mpz_add_ui(num_gmp, num_gmp, num);

		mpz_mod(a[i], num_gmp, Q_gmp);
	}

	mpz_clear(Q_gmp);
	mpz_clear(num_gmp);
}


void
sample_m
(uint32_t a[SIFE_N])
{
	int i;
	uint64_t num;
	for (i = 0; i < SIFE_N; ++i) {
		do {
			num = rand() & (SIFE_B_x-1);	//only works for power-of-two numbers to be changed
		} while (num >= SIFE_B_x);
		a[i] = num;
	}
}

void
sample_x
(uint32_t a[SIFE_L])
{
	int i;
	uint64_t num;
	for (i = 0; i < SIFE_L; ++i) {
		do {
			num = rand() & (SIFE_B_x-1);	//only works for power-of-two numbers to be changed
		} while (num >= SIFE_B_x);
		a[i] = num;
	}
}

void
sample_y
(uint32_t a[SIFE_L])
{
	int i;
	uint64_t num;
	for (i = 0; i < SIFE_L; ++i) {
		do {
			num = rand() & (SIFE_B_y-1);
		} while (num >= SIFE_B_y);
		a[i] = num;
	}
}
