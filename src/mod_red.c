#include<stdio.h>
#include<stdint.h>
#include "mod_red.h"

#define k1_q1 24
#define k2_q1 14

#define k1_q2 26
#define k2_q2 16

#define k1_q3 31
#define k2_q3 24


#define k1_q1_minus_one ( (1UL<<k1_q1)-1 )
#define k1_q2_minus_one ( (1UL<<k1_q2)-1 )
#define k1_q3_minus_one ( (1UL<<k1_q3)-1 )

uint32_t add_mod_ntt(uint32_t a, uint32_t b, uint32_t sel){

	uint64_t c;
	c=a+b;

	if(c>=MIFE_MOD_Q_I[sel]){
			c=c-MIFE_MOD_Q_I[sel];
	}

	return (uint32_t)c;

}

uint32_t sub_mod_ntt(uint32_t a, uint32_t b, uint32_t sel){	//returns a-b Mod Q

	uint64_t c;
	c=a+MIFE_MOD_Q_I[sel]-b;

	if(c>=MIFE_MOD_Q_I[sel]){
			c=c-MIFE_MOD_Q_I[sel];
	}

	return (uint32_t)c;

}

uint32_t mul_mod_ntt(uint32_t a, uint32_t b, uint32_t sel){


	uint64_t m;
	

	m=(uint64_t)a*(uint64_t)b;
	while(m>(2*MIFE_MOD_Q_I[sel]) ){
		if(sel==0){
			m=( m& (k1_q1_minus_one) ) + ( ((m>>k1_q1)<<k2_q1) - (m>>k1_q1) );
		}
		else if(sel==1){
			m=( m& (k1_q2_minus_one) ) + ( ((m>>k1_q2)<<k2_q2) - (m>>k1_q2) );
		}
		else{
			m=( m& (k1_q3_minus_one) ) + ( ((m>>k1_q3)<<k2_q3) - (m>>k1_q3) );
		}

	}
	if(m>=MIFE_MOD_Q_I[sel]){
		m=m-MIFE_MOD_Q_I[sel];
	}		


	return (uint32_t)m;
}

uint32_t mod_prime(uint64_t m, uint32_t sel){

	while(m>(2*MIFE_MOD_Q_I[sel]) ){
		if(sel==0){
			m=( m& (k1_q1_minus_one) ) + ( ((m>>k1_q1)<<k2_q1) - (m>>k1_q1) );
		}
		else if(sel==1){
			m=( m& (k1_q2_minus_one) ) + ( ((m>>k1_q2)<<k2_q2) - (m>>k1_q2) );
		}
		else{
			m=( m& (k1_q3_minus_one) ) + ( ((m>>k1_q3)<<k2_q3) - (m>>k1_q3) );
		}

	}
	if(m>=MIFE_MOD_Q_I[sel]){
		m=m-MIFE_MOD_Q_I[sel];
	}		


	return (uint32_t)m;
}
