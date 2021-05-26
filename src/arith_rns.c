#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "arith_rns.h"

#if SEC_LEVEL==0

/*
#define k1_q1 16
#define k2_q1 12

#define k1_q2 17
#define k2_q2 14

#define k1_q3 31
#define k2_q3 17
*/
#define k1_q1 14
#define k2_q1 12

#define k1_q2 23
#define k2_q2 17

#define k1_q3 29
#define k2_q3 18

#elif SEC_LEVEL==1

#define k1_q1 24
#define k2_q1 14

//#define k1_q2 26
//#define k2_q2 16

#define k1_q2 31
#define k2_q2 17

#define k1_q3 31
#define k2_q3 24


#elif SEC_LEVEL==2

/*
#define k1_q1 31
#define k2_q1 17

#define k1_q2 31
#define k2_q2 19

#define k1_q3 32
#define k2_q3 20
*/

#define k1_q1 17
#define k2_q1 14

#define k1_q2 20
#define k2_q2 14

#define k1_q3 32
#define k2_q3 20

#define k1_q4 32
#define k2_q4 30

#endif



#define k1_q1_minus_one ( (1UL<<k1_q1)-1 )
#define k1_q2_minus_one ( (1UL<<k1_q2)-1 )
#define k1_q3_minus_one ( (1UL<<k1_q3)-1 )

#if SEC_LEVEL==2
#define k1_q4_minus_one ( (1UL<<k1_q4)-1 )
#endif

uint32_t
add_mod_ntt
(uint32_t a, uint32_t b, uint32_t sel)
{
	uint64_t c;

	c = (uint64_t)a + (uint64_t)b;

	if (c >= SIFE_MOD_Q_I[sel]) {
		c -= SIFE_MOD_Q_I[sel];
	}
	return (uint32_t)c;
}

uint32_t
sub_mod_ntt
(uint32_t a, uint32_t b, uint32_t sel)//returns a-b Mod Q
{
	uint64_t c;

	c = (uint64_t)a + (uint64_t)SIFE_MOD_Q_I[sel] - (uint64_t)b;

	if (c >= SIFE_MOD_Q_I[sel]) {
		c -= SIFE_MOD_Q_I[sel];
	}
	return (uint32_t)c;
}

uint32_t
mul_mod_ntt
(uint32_t a, uint32_t b, uint32_t sel)
{
	uint64_t m;

	m = (uint64_t)a*(uint64_t)b;

	return (uint32_t)(mod_prime(m,sel));
}

uint32_t
mod_prime
(uint64_t m, uint32_t sel)

{

	#if SEC_LEVEL==2


		while ( m > (2*(uint64_t)SIFE_MOD_Q_I[sel]) ) {
			if (sel == 0) {
				m = ( m& (k1_q1_minus_one) ) + ( ((m>>k1_q1)<<k2_q1) - (m>>k1_q1) );
			}
			else if (sel == 1) {
				m = ( m& (k1_q2_minus_one) ) + ( ((m>>k1_q2)<<k2_q2) - (m>>k1_q2) );
			}
			else if (sel == 2){
				m = ( m& (k1_q3_minus_one) ) + ( ((m>>k1_q3)<<k2_q3) - (m>>k1_q3) );
			}
			else{
				m = ( m& (k1_q4_minus_one) ) + ( ((m>>k1_q4)<<k2_q4) - (m>>k1_q4) );
			}
		}

	#else
		while ( m > (2*(uint64_t)SIFE_MOD_Q_I[sel]) ) {
			if (sel == 0) {
				m = ( m& (k1_q1_minus_one) ) + ( ((m>>k1_q1)<<k2_q1) - (m>>k1_q1) );
			}
			else if (sel == 1) {
				m = ( m& (k1_q2_minus_one) ) + ( ((m>>k1_q2)<<k2_q2) - (m>>k1_q2) );
			}
			else{
				m = ( m& (k1_q3_minus_one) ) + ( ((m>>k1_q3)<<k2_q3) - (m>>k1_q3) );
			}
		}
	#endif

	if (m >= SIFE_MOD_Q_I[sel]) {
		m = m - SIFE_MOD_Q_I[sel];
	}
	return (uint32_t)m;
}

void
poly_add_mod
(const uint32_t a[SIFE_N], const uint32_t b[SIFE_N], uint32_t c[SIFE_N], const uint32_t sel)
{
	int i;
	for (i = 0; i < SIFE_N; ++i) {
		c[i] = add_mod_ntt(a[i], b[i], sel);
	}
}

void
poly_sub_mod
(const uint32_t a[SIFE_N], const uint32_t b[SIFE_N], uint32_t c[SIFE_N], const uint32_t sel)
{
	int i;
	for (i = 0; i < SIFE_N; ++i) {
		c[i] = sub_mod_ntt(a[i], b[i], sel);
	}
}
