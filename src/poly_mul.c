#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "modred.h"
#include "params.h"
#include "poly_mul.h"

void
add_mod
(const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod)
{
	int i;
	uint64_t sum;
	for (i = 0; i < MIFE_N; ++i) {
		sum = (a[i] + b[i]);
		c[i] = mod_red(sum, mod);
	}
}

void
sub_mod
(const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod)
{
	int i;
	uint64_t sub;
	for (i = 0; i < MIFE_N; ++i) {
		sub = (a[i] - b[i]);
		c[i] = mod_red(sub, mod);
	}
}

void
poly_mul_mod
(const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod)
{
	int i;
	for (i = 0; i < MIFE_N; ++i) {
		c[i] = 0;
	}
	poly_mul_mac_mod(a, b, c, mod);
}

void
poly_mul_mac_mod
(const uint32_t a[MIFE_N], const uint32_t b[MIFE_N], uint32_t c[MIFE_N], const uint32_t mod)
{
	int i, j;
	uint64_t mac;
	uint64_t res[2*MIFE_N-1] = {0};
	for (i = 0; i < MIFE_N; ++i) {
		for (j = 0; j < MIFE_N; ++j) {
			//printf("a = %ld --- b = %ld \n", a[i], b[j]);
			mac = (uint64_t)a[i]*b[j];
			//printf("mac = %ld \n", mac);
			mac = mac + res[i+j];
			//printf("mac = %ld \n", mac);
			res[i+j] = mod_red(mac, mod);
			//printf("res[%ld] = %ld \n", i+j, res[i+j]);
		}
	}
//	printf("INTERMEDIATE RESULT = ");
//		for (i = 0; i < 2*MIFE_N-1; ++i) {
//		printf("Mod(%ld,%ld)*x^%d + ",res[i],mod,i);
//	}
	for (i = MIFE_N; i < 2*MIFE_N-1; ++i) {
		mac = res[i - MIFE_N] - res[i];
		mac = mac + c[i - MIFE_N];
		c[i - MIFE_N] = mod_red(mac, mod);
	}
	mac = (c[MIFE_N-1] + res[MIFE_N-1]);
	c[MIFE_N-1] =  mod_red(mac, mod);
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