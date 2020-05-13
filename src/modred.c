#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "modred.h"

#include <gmp.h>

uint32_t
mod_red
(uint64_t a, uint32_t mod)
{
	uint64_t mul;
	uint32_t minus_one = mod-1;

	uint32_t result = 0;

	if ( ((a>>63) & 1) == 1 ) {
		a = a*(-1);
		mul = (uint32_t)(a % (uint64_t)mod);
		mul = mul * minus_one;
		result = (uint32_t)(mul % (uint64_t)mod);
	} else {
		result = (uint32_t)(a % (uint64_t)mod);		
	}

	return result;
}

uint32_t
mod_red_gmp
(uint64_t a, uint32_t mod)
{
	uint32_t result;
	mpz_t aa, mmod;
	mpz_init(aa);
	mpz_init(mmod);
	mpz_set_ui(aa, a);
	mpz_set_ui(mmod, mod);
	mpz_mod(aa, aa, mmod);
	result = mpz_get_ui(aa);
	mpz_clear(aa);
	mpz_clear(mmod);
	return result;
}
