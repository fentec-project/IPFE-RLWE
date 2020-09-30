#ifndef MOD_RED_H
#define MOD_RED_H

#include<stdint.h>
#include "params.h"

uint32_t add_mod_ntt(uint32_t a, uint32_t b, uint32_t sel);
uint32_t sub_mod_ntt(uint32_t a, uint32_t b, uint32_t sel);
uint32_t mul_mod_ntt(uint32_t a, uint32_t b, uint32_t sel);

uint32_t mod_prime(uint64_t m, uint32_t sel);


#endif
