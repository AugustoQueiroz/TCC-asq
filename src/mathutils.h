#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define IS_EVEN( N ) ( ((N) % 2) == 0 )
#define DIVCEIL(NUM, DEN) (ceil(((double)(NUM))/((double)(DEN))))

uint64_t mod_sum(uint64_t a, uint64_t b, uint64_t m);
uint64_t mod_mult(uint64_t a, uint64_t b, uint64_t m);
uint64_t mod_pow(uint64_t b, uint64_t e, uint64_t m);
bool naive_is_prime(uint64_t val);
bool is_prime(uint64_t n);
uint64_t prime_succ(uint64_t n);

#endif