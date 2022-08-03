#include <assert.h>

#include "mathutils.h"

/**
 * @brief Modular sum of a and b modulo m. Guaranteed to not be clipped by overflow.
 * 
 * @param a 
 * @param b 
 * @param m 
 * @return uint64_t 
 */
uint64_t mod_sum(uint64_t a, uint64_t b, uint64_t m)
{
	a %= m;
	b %= m;
	if ( a <= UINT64_MAX - b) {
		return (a + b) % m;
	} else {
		return b - ( m - a);
	}

}

/**
 * @brief Modular multiplication of a and b modulo m. Guaranteed to not be clipped by overflow.
 * 
 * @param a 
 * @param b 
 * @param m 
 * @return uint64_t 
 */
uint64_t mod_mult(uint64_t a, uint64_t b, uint64_t m)
{
	uint64_t res = 0;

    // Make b the smaller number
    if (a < b) {
        uint64_t tmp = a;
        a = b;
        b = tmp;
    }

	a = a % m;
	while (b > 0) {
		if ( b & 1 ) { // b is odd
			res = mod_sum(res, a, m); // res = res + a mod m
		}
		a = mod_sum(a, a, m); // a = 2*a mod m
		b /= 2;
	}
	return res % m;
}

/**
 * @brief Modular potentiation of a to the power of b modulo m. Guaranteed to not be clipped by overflow.
 * 
 * @param b 
 * @param e 
 * @param m 
 * @return uint64_t 
 */
uint64_t mod_pow(uint64_t b, uint64_t e, uint64_t m)
{
	b = b % m;
	if (b == 0) return 0;
	uint64_t res = 1;
	while (e) {
		if (e & 1)
			res = mod_mult(res, b, m); // res = (res*b) % m;
		e >>= 1;
		b = mod_mult(b, b, m); // b = b^2 % m;
	}
	return res;
}

/**
 * @brief Determines whether or not n is prime.
 * 
 * @param n 
 * @return true 
 * @return false 
 */
bool is_prime(uint64_t n)
{
	uint64_t a[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
	uint64_t limits[] = {2046llu, 1373652llu, 25326000llu, 3215031750llu, 2152302898746llu, 3474749660382llu, 341550071728320llu, 3825123056546413050llu, 18446744073709551615llu};
	size_t nwitness[] = {1, 2, 3, 4, 5, 6, 7, 9, 12};
	if (n < 2) return false;
	if (n == 2) return true;
	if (IS_EVEN(n)) return false;

	size_t nwit = 1;
	for (size_t j = 0; j < 9 && n >= limits[j]; nwit = nwitness[++j]);

	uint64_t d = n-1;
	uint64_t r = 0;
	while (IS_EVEN(d)) {
		d >>= 1;
		r += 1;
	}

	bool prime = true;
	for (size_t i = 0; prime && i < nwit; i++) {
		prime = false;
		uint64_t x = mod_pow(a[i], d, n); // (a[i]^d) mod n
		if ( x == 1 || x == (n-1) ) {
			prime = true;
			continue;
		}
		for (size_t _j=0; _j < r-1; _j++) {
			x = mod_mult(x, x, n); // (x*x) % n;
			if ( x == n-1 ) {
				prime = true;
				break;
			}
		}
	}
	return prime;
}

/**
 * @brief Finds the prime successor of n.
 * 
 * @param n 
 * @return uint64_t 
 */
uint64_t prime_succ(uint64_t n)
{
	uint64_t firstfew[] = {2,2,2,3,5,5,7};
	if (n<=6) return firstfew[n];
	uint64_t k = (uint64_t) DIVCEIL(n, 6);
	assert (n <= k*6);
	uint64_t ret = 6 * k - 1;
	bool pm = true;
	if (ret < n) {
		ret += 2;
		pm = false;
	}
	while ( !is_prime(ret) ) {
		ret += (pm)?2:4;
		pm = !pm;
	}
	return ret;
}