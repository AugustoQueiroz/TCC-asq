#include <stdlib.h>
#include <stdio.h>

#include "Hashing.h"
#include "mathutils.h"

/**
 * @brief hash(A, B, W, key) = (A * key + B) mod p mod W, where p is a large prime
 * 
 * @param A 
 * @param B 
 * @param indexSize 
 * @param key 
 * @return The hashing result in the range 0 to indexSize-1
 */
size_t polynomialHashing(size_t A, size_t B, size_t indexSize, size_t key) {
    return mod_sum(mod_mult(A, key, LARGE_PRIME), B, LARGE_PRIME) % indexSize;
}

size_t fibonacciHash(size_t key, size_t hashSize) {
    return (key * LARGE_FIB_NUMBER) % hashSize;
}

size_t fibonacciFingerprint(size_t key) {
    return (key * LARGE_FIB_NUMBER) >> 61;
}