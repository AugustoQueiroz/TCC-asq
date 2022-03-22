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

/**
 * @brief hash(key) = (key * LARGE_FIB_NUMBER) mod hashSize.
 * 
 * @param key 
 * @param hashSize 
 * @return size_t An integer between 0 and hashSize-1
 */
size_t fibonacciHash(size_t key, size_t hashSize) {
    return (key * LARGE_FIB_NUMBER) % hashSize;
}

// TODO - Update the fingerprint method to take the number of bits of the output as a parameter
/**
 * @brief The fibonacci fingerprint of a number are the 3 most significant bits of the fibonacci hash of the number before the mod / shift operation.
 * 
 * @param key 
 * @return size_t 
 */
size_t fibonacciFingerprint(size_t key) {
    return (key * LARGE_FIB_NUMBER) >> 61;
}