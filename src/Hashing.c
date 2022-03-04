#include <stdlib.h>
#include <stdio.h>

#include "Hashing.h"

#define LARGE_PRIME 18446744073709551557

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
    return (A * key + B) % LARGE_PRIME % indexSize;
}

size_t incrementalFibonacciHash(size_t key, size_t hashSize, size_t increment) {
    return ((increment+1) * key * 11400714819323198485llu) % hashSize;
}

size_t fibonacciHash(size_t key, size_t hashSize) {
    return incrementalFibonacciHash(key, hashSize, 0);
}

size_t fibonacciFingerprint(size_t key) {
    return (key * 11400714819323198485llu) >> 61;
}