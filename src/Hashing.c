#include <stdlib.h>
#include <stdio.h>

#include "Hashing.h"

size_t incrementalFibonacciHash(size_t key, size_t hashSize, size_t increment) {
    return ((increment+1) * key * 11400714819323198485llu) % hashSize;
}

size_t fibonacciHash(size_t key, size_t hashSize) {
    return incrementalFibonacciHash(key, hashSize, 0);
}

size_t fibonacciFingerprint(size_t key) {
    return (key * 11400714819323198485llu) >> 61;
}