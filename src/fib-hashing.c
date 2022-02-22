#include <stdlib.h>
#include <stdio.h>

#include "fib-hashing.h"

size_t fibonacciHash(size_t key, size_t hashSize) {
    return (key * 11400714819323198485llu) >> (64 - hashSize);
}

size_t fibonacciFingerprint(size_t key) {
    return (key * 11400714819323198485llu) & 0b111;
}