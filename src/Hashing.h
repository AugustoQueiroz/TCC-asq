#ifndef FIB_HASHING_H
#define FIB_HASHING_H

#include <stdlib.h>

#define LARGE_PRIME 18446744073709551557llu
#define LARGE_FIB_NUMBER 11400714819323198485llu

size_t polynomialHashing(size_t A, size_t B, size_t indexSize, size_t key);
size_t fibonacciHash(size_t key, size_t table_size);
size_t fibonacciFingerprint(size_t key);

#endif