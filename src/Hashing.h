#ifndef FIB_HASHING_H
#define FIB_HASHING_H

#include <stdlib.h>

size_t incrementalFibonacciHash(size_t key, size_t indexSize, size_t increment);
size_t fibonacciHash(size_t key, size_t table_size);
size_t fibonacciFingerprint(size_t key);

#endif