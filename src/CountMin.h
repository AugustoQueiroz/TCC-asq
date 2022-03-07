#ifndef COUNTMIN_H
#define COUNTMIN_H

#include <stdlib.h>
#include <stdio.h>

#define COUNTER_MASK 0xFFFFFFFFFFFFFFF

struct DeBruijnCountMin {
    size_t W, D;
    uint64_t** table; // 4 bits for out-edges followed by 60 bit counter
    size_t** hashFunctionCoefficients; // An array of size D of pairs of numbers (A, B) that are the coefficients for the hash function
    size_t(*hashFunction)(size_t, size_t, size_t, size_t); // The hash function received the key, W, and the coefficients for that row
};

struct DeBruijnCountMin* createDeBruijnCountMinSketch(size_t W, size_t D);
void incrementDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key);
void updateDeBruijnCountMinOutEdges(struct DeBruijnCountMin* dBCM, size_t key, uint8_t outEdges);
uint64_t queryDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key);
void saveDeBruijnCountMin(struct DeBruijnCountMin* dBCM, FILE* outputFile);

#endif