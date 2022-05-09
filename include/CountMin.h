#ifndef COUNTMIN_H
#define COUNTMIN_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define COUNTER_MASK 0xFFF

struct DeBruijnCountMin {
    size_t W, D;
    uint16_t** table; // 4 bits for out-edges followed by 60 bit counter
    size_t** hashFunctionCoefficients; // An array of size D of pairs of numbers (A, B) that are the coefficients for the hash function
    size_t(*hashFunction)(size_t, size_t, size_t, size_t); // The hash function received the key, W, and the coefficients for that row
};

struct DeBruijnCountMin* createDeBruijnCountMinSketch(size_t W, size_t D);
void deleteDeBruijnCountMinSketch(struct DeBruijnCountMin* dBCM);
size_t* getHashes(struct DeBruijnCountMin* dBCM, size_t key);
void incrementDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key);
void updateDeBruijnCountMinOutEdges(struct DeBruijnCountMin* dBCM, size_t key, uint8_t outEdges);
uint16_t queryDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key);
void dumpTable(struct DeBruijnCountMin* dBCM, FILE* dumpFile);
void saveDeBruijnCountMin(struct DeBruijnCountMin* dBCM, FILE* outputFile);
struct DeBruijnCountMin* loadDeBruijnCountMin(const char* inputFilePath);
bool isMemberOfDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t kmerCode, size_t presence_threshold);

#endif
