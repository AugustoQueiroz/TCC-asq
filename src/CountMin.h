#ifndef COUNTMIN_H
#define COUNTMIN_H

#include <stdlib.h>

struct DeBruijnCountMin {
    size_t W, D;
    uint64_t** table; // 4 bits for out-edges followed by 60 bit counter
    size_t(*hashFunctionWithLevel)(size_t, size_t, size_t);
};

void incrementDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key);
void updateDeBruijnCountMinOutEdges(struct DeBruijnCountMin* dBCM, size_t key, uint8_t outEdges);
uint64_t queryDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key);

#endif