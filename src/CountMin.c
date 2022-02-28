#include "CountMin.h"

/**
 * @brief Get the hashes for a given key for all the rows in the CountMin sketch
 * 
 * @param dBCM A pointer to the de Bruijn CountMin sketch.
 * @param key The key whose hashes are to be calculated.
 * @return size_t* An array of hashes, one for each row in the CountMin sketch.
 */
size_t* getHashes(struct DeBruijnCountMin* dBCM, size_t key) {
    size_t* hashes = malloc(dBCM->D * sizeof(size_t));
    for (size_t i = 0; i < dBCM->D; i++) {
        hashes[i] = dBCM->hashFunctionWithLevel(key, dBCM->W, i);
    }
    return hashes;
}

/**
 * @brief Increment the count for a given key in the CountMin sketch.
 * 
 * @param dBCM A pointer to the de Bruijn CountMin sketch.
 * @param key The key whose count is to be incremented.
 */
void incrementDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key) {
    size_t* hashes = getHashes(dBCM, key);
    for (size_t i = 0; i < dBCM->D; i++) {
        size_t hash = hashes[i];
        dBCM->table[i][hash]++;
    }
    free(hashes);
}

/**
 * @brief Update the out-edge information for the given key.
 * 
 * @param dBCM A pointer to the de Bruijn CountMin sketch.
 * @param key The key whose out-edge information is to be updated.
 * @param outEdges The new out-edge information to be added.
 */
void updateDeBruijnCountMinOutEdges(struct DeBruijnCountMin* dBCM, size_t key, uint8_t outEdges) {
    size_t* hashes = getHashes(dBCM, key);
    for (size_t i = 0; i < dBCM->D; i++) {
        size_t hash = hashes[i];
        dBCM->table[i][hash] |= ((uint64_t) outEdges) << 60;
    }
    free(hashes);
}

/**
 * @brief Query the de Bruijn CountMin sketch for the given key.
 * 
 * @param dBCM A pointer to the de Bruijn CountMin sketch.
 * @param key The key to be queried
 * @return uint64_t A 64-bit integer representing the count as well as the out-edges for the given key. The 4 most significant bits represent the out-edges {A, C, G, T}, and the other 60 bits are the counter.
 */
uint64_t queryDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key) {
    size_t* hashes = getHashes(dBCM, key);
    uint64_t count = 0;
    uint64_t outEdges = -1;
    for (size_t i = 0; i < dBCM->D; i++) {
        size_t hash = hashes[i];
        count = dBCM->table[i][hash] < count ? dBCM->table[i][hash] : count; // Get the minimum count
        outEdges &= dBCM->table[i][hash] >> 60; // And the intersection of out-edges
    }
    free(hashes);
    return count | (outEdges << 60);
}