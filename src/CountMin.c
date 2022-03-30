#include <time.h>

#include "CountMin.h"
#include "Hashing.h"

/**
 * @brief Create a new de Bruijn CountMin sketch
 * 
 * @param W The width of the sketch (i.e. the number of columns)
 * @param D The depth of the sketch (i.e. the number of rows, also the number of hash functions used)
 * @return A pointer to the new de Bruijn CountMin sketch (dynamically allocated, must be freed!)
 */
struct DeBruijnCountMin* createDeBruijnCountMinSketch(size_t W, size_t D) {
    srand(time(NULL));
    struct DeBruijnCountMin* sketch = malloc(sizeof(struct DeBruijnCountMin));
    sketch->W = W;
    sketch->D = D;
    sketch->table = malloc(sizeof(uint16_t*) * W);
    sketch->hashFunctionCoefficients = malloc(D * sizeof(size_t*));
    sketch->hashFunction = polynomialHashing;
    for (size_t i = 0; i < D; i++) {
        // Allocate and initialize the rows of the table to all 0s
        sketch->table[i] = calloc(W, sizeof(uint16_t));

        // Set the hash function coefficients
        sketch->hashFunctionCoefficients[i] = malloc(2 * sizeof(size_t));
        sketch->hashFunctionCoefficients[i][0] = rand() % LARGE_PRIME;
        sketch->hashFunctionCoefficients[i][0] = sketch->hashFunctionCoefficients[i][0] == 0 ? 1 : sketch->hashFunctionCoefficients[i][0];
        sketch->hashFunctionCoefficients[i][1] = rand() % LARGE_PRIME;
    }
    return sketch;
}

/**
 * @brief Frees the memory allocated to the de Bruijn CountMin sketch
 * 
 * @param sketch 
 */
void deleteDeBruijnCountMinSketch(struct DeBruijnCountMin* sketch) {
    for (size_t i = 0; i < sketch->D; i++) {
        free(sketch->table[i]);
        free(sketch->hashFunctionCoefficients[i]);
    }
    free(sketch->table);
    free(sketch->hashFunctionCoefficients);
    free(sketch);
}

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
        hashes[i] = dBCM->hashFunction(dBCM->hashFunctionCoefficients[i][0], dBCM->hashFunctionCoefficients[i][1], dBCM->W, key);
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
        dBCM->table[i][hashes[i]]++;
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
        dBCM->table[i][hashes[i]] |= ((uint16_t) outEdges) << 12;
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
uint16_t queryDeBruijnCountMin(struct DeBruijnCountMin* dBCM, size_t key) {
    size_t* hashes = getHashes(dBCM, key);
    uint16_t count = 0xFFFF;
    uint16_t outEdges = 0xF;
    for (size_t i = 0; i < dBCM->D; i++) {
        count = dBCM->table[i][hashes[i]] < count ? dBCM->table[i][hashes[i]] : count; // Get the minimum count
        outEdges &= (dBCM->table[i][hashes[i]] >> 12); // And the intersection of out-edges
    }
    free(hashes);
    return count | (outEdges << 12);
}

/**
 * @brief Writes the counter table to a file.
 * 
 * @param dBCM The sketch to be printed.
 * @param dumpFile The file that the table should be written to.
 */
void dumpTable(struct DeBruijnCountMin* dBCM, FILE* dumpFile) {
    for (size_t i = 0; i < dBCM->D; i++) {
        for (size_t j = 0; j < dBCM->W; j++) {
            fprintf(dumpFile, "%hu ", dBCM->table[i][j]);
        }
        fprintf(dumpFile, "\n");
    }
}

/**
 * @brief Save the de Bruijn CountMin sketch to disk
 * 
 * @param dBCM A pointer to the de Bruijn CountMin sketch.
 * @param outputFile The file pointer to which the sketch is to be written.
 */
void saveDeBruijnCountMin(struct DeBruijnCountMin* dBCM, FILE* outputFile) {
    fwrite(&dBCM->W, sizeof(size_t), 1, outputFile);
    fwrite(&dBCM->D, sizeof(size_t), 1, outputFile);
    //fwrite(&dBCM->hashFunction, sizeof(void*), 1, outputFile);
    for (size_t i = 0; i < dBCM->D; i++) {
        fwrite(dBCM->hashFunctionCoefficients[i], sizeof(size_t), 2, outputFile);
    }
    for (size_t i = 0; i < dBCM->D; i++) {
        fwrite(dBCM->table[i], sizeof(uint16_t), dBCM->W, outputFile);
    }
}