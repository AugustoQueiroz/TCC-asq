#ifndef SEQUENCING_SIMULATOR_H
#define SEQUENCING_SIMULATOR_H

#include <stdlib.h>

/**
 * @brief Generates a random DNA sequence of the given length.
 * 
 * @param length The length, in bases, of the sequence to generate.
 * @return char* The generated DNA sequence as a string. (Dynamically allocated, must be freed)
 */
char* generateRandomSequence(size_t length);

/**
 * @brief Get a read of size @p readLength starting at a random position in the given sequence.
 * 
 * @param sequence The DNA sequence the read is to be taken from.
 * @param sequenceLength The length of the sequence.
 * @param readLength The length of the read.
 * @return char* A substring of @p sequence of length @p readLength, starting at a random position. (Dynamically allocated, must be freed)
 */
char* getRandomRead(char* sequence, size_t sequenceLength, size_t readLength);

/**
 * @brief Gets a substring of length k (@p kmerLength) starting the given position (@p start) in the given read.
 * 
 * @param read The read from which to take the k-mer.
 * @param start The first position to be included in the k-mer.
 * @param kmerLength k; the length of the k-mer.
 * @return char* The k-mer.
 */
char* getKMerStartingAt(char* read, size_t start, size_t kmerLength);

#endif