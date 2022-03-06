#include <stdlib.h>
#include <stdio.h>

#include "SequencingSimulator.h"

#define KMER_LENGTH 32
#define COUNTER_MASK 0xFFFFFFFFFFFFFFF
#define PRESENCE_THRESHOLD 25

/**
 * @brief Generates a random DNA sequence of the given length.
 * 
 * @param length The length, in bases, of the sequence to generate.
 * @return char* The generated DNA sequence as a string. (Dynamically allocated, must be freed)
 */
char* generateRandomSequence(size_t length) {
    char bases[] = "ACGT";
    char* sequence = malloc(length + 1);
    for (size_t i = 0; i < length; i++) {
        sequence[i] = bases[rand() % 4];
    }
    sequence[length] = '\0';
    return sequence;
}

/**
 * @brief Get a read of size @p readLength starting at a random position in the given sequence.
 * 
 * @param sequence The DNA sequence the read is to be taken from.
 * @param sequenceLength The length of the sequence.
 * @param readLength The length of the read.
 * @return char* A substring of @p sequence of length @p readLength, starting at a random position. (Dynamically allocated, must be freed)
 */
char* getRandomRead(char* sequence, size_t sequenceLength, size_t readLength) {
    char* read = malloc(readLength + 1);
    size_t start = rand() % (sequenceLength - readLength);
    for (size_t i = 0; i < readLength; i++) {
        read[i] = sequence[start + i];
    }
    read[readLength] = '\0';
    return read;
}

/**
 * @brief Gets a substring of length k (@p kmerLength) starting the given position (@p start) in the given read.
 * 
 * @param read The read from which to take the k-mer.
 * @param start The first position to be included in the k-mer.
 * @param kmerLength k; the length of the k-mer.
 * @return char* The k-mer.
 */
char* getKMerStartingAt(char* read, size_t start, size_t kmerLength) {
    char* kmer = malloc(kmerLength + 1);
    for (size_t i = 0; i < kmerLength; i++) {
        kmer[i] = read[start + i];
    }
    kmer[kmerLength] = '\0';
    return kmer;
}