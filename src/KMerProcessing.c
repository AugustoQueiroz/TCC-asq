#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "KMerProcessing.h"

size_t mapKMer(char* kmer) {
    size_t kmerCode = 0;
    for (size_t i = 0; kmer[i] != '\0'; i++) {
        kmerCode <<= 2;
        switch (kmer[i]) {
            case 'A':
                kmerCode |= 0;
                break;
            case 'C':
                kmerCode |= 1;
                break;
            case 'G':
                kmerCode |= 2;
                break;
            case 'T':
                kmerCode |= 3;
                break;
        }
    }

    return kmerCode;
}

char* kMerFromCode(size_t kmerCode, size_t kmerLength) {
    char* kmer = calloc(kmerLength+1, sizeof(char));
    for (size_t i = 0; i < kmerLength; i++) {
        kmer[kmerLength - i - 1] = "ACGT"[kmerCode & 0b11];
        kmerCode >>= 2;
    }
    kmer[kmerLength] = '\0';

    return kmer;
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