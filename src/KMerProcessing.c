#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "KMerProcessing.h"

/**
 * @brief Gets a substring of length k (@p kmerLength) starting the given position (@p start) in the given read.
 * 
 * @param read The read from which to take the k-mer.
 * @param start The first position to be included in the k-mer.
 * @param kmerLength k; the length of the k-mer.
 * @return char* The k-mer.
 */
char* getKMerStartingAt(char* read, size_t start, size_t kmerLength) {
    char* kmer = (char*) malloc(kmerLength + 1);
    for (size_t i = 0; i < kmerLength; i++) {
        kmer[i] = read[start + i];
    }
    kmer[kmerLength] = '\0';
    return kmer;
}

/**
 * @brief Converts the k-mer string to a number. A k-mer can be thought of as a base-4 number, with A=0, C=1, G=2, T=3.
 * 
 * @param kmer The k-mer string.
 * @return size_t The number it encodes.
 */
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

/**
 * @brief Generates the k-mer string from the k-mer code (the number encoded by the string). Turns the number into a base-4 string, with A=0, C=1, G=2, T=3.
 * 
 * @param kmerCode Number representing the k-mer
 * @param kmerLength The length of the k-mer
 * @return char* A string of size k, ending with '\0', representing the k-mer
 */
char* kMerFromCode(size_t kmerCode, size_t kmerLength) {
    char* kmer = (char*) calloc(kmerLength+1, sizeof(char));
    for (size_t i = 0; i < kmerLength; i++) {
        kmer[kmerLength - i - 1] = "ACGT"[kmerCode & 0b11];
        kmerCode >>= 2;
    }
    kmer[kmerLength] = '\0';

    return kmer;
}

/**
 * @brief Generates the reverse complement of the given k-mer. I.e. the k-mer is reversed and each base is substituted by it's complement (A = T, C = G).
 * 
 * @param kmer
 * @return char* The reverse complement of the given k-mer. Dynamically allocated, **must be freed**.
 */
char* reverseComplement(char* kmer) {
    size_t kmerLength = strlen(kmer);
    char* reverseComplementKMer = (char*) malloc((kmerLength + 1) * sizeof(char));

    char complements[26];
    complements['A' - 'A'] = 'T';
    complements['C' - 'A'] = 'G';
    complements['G' - 'A'] = 'C';
    complements['T' - 'A'] = 'A';

    for (size_t i = kmerLength; i > 0; i--) {
        reverseComplementKMer[kmerLength - i] = complements[kmer[i-1] - 'A'];
    }

    reverseComplementKMer[kmerLength] = '\0';

    return reverseComplementKMer;
}