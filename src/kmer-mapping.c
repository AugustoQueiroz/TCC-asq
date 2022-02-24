#include <stdlib.h>
#include <stdio.h>
#include <string.h>

size_t mapKMer(char* kmer, size_t kmerLength) {
    size_t kmerCode = 0;
    for (size_t i = 0; kmer[i] != '\0'; i++) {
        kmerCode <<= 2;
        if (kmer[i] == 'A') {
            kmerCode += 0;
        } else if (kmer[i] == 'C') {
            kmerCode += 1;
        } else if (kmer[i] == 'G') {
            kmerCode += 2;
        } else if (kmer[i] == 'T') {
            kmerCode += 3;
        } else {
            // Error
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