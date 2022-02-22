#include <stdlib.h>
#include <stdio.h>
#include <string.h>

size_t mapKMer(char* kmer, size_t kmerLength) {
    size_t kmerCode = 0;
    for (size_t i = 0; i < kmerLength; i++) {
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