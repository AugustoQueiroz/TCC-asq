#ifndef KMER_MAPPING_H
#define KMER_MAPPING_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

char* getKMerStartingAt(char* read, size_t start, size_t kmerLength);
size_t mapKMer(char* kmer);
char* kMerFromCode(size_t kmerCode, size_t kmerLength);
char* reverseComplement(char* kmer);

#endif
