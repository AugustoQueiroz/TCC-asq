#ifndef KMER_MAPPING_H
#define KMER_MAPPING_H

#include <stdlib.h>
#include <string.h>

size_t mapKMer(char* kmer);
char* kMerFromCode(size_t kmerCode, size_t kmerLength);

#endif