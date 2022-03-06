#ifndef SEQUENCING_SIMULATOR_H
#define SEQUENCING_SIMULATOR_H

#include <stdlib.h>

char* generateRandomSequence(size_t length);
char* getRandomRead(char* sequence, size_t sequenceLength, size_t readLength);
char* getKMerStartingAt(char* read, size_t start, size_t kmerLength);

#endif