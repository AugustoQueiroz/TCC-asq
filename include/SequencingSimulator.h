#ifndef SEQUENCING_SIMULATOR_H
#define SEQUENCING_SIMULATOR_H

#include <stdlib.h>

char* generateRandomSequence(size_t length);
char* getRandomRead(char* sequence, size_t sequenceLength, size_t readLength);

#endif