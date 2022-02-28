#include <stdlib.h>
#include <stdio.h>

#include "SequencingSimulator.h"

#define KMER_LENGTH 32
#define COUNTER_MASK 0xFFFFFFFFFFFFFFF
#define PRESENCE_THRESHOLD 25

char* generateRandomSequence(size_t length) {
    char bases[] = "ACGT";
    char* sequence = malloc(length + 1);
    for (size_t i = 0; i < length; i++) {
        sequence[i] = bases[rand() % 4];
    }
    sequence[length] = '\0';
    return sequence;
}

char* getRandomRead(char* sequence, size_t sequenceLength, size_t readLength) {
    char* read = malloc(readLength + 1);
    size_t start = rand() % (sequenceLength - readLength);
    for (size_t i = 0; i < readLength; i++) {
        read[i] = sequence[start + i];
    }
    read[readLength] = '\0';
    return read;
}

char* getKMerStartingAt(char* read, size_t start, size_t kmerLength) {
    char* kmer = malloc(kmerLength + 1);
    for (size_t i = 0; i < kmerLength; i++) {
        kmer[i] = read[start + i];
    }
    kmer[kmerLength] = '\0';
    return kmer;
}