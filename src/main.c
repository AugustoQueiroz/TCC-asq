// $ gcc main.c fib-hashing.c kmer-mapping.c hashTable.c
// $ ./a.out <alpha> <K>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Hashing.h"
#include "SequencingSimulator.h"
#include "HashTable.h"
#include "kmer-mapping.h"

uint8_t alreadyExplored(size_t* exploredKMers, size_t lastPos, size_t kMerCode);
uint8_t inArray(size_t value, size_t* array, size_t size);

int main(int argc, char** argv) {
    srand(0); // Seed for random number generator

    double alpha = atof(argv[1]); // Load factor received from command line
    size_t K = sscanf(argv[2], "%zu", &K); // K-mer size received from command line

    // Generate the random ("secret") sequence
    FILE* sequenceFile = fopen("log/sequence.log", "w");
    size_t sequenceLength = 1 << 20;
    char* sequence = generateRandomSequence(sequenceLength);
    fprintf(sequenceFile, "%s\n", sequence);
    fclose(sequenceFile);

    // Initialize the hash table
    struct HashTable hashTable;
    hashTable.indexSize = (size_t) ceil(sequenceLength / alpha);
    hashTable.table = calloc(hashTable.indexSize, sizeof(uint8_t));
    hashTable.hashFunction = fibonacciHash;
    hashTable.fingerprintFunction = fibonacciFingerprint;

    // Generate random reads
    printf("Generating reads and adding them to the hash table...\n");
    FILE* readsFile = fopen("log/reads.log", "w");
    size_t readLength = 100;
    size_t totalReads = ((sequenceLength / 100) * 50);
    char** initialKMers = calloc(totalReads, sizeof(char*));
    for (size_t i = 0; i < totalReads; i++) {
        initialKMers[i] = calloc(K + 1, sizeof(char));
        printf("%lu / %lu\r", i, totalReads);
        char* read = getRandomRead(sequence, sequenceLength, readLength);
        fprintf(readsFile, "%s\n", read);

        size_t prevCode = 0;
        for (size_t j = 0; j < readLength - K + 1; j++) {
            char* kmer = getKMerStartingAt(read, j, K);
            size_t kmerCode = mapKMer(kmer, K);

            insertIntoHashTable(&hashTable, kmerCode);

            if (j > 0) {
                switch (kmer[K-1]) {
                    case 'A':
                        updateHashTableWithEdges(&hashTable, prevCode, 0b1000);
                        break;
                    case 'C':
                        updateHashTableWithEdges(&hashTable, prevCode, 0b0100);
                        break;
                    case 'G':
                        updateHashTableWithEdges(&hashTable, prevCode, 0b0010);
                        break;
                    case 'T':
                        updateHashTableWithEdges(&hashTable, prevCode, 0b0001);
                        break;
                }
            } else if (j == 0) {
                strcpy(initialKMers[i], kmer);
            }

            prevCode = kmerCode;
            free(kmer);
        }

        free(read);
    }
    fclose(readsFile);

    // Explore all possible k-mers
    printf("Exploring all possible k-mers\n");
    FILE* allKMersResultFile = fopen("log/all-kmers.log", "w");
    for (size_t kmerCode = 0; kmerCode < ((size_t) 1) << (2*K); kmerCode++) {
        printf("%zu / %zu\r", kmerCode, ((size_t) 1) << (2*K));
        char* kmer = kMerFromCode(kmerCode, K);
        uint8_t outEdges = queryHashTable(&hashTable, kmerCode);
        fprintf(allKMersResultFile, "%s: %hhu\n", kmer, outEdges);
        free(kmer);
    }

    free(sequence);
    return 0;
}