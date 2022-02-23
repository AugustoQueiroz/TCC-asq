// $ gcc main.c fib-hashing.c kmer-mapping.c hashTable.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fib-hashing.h"
#include "kmer-mapping.h"
#include "hashTable.h"

char* getKMerStartingAt(char* sequence, size_t start, size_t kmerLength) {
    char* kmer = malloc(kmerLength + 1);
    for (size_t i = 0; i < kmerLength; i++) {
        kmer[i] = sequence[start + i];
    }
    kmer[kmerLength] = '\0';
    return kmer;
}

int main() {
    // // k-mer mapping example/test
    // char* kmer = "ACGTGACATGACATAGCAGACATTA";
    // size_t kmerLength = strlen(kmer);
    // size_t kmerCode = mapKMer(kmer, kmerLength);
    
    // for (size_t i = 0; i < kmerLength; i++) {
    //     printf("%c: %lx\n", kmer[i], kmerCode >> (2 * (kmerLength - i - 1)) & 0x3);
    // }

    // // Fibonacci hashing/fingerprint example/test
    // int hashSize = 8;
    // int tableSize = 1 << hashSize;
    // int *table = calloc(tableSize, sizeof(int));
    // FILE* hashingLog = fopen("log/hashing.log", "w");
    // FILE* fingerprintLog = fopen("log/fingerprint.log", "w");
    // FILE* collisionLog = fopen("log/collision.log", "w");

    // for (int i = 0; i < 1024; i++) {
    //     size_t key = i;
    //     size_t hash = fibonacciHash(key, hashSize);
    //     size_t fingerprint = fibonacciFingerprint(key);

    //     fprintf(hashingLog, "%d: %zu\n", i, hash); // Print the hashing results
    //     fprintf(fingerprintLog, "%d: %zu\n", i, fingerprint); // Print the fingerprint results
    //     table[hash]++; // Count the number of collisions
    // }

    // for (int i = 0; i < tableSize; i++) {
    //     fprintf(collisionLog, "%d: %d\n", i, table[i]); // Print the collision counts
    // }

    // fclose(hashingLog);
    // fclose(collisionLog);
    // fclose(fingerprintLog);

    // struct HashTable hashTable;
    // hashTable.indexSize = 8;
    // hashTable.table = calloc(hashTable.indexSize, sizeof(size_t));
    // hashTable.hashFunction = fibonacciHash;
    // hashTable.fingerprintFunction = fibonacciFingerprint;
    
    // for (int i = 0; i < 1024; i++) {
    //     size_t key = i;
    //     size_t value = 1 << 7;
    //     printf("Inserting %zu into table\n", key);
    //     insertIntoHashTable(&hashTable, key);
    // }

    // updateHashTableWithEdges(&hashTable, 0, 0b0101);
    // updateHashTableWithEdges(&hashTable, 1, 0b1001);

    double alpha = 0;
    scanf("%lf", &alpha);

    // Generate random reading
    char bases[] = "ACGT";
    size_t readLength = 1024;
    char* read = calloc(readLength, sizeof(char));
    for (int i = 0; i < readLength; i++) {
        read[i] = bases[rand() % 4];
    }

    // Read k-mers
    int kmerLength = 8;

    struct HashTable hashTable;
    size_t tableSize = ceil(readLength*(1/alpha));
    printf("Table size: %zu\n", tableSize);
    hashTable.indexSize = tableSize;
    hashTable.table = calloc(hashTable.indexSize, sizeof(uint8_t));
    hashTable.hashFunction = fibonacciHash;
    hashTable.fingerprintFunction = fibonacciFingerprint;
    FILE* readResultsLog = fopen("log/read-results.log", "w");

    size_t prevKMerCode = 0;
    for (int i = 0; i < readLength - kmerLength; i++) {
        char* kmer = getKMerStartingAt(read, i, kmerLength);
        size_t kmerCode = mapKMer(kmer, kmerLength);
        fprintf(readResultsLog, "%s\n", kmer);

        insertIntoHashTable(&hashTable, kmerCode);
        if (i > 0) {
            uint8_t outgoingEdge = 0;
            switch(kmer[kmerLength-1]) {
                case 'A':
                    outgoingEdge = 0b1000;
                    break;
                case 'C':
                    outgoingEdge = 0b0100;
                    break;
                case 'G':
                    outgoingEdge = 0b0010;
                    break;
                case 'T':
                    outgoingEdge = 0b0001;
                    break;
                default:
                    printf("Error: Invalid base\n");
                    break;
            }
            updateHashTableWithEdges(&hashTable, prevKMerCode, outgoingEdge);
        }
        prevKMerCode = kmerCode;
    }

    // Try all possible k-mers
    FILE* allKMersLog = fopen("log/all-kmers.log", "w");
    for (size_t kmerCode = 0; kmerCode < (1 << (2*kmerLength+1)); kmerCode++) {
        char* kmer = kMerFromCode(kmerCode, kmerLength);
        fprintf(allKMersLog, "%s: %zu - %zu - %zu - %hhu\n", kmer, kmerCode, fibonacciHash(kmerCode, hashTable.indexSize), fibonacciFingerprint(kmerCode), queryHashTable(&hashTable, kmerCode));
    }

    return 0;
}