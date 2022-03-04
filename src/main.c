// $ gcc main.c fib-hashing.c kmer-mapping.c hashTable.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Hashing.h"
#include "kmer-mapping.h"
#include "CountMin.h"
#include "SequencingSimulator.h"

#define PRESENCE_THRESHOLD 25

int main(int argc, char** argv) {
    if (argc != 3) {
        printf("Usage: ./main <K> <Read Length>\n");
        return 1;
    }
    size_t K; sscanf(argv[1], "%zu", &K); // Get K from CLI
    size_t readLength; sscanf(argv[2], "%zu", &readLength); // Get read length from CLI

    // Initialize the CountMin sketch
    size_t W = 1 << 20,
           D = 8;
    struct DeBruijnCountMin* sketch = createDeBruijnCountMinSketch(1 << 20, 8);

    // Start receiving the reads
    char* read = malloc(readLength + 1); // Allocate memory for read
    while (scanf("%s", read) != EOF) {
        
        // Get the k-mers from the read
        size_t previousKMer = 0;
        uint8_t previousIsPresent = 0;
        for (size_t i = 0; i < readLength - K + 1; i++) {
            char* kmer = getKMerStartingAt(read, i, K);
            size_t kmerCode = mapKMer(kmer, K);

            // Increment the count for the k-mer
            incrementDeBruijnCountMin(sketch, kmerCode);

            // If the k-mer is considered to be present, add out edge to previous k-mer
            uint8_t kmerIsPresent = (queryDeBruijnCountMin(sketch, kmerCode) & COUNTER_MASK) > PRESENCE_THRESHOLD;
            if (i > 0 && kmerIsPresent) {
                switch (kmer[K-1]) {
                    case 'A':
                        updateDeBruijnCountMinOutEdges(sketch, previousKMer, 0b1000);
                        break;
                    case 'C':
                        updateDeBruijnCountMinOutEdges(sketch, previousKMer, 0b0100);
                        break;
                    case 'G':
                        updateDeBruijnCountMinOutEdges(sketch, previousKMer, 0b0010);
                        break;
                    case 'T':
                        updateDeBruijnCountMinOutEdges(sketch, previousKMer, 0b0001);
                        break;
                }
            }

            previousKMer = kmerCode;
            previousIsPresent = kmerIsPresent;
            free(kmer);
        }
    }

    // Print the sketch
    for (size_t i = 0; i < D; i++) {
        for (size_t j = 0; j < W; j++) {
            printf("%llu ", sketch->table[i][j]);
        }
        printf("\n");
    }

    // Ending the program
    free(read);
    return 0;
}