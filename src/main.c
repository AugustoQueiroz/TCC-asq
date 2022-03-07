// $ gcc main.c Hashing.c kmer-mapping.c CountMin.c SequencingSimulator.c
// $ ./a.out <K> <Read Length>
// Reads can then be passed in with line break between each read, and EOF at the end.
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
           D = 16;
    struct DeBruijnCountMin* sketch = createDeBruijnCountMinSketch(W, D);

    // Start receiving the reads
    char* read = malloc(readLength + 1); // Allocate memory for read
    while (scanf("%s", read) != EOF) {
        
        // Get the k-mers from the read
        size_t previousKMer = 0;
        uint8_t previousIsPresent = 0;
        for (size_t i = 0; i < readLength - K + 1; i++) {
            char* kmer = getKMerStartingAt(read, i, K);
            size_t kmerCode = mapKMer(kmer);

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

    // Save the sketch
    printf("Saving sketch to sketch.bin\n");
    FILE* outputFile = fopen("sketch.bin", "wb");
    saveDeBruijnCountMin(sketch, outputFile);

    // // Test all k-mers
    // printf("Exploring all possible k-mers\n");
    // FILE* allKMersResultFile = fopen("log/all-kmers.log", "w");
    // for (size_t kmerCode = 0; kmerCode < ((size_t) 1) << (2*K); kmerCode++) {
    //     printf("%zu / %zu\r", kmerCode, ((size_t) 1) << (2*K));
    //     char* kmer = kMerFromCode(kmerCode, K);
    //     uint64_t queryResult = queryDeBruijnCountMin(sketch, kmerCode);
    //     uint8_t outEdges = queryResult >> 60;
    //     uint64_t count = queryResult & COUNTER_MASK;
    //     if (count < PRESENCE_THRESHOLD) {
    //         outEdges = -1;
    //     }
    //     fprintf(allKMersResultFile, "%s: %hhu\n", kmer, outEdges);
    //     free(kmer);
    // }

    // Ending the program
    free(read);
    return 0;
}