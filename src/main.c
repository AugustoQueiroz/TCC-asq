// $ gcc main.c Hashing.c kmer-mapping.c CountMin.c SequencingSimulator.c
// $ ./a.out <K> <Read Length>
// Reads can then be passed in with line break between each read, and EOF at the end.
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "Hashing.h"
#include "KMerProcessing.h"
#include "CountMin.h"
#include "mathutils.h"


void processRead(struct DeBruijnCountMin* sketch, char* read, size_t readLength, size_t K, size_t PRESENCE_THRESHOLD, FILE* startingKMersFile) {
	size_t previousKMer = 0; // The k-mer is stored in its base 4 representation
	for (size_t i = 0; i < readLength - K + 1; i++) {
		char* kmer = getKMerStartingAt(read, i, K);
		size_t kmerCode = mapKMer(kmer);

		incrementDeBruijnCountMin(sketch, kmerCode);

		bool kmerIsPresent = (queryDeBruijnCountMin(sketch, kmerCode) & COUNTER_MASK) >= PRESENCE_THRESHOLD;
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
		} else if (kmerIsPresent) {
			fprintf(startingKMersFile, "%s\n", kmer);
		}
		previousKMer = kmerCode;
		free(kmer);
	}
}

int main(int argc, char** argv) {
    if (argc != 6) {
        printf("Usage: ./main <K> <Read Length> <W> <D> <Presence Threshold>\n");
        return 1;
    }
    size_t K; sscanf(argv[1], "%zu", &K); // Get K from CLI
    size_t readLength; sscanf(argv[2], "%zu", &readLength); // Get read length from CLI
    size_t W; sscanf(argv[3], "%zu", &W); // Get W from CLI
    size_t D; sscanf(argv[4], "%zu", &D); // Get D from CLI
    size_t PRESENCE_THRESHOLD; sscanf(argv[5], "%zu", &PRESENCE_THRESHOLD); // Get presence threshold from CLI

    // Initialize the CountMin sketch
    // W = prime_succ(1 << 16); // The width of the sketch is the smallest prime number greater than the expected number of k-mers
    struct DeBruijnCountMin* sketch = createDeBruijnCountMinSketch(W, D);

    // Start receiving the reads
    char* read = malloc(readLength + 1); // Allocate memory for read
    size_t readCount = 0;
    FILE* startingKMersFile = fopen("starting-kmers.txt", "w");
    while (scanf("%s", read) != EOF) {
        printf("Current Read: %zu\r", readCount++);

		processRead(sketch, read, readLength, K, PRESENCE_THRESHOLD, startingKMersFile);
		char* reverseComplementRead = reverseComplement(read);
		processRead(sketch, reverseComplementRead, readLength, K, PRESENCE_THRESHOLD, startingKMersFile);
		free(reverseComplementRead);
    }
    free(read);

    // Save the sketch
    printf("Saving sketch to sketch.bin\n");
    FILE* outputFile = fopen("sketch.bin", "wb");
    saveDeBruijnCountMin(sketch, outputFile);
    fclose(outputFile);

    // Ending the program
    deleteDeBruijnCountMinSketch(sketch);
    return 0;
}