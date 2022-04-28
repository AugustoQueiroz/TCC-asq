#include <stdio.h>
#include <stdlib.h>

#include "SequencingSimulator.h"

int main(int argc, char** argv) {
    if (argc != 4) {
        printf("Usage: ./main <Sequence Length> <Read Length> <Coverage>\n");
        return 1;
    }

    size_t sequenceLength; sscanf(argv[1], "%zu", &sequenceLength); // Get sequence length from CLI
    size_t readLength; sscanf(argv[2], "%zu", &readLength); // Get read length from CLI
    size_t coverage; sscanf(argv[3], "%zu", &coverage); // Get coverage from CLI

    // Generate sequence
    char* sequence = generateRandomSequence(sequenceLength);

    // Log sequence
    FILE* sequenceFile = fopen("sequence.txt", "w");
    fprintf(sequenceFile, "%s", sequence);
    fclose(sequenceFile);

    // Generate reads (and log them immediately)
    size_t numberOfReads = sequenceLength * coverage / readLength;
    FILE* readsFile = fopen("reads.txt", "w");
    for (size_t i = 0; i < numberOfReads; i++) {
        char* read = getRandomRead(sequence, sequenceLength, readLength, rand() & 1);
        fprintf(readsFile, "%s\n", read);
        free(read);
    }
    fclose(readsFile);

    // Closing statements
    free(sequence);

    return 0;
}