#include <stdlib.h>
#include <stdio.h>

size_t fibonacciHash(size_t key, size_t hashSize) {
    return (key * 11400714819323198485llu) >> (64 - hashSize);
}

int main() {
    int hashSize = 8;
    int tableSize = 1 << hashSize;
    int *table = calloc(tableSize, sizeof(int));
    FILE* hashingLog = fopen("log/hashing.log", "w");
    FILE* collisionLog = fopen("log/collision.log", "w");

    for (int i = 0; i < 1024; i++) {
        size_t key = i;
        size_t hash = fibonacciHash(key, hashSize);

        fprintf(hashingLog, "%d: %zu\n", i, hash); // Print the hashing results
        table[hash]++; // Count the number of collisions
    }

    for (int i = 0; i < tableSize; i++) {
        fprintf(collisionLog, "%d: %d\n", i, table[i]); // Print the collision counts
    }

    fclose(hashingLog);
    fclose(collisionLog);

    return 0;
}