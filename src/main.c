#include <stdio.h>
#include <stdlib.h>

#include "fib-hashing.h"

int main() {
    int hashSize = 8;
    int tableSize = 1 << hashSize;
    int *table = calloc(tableSize, sizeof(int));
    FILE* hashingLog = fopen("log/hashing.log", "w");
    FILE* fingerprintLog = fopen("log/fingerprint.log", "w");
    FILE* collisionLog = fopen("log/collision.log", "w");

    for (int i = 0; i < 1024; i++) {
        size_t key = i;
        size_t hash = fibonacciHash(key, hashSize);
        size_t fingerprint = fibonacciFingerprint(key);

        fprintf(hashingLog, "%d: %zu\n", i, hash); // Print the hashing results
        fprintf(fingerprintLog, "%d: %zu\n", i, fingerprint); // Print the fingerprint results
        table[hash]++; // Count the number of collisions
    }

    for (int i = 0; i < tableSize; i++) {
        fprintf(collisionLog, "%d: %d\n", i, table[i]); // Print the collision counts
    }

    fclose(hashingLog);
    fclose(collisionLog);
    fclose(fingerprintLog);

    return 0;
}