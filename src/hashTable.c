#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fib-hashing.h"
#include "hashTable.h"

void insertIntoHashTable(struct HashTable* hashTable, size_t key) {
    size_t hash = hashTable->hashFunction(key, hashTable->indexSize);
    size_t fingerprint = hashTable->fingerprintFunction(key);

    size_t index = hash;
    while (hashTable->table[index] & (1 << 7)) { // Dirty bit is set
        if (((hashTable->table[index] >> 4) & 0b111) == fingerprint) {
            // The fingerprint matches, so element is already in hashtable
            // printf("Element is already in hashtable\n");
            return;
        }
        index = (index + 1) % (hashTable->indexSize);
        
        if (index == hash) { // Looped around the table
            // Table is full
            // printf("Table is full, %zu not inserted\n", key);
            return;
        }
    }

    size_t value = (1 << 7) | (fingerprint << 4);
    // printf("%zu: %zu\n", index, value);
    hashTable->table[index] = value;
}

void updateHashTableWithEdges(struct HashTable* hashTable, size_t key, uint8_t outEdges) {
    size_t hash = hashTable->hashFunction(key, hashTable->indexSize);
    size_t fingerprint = hashTable->fingerprintFunction(key);

    size_t index = hash;
    while ((hashTable->table[index] & (1 << 7)) && ((hashTable->table[index] >> 4) & 0b111) != fingerprint) { // Dirty bit is set but fingerprint doesn't match
        index = (index + 1) % (hashTable->indexSize);
        
        if (index == hash) { // Looped around the table
            // Table is full
            // printf("Item is not in table\n");
            return;
        }
    }

    if (!(hashTable->table[index] & (1 << 7))) { // Dirty bit is not set
        // printf("Item is not in table, cannot update\n");
        return;
    }

    hashTable->table[index] |= outEdges;
    // printf("%zu: %hhu\n", index, hashTable->table[index]);
}

uint8_t queryHashTable(struct HashTable* hashTable, size_t key) {
    size_t hash = hashTable->hashFunction(key, hashTable->indexSize);
    size_t fingerprint = hashTable->fingerprintFunction(key);

    size_t index = hash;
    while ((hashTable->table[index] & (1 << 7)) && ((hashTable->table[index] >> 4) & 0b111) != fingerprint) { // Dirty bit is set but fingerprint doesn't match
        index = (index + 1) % (hashTable->indexSize);
        
        if (index == hash) { // Looped around the table
            // Table is full
            // printf("Item is not in table\n");
            return -1;
        }
    }

    if (!(hashTable->table[index] & (1 << 7))) { // Dirty bit is not set
        // printf("Item is not in table\n");
        return -1;
    }

    return hashTable->table[index] & 0b1111;
}