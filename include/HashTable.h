#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include <stdint.h>
#include <stdlib.h>

struct HashTable {
    size_t indexSize; // # of bits of the index
    uint8_t* table;
    size_t(*hashFunction)(size_t, size_t);
    size_t(*fingerprintFunction)(size_t);
};

void insertIntoHashTable(struct HashTable* hashTable, size_t key);
void updateHashTableWithEdges(struct HashTable* hashTable, size_t key, uint8_t outEdges);
uint8_t queryHashTable(struct HashTable* hashTable, size_t key);

#endif
