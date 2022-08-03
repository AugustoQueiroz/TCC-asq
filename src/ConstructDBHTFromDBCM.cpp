#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <fstream>
#include <bitset>

#include <chrono>

#include <stdlib.h>

extern "C" {
#include "CountMin.h"
#include "HashTable.h"
#include "KMerProcessing.h"
}

std::set<size_t> loadStartingKMers(const char* startingKMersFilePath) {
    std::ifstream startingKMersFile(startingKMersFilePath);
    std::set<size_t> startingKMers = std::set<size_t>();

    std::string startingKMer;
    while (startingKMersFile >> startingKMer) {
        startingKMers.insert(mapKMer(startingKMer.c_str()));
    }

    startingKMersFile.close();
    return startingKMers;
}

// TODO - Move into KMerProcessing.c
size_t extendKMer(size_t currentKMer, char nextBase, int K) {
    size_t kmerMask = ((size_t) 1 << (2*K)) - 1;
    size_t nextKMer = (currentKMer << 2);
    switch (nextBase) {
        case 'A':
        nextKMer |= 0b00;
        break;
        case 'C':
        nextKMer |= 0b01;
        break;
        case 'G':
        nextKMer |= 0b10;
        break;
        case 'T':
        nextKMer |= 0b11;
        break;
    }
    return nextKMer & kmerMask;
}

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " <k> <hash_table_size> <sketch_file> <starting_kmers_file> <dbcm_output_file> <dbht_output_file> <presence_threshold>" << std::endl;
        return 1;
    }
    size_t K = std::stoi(argv[1]);
    size_t tableSize = std::stoi(argv[2]);
    std::string sketch_fp = argv[3];
    std::ifstream startingKMersFile(argv[4]);
    std::ofstream dBCMOutputFile(argv[5]);
    std::ofstream dBHTOutputFile(argv[6]);
    size_t presence_threshold = std::stoi(argv[7]);

    // Load the sketch from the file
    std::cout << "Loading sketch from " << sketch_fp << std::endl;
    struct DeBruijnCountMin* sketch = loadDeBruijnCountMin(sketch_fp.c_str());

    // Setup for graph traversal
    std::cout << "Setting up for traversal" << std::endl;
    std::queue<size_t> toVisit;
    std::set<size_t> visited;

    // Load the starting kmers
    std::set<size_t> startingKMers = loadStartingKMers(argv[4]);
    std::cout << "Loading starting kmers" << std::endl;
    for (auto starting_it = startingKMers.begin(); starting_it != startingKMers.end(); ++starting_it) {
        toVisit.push(*starting_it);
    }

    std::cout << tableSize << std::endl;
    struct HashTable* dBHT = createHashTable(tableSize);

    // Traversing the graph
    std::cout << "Traversing DBCM" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    size_t currentKMer = 0;
    while (!toVisit.empty()) {
        currentKMer = toVisit.front();
        toVisit.pop();
        visited.insert(currentKMer);

        uint16_t queryResult = queryDeBruijnCountMin(sketch, currentKMer);
        dBCMOutputFile << kMerFromCode(currentKMer, K) << "," << (queryResult & COUNTER_MASK) << "," << std::bitset<4>(queryResult >> 12) << std::endl;

        if (queryResult >= presence_threshold) {
            insertIntoHashTable(dBHT, currentKMer);
            uint8_t outEdges = queryResult >> 12;
            if (outEdges & 0b1000) {
                size_t nextKMer = extendKMer(currentKMer, 'A', K);
                if (isMemberOfDeBruijnCountMin(sketch, nextKMer, presence_threshold)) {
                    updateHashTableWithEdges(dBHT, currentKMer, 0b1000);
                    if (queryHashTable(dBHT, nextKMer))
                        toVisit.push(nextKMer);
                }
            }
            if (outEdges & 0b0100) {
                size_t nextKMer = extendKMer(currentKMer, 'C', K);
                if (isMemberOfDeBruijnCountMin(sketch, nextKMer, presence_threshold)) {
                    updateHashTableWithEdges(dBHT, currentKMer, 0b0100);
                    if (queryHashTable(dBHT, nextKMer))
                        toVisit.push(nextKMer);
                }
            }
            if (outEdges & 0b0010) {
                size_t nextKMer = extendKMer(currentKMer, 'G', K);
                if (isMemberOfDeBruijnCountMin(sketch, nextKMer, presence_threshold)) {
                    updateHashTableWithEdges(dBHT, currentKMer, 0b0010);
                    if (queryHashTable(dBHT, nextKMer))
                        toVisit.push(nextKMer);
                }
            }
            if (outEdges & 0b0001) {
                size_t nextKMer = extendKMer(currentKMer, 'T', K);
                if (isMemberOfDeBruijnCountMin(sketch, nextKMer, presence_threshold)) {
                    updateHashTableWithEdges(dBHT, currentKMer, 0b0001);
                    if (queryHashTable(dBHT, nextKMer))
                        toVisit.push(nextKMer);
                }
            }
        }
    }
    dBCMOutputFile.close();
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to construct DBHT through traversal of DBCM (in microseconds): " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;

    // Load the starting kmers
    std::cout << "Loading starting kmers" << std::endl;
    for (auto starting_it = startingKMers.begin(); starting_it != startingKMers.end(); ++starting_it) {
        toVisit.push(*starting_it);
    }

    visited = std::set<size_t>();

    // Traversing the graph
    std::cout << "Traversing DBHT" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    currentKMer = 0;
    while (!toVisit.empty()) {
        currentKMer = toVisit.front();
        toVisit.pop();
        visited.insert(currentKMer);

        uint8_t queryResult = queryHashTable(dBHT, currentKMer);
        dBHTOutputFile << kMerFromCode(currentKMer, K) << "," << (int) queryResult << std::endl;

        if (queryResult != (uint8_t) -1) {
            if (queryResult & 0b1000) {
                size_t nextKMer = extendKMer(currentKMer, 'A', K);
                if (queryHashTable(dBHT, nextKMer) != -1) {
                    if (visited.find(nextKMer) == visited.end())
                        toVisit.push(nextKMer);
                }
            }
            if (queryResult & 0b0100) {
                size_t nextKMer = extendKMer(currentKMer, 'C', K);
                if (queryHashTable(dBHT, nextKMer) != -1) {
                    if (visited.find(nextKMer) == visited.end())
                        toVisit.push(nextKMer);
                }
            }
            if (queryResult & 0b0010) {
                size_t nextKMer = extendKMer(currentKMer, 'G', K);
                if (queryHashTable(dBHT, nextKMer) != -1) {
                    if (visited.find(nextKMer) == visited.end())
                        toVisit.push(nextKMer);
                }
            }
            if (queryResult & 0b0001) {
                size_t nextKMer = extendKMer(currentKMer, 'T', K);
                if (queryHashTable(dBHT, nextKMer) != -1) {
                    if (visited.find(nextKMer) == visited.end())
                        toVisit.push(nextKMer);
                }
            }
        }
    }
    dBHTOutputFile.close();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to traverse DBHT (in microseconds): " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;

    deleteDeBruijnCountMinSketch(sketch);
    deleteHashTable(dBHT);

    return 0;
}
