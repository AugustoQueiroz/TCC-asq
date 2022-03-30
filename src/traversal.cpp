#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <fstream>
#include <bitset>

#include <stdlib.h>

extern "C" {
#include "CountMin.h"
#include "KMerProcessing.h"
}

#define PRESENCE_THRESHOLD 25

std::string extendKMer(std::string currentKMer, char nextBase) {
  std::string nextKMer = currentKMer.substr(1);
  nextKMer.push_back(nextBase);
  return nextKMer;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <k> <sketch_file> <starting_kmers_file> <output_file>" << std::endl;
        return 1;
    }
    size_t K = std::stoi(argv[1]);
    std::string sketch_fp = argv[2];
    std::ifstream startingKMersFile(argv[3]);
    std::ofstream outputFile(argv[4]);

    // Load the sketch from the file
    std::cout << "Loading sketch from " << sketch_fp << std::endl;
    FILE* sketch_file = fopen(sketch_fp.c_str(), "r");
    struct DeBruijnCountMin* sketch = loadDeBruijnCountMin(sketch_file);
    fclose(sketch_file);

    // Setup for graph traversal
    std::cout << "Setting up for traversal" << std::endl;
    std::queue<std::string> toVisit;
    std::set<std::string> visited;

    // Load the starting kmers
    std::cout << "Loading starting kmers" << std::endl;
    std::string startingKMer;
    while (startingKMersFile >> startingKMer) {
        toVisit.push(startingKMer);
    }
    startingKMersFile.close();

    // Traversing the graph
    std::cout << "Traversing graph" << std::endl;
    char* currentKMer = (char*) malloc((K+1) * sizeof(char));
    currentKMer[K] = '\0';
    while (!toVisit.empty()) {
        visited.insert(toVisit.front());
        strcpy(currentKMer, toVisit.front().c_str());
        currentKMer[K] = '\0';
        toVisit.pop();
        size_t currentKMerCode = mapKMer(currentKMer);

        uint16_t queryResult = queryDeBruijnCountMin(sketch, currentKMerCode);
        outputFile << currentKMer << ": " << (queryResult & COUNTER_MASK) << " " << std::bitset<4>(queryResult >> 12) << std::endl;

        if (queryResult >= PRESENCE_THRESHOLD) {
            uint8_t outEdges = queryResult >> 12;
            if (outEdges & 0b1000) {
                std::string nextKMer = extendKMer(currentKMer, 'A');
                if (visited.find(nextKMer) == visited.end()) {
                    toVisit.push(nextKMer);
                }
            }
            if (outEdges & 0b0100) {
                std::string nextKMer = extendKMer(currentKMer, 'C');
                if (visited.find(nextKMer) == visited.end()) {
                    toVisit.push(nextKMer);
                }
            }
            if (outEdges & 0b0010) {
                std::string nextKMer = extendKMer(currentKMer, 'G');
                if (visited.find(nextKMer) == visited.end()) {
                    toVisit.push(nextKMer);
                }
            }
            if (outEdges & 0b0001) {
                std::string nextKMer = extendKMer(currentKMer, 'T');
                if (visited.find(nextKMer) == visited.end()) {
                    toVisit.push(nextKMer);
                }
            }
        }
    }
    free(currentKMer);
    outputFile.close();

    return 0;
}