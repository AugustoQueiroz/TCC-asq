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
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <k> <sketch_file> <starting_kmers_file> <output_file> <presence_threshold>" << std::endl;
        return 1;
    }
    size_t K = std::stoi(argv[1]);
    std::string sketch_fp = argv[2];
    std::ifstream startingKMersFile(argv[3]);
    std::ofstream outputFile(argv[4]);
    size_t presence_threshold = std::stoi(argv[5]);

    // Load the sketch from the file
    std::cout << "Loading sketch from " << sketch_fp << std::endl;
    struct DeBruijnCountMin* sketch = loadDeBruijnCountMin(sketch_fp.c_str());

    // Setup for graph traversal
    std::cout << "Setting up for traversal" << std::endl;
    std::queue<size_t> toVisit;
    std::set<size_t> visited;

    // Load the starting kmers
    std::cout << "Loading starting kmers" << std::endl;
    std::string startingKMer;
    while (startingKMersFile >> startingKMer) {
        toVisit.push(mapKMer(startingKMer.c_str()));
    }
    startingKMersFile.close();

    // Traversing the graph
    std::cout << "Traversing graph" << std::endl;
    size_t currentKMer = 0;
    while (!toVisit.empty()) {
        currentKMer = toVisit.front();
        toVisit.pop();
        visited.insert(currentKMer);

        uint16_t queryResult = queryDeBruijnCountMin(sketch, currentKMer);
        outputFile << kMerFromCode(currentKMer, K) << ": " << (queryResult & COUNTER_MASK) << " " << std::bitset<4>(queryResult >> 12) << std::endl;

        if (queryResult >= presence_threshold) {
            uint8_t outEdges = queryResult >> 12;
            if (outEdges & 0b1000) {
                size_t nextKMer = extendKMer(currentKMer, 'A', K);
				size_t reverseComplementKMer = mapKMer(reverseComplement(kMerFromCode(nextKMer, K)));
                if (visited.find(nextKMer) == visited.end())
                    toVisit.push(nextKMer);
				if (visited.find(reverseComplementKMer) == visited.end())
					toVisit.push(reverseComplementKMer);
            }
            if (outEdges & 0b0100) {
                size_t nextKMer = extendKMer(currentKMer, 'C', K);
				size_t reverseComplementKMer = mapKMer(reverseComplement(kMerFromCode(nextKMer, K)));
                if (visited.find(nextKMer) == visited.end())
                    toVisit.push(nextKMer);
				if (visited.find(reverseComplementKMer) == visited.end())
					toVisit.push(reverseComplementKMer);
            }
            if (outEdges & 0b0010) {
                size_t nextKMer = extendKMer(currentKMer, 'G', K);
				size_t reverseComplementKMer = mapKMer(reverseComplement(kMerFromCode(nextKMer, K)));
                if (visited.find(nextKMer) == visited.end())
                    toVisit.push(nextKMer);
				if (visited.find(reverseComplementKMer) == visited.end())
					toVisit.push(reverseComplementKMer);
            }
            if (outEdges & 0b0001) {
                size_t nextKMer = extendKMer(currentKMer, 'T', K);
				size_t reverseComplementKMer = mapKMer(reverseComplement(kMerFromCode(nextKMer, K)));
                if (visited.find(nextKMer) == visited.end())
                    toVisit.push(nextKMer);
				if (visited.find(reverseComplementKMer) == visited.end())
					toVisit.push(reverseComplementKMer);
            }
        }
    }
    outputFile.close();

    return 0;
}
