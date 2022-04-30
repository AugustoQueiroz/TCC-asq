#include <map>
#include <fstream>
#include <string>
#include <iostream>

extern "C" {
#include "CountMin.h"
#include "KMerProcessing.h"
}

int main(int argc, char** argv) {
    if (argc != 5) {
        printf("usage: %s <K> <read file> <sketch file> <output file>", argv[0]);
        return 1;
    }
    std::ifstream readsFile(argv[2]);
    int K = std::stoi(argv[1]);

    struct DeBruijnCountMin* sketch = loadDeBruijnCountMin(argv[3]);

    std::map<std::string, int> kmerCounts = std::map<std::string, int>();
    std::string read;

    while (readsFile >> read) {
        for (size_t i = 0; i < read.length() - K + 1; i++) {
            std::string kmer = read.substr(i, K);
            kmerCounts[kmer]++;
        }
    }

    std::ofstream outputFile(argv[4]);
    for (auto counters_it = kmerCounts.begin(); counters_it != kmerCounts.end(); ++counters_it) {
        size_t kmerCode = mapKMer(counters_it->first.c_str());
        uint16_t sketchCount = queryDeBruijnCountMin(sketch, kmerCode) & COUNTER_MASK;
        outputFile << counters_it->first << "," << counters_it->second << "," << sketchCount << std::endl;
    }

    return 0;
}