from functools import reduce
from random import randint
from matplotlib import pyplot as plt

def map_kmer(kmer):
    kmer_code = 0
    for base in kmer:
        kmer_code <<= 2
        if base == 'A':
            kmer_code += 0
        elif base == 'C':
            kmer_code += 1
        elif base == 'G':
            kmer_code += 2
        elif base == 'T':
            kmer_code += 3
    return kmer_code

class DeBruijnCountMin:
    counter_mask = ((1 << 60) - 1) & 0xFFFFFFFFFFFFFFFF
    out_edges_mask = ~counter_mask & 0xFFFFFFFFFFFFFFFF
    large_prime = 4294967311

    @staticmethod
    def from_file(filename):
        with open(filename, 'rb') as f:
            W = int.from_bytes(f.read(8), byteorder='little')
            D = int.from_bytes(f.read(8), byteorder='little')
            sketch = DeBruijnCountMin(W, D)
            hashCoefficients = []
            for i in range(D):
                hashCoefficients.append((int.from_bytes(f.read(8), byteorder='little'), int.from_bytes(f.read(8), byteorder='little')))
            sketch.hash_coefficients = hashCoefficients
            counters = []
            for i in range(D):
                row = []
                for j in range(W):
                    row.append(int.from_bytes(f.read(8), byteorder='little'))
                counters.append(row)

            sketch.table = counters
        return sketch

    def __init__(self, W, D):
        self.W = W
        self.D = D
        self.hash_coefficients = [(randint(0, W), randint(0, W)) for _ in range(D)]
        self.table = [[0 for _ in range(W)] for _ in range(D)]

    def hash(self, x, i):
        return ((x * self.hash_coefficients[i][0]) + self.hash_coefficients[i][1]) % self.large_prime % self.W
    
    def query(self, x):
        hashes = [self.hash(x, i) for i in range(self.D)]
        count = min(self.table[i][hashes[i]] & self.counter_mask for i in range(self.D))
        out_edges = reduce(lambda x, y: x & y, [self.table[i][hashes[i]] & self.out_edges_mask for i in range(self.D)], 15 << 60) >> 60
        return count, out_edges

    def get_neighbors(self, kmer):
        kmer_code = map_kmer(kmer)
        count, out_edges = self.query(kmer_code)
        base_value = {
            'A': 0,
            'C': 1,
            'G': 2,
            'T': 3
        }
        if count > 25:
            neighbors = set(kmer[1:] + base for base in base_value if out_edges & (0b1000 >> base_value[base]))
            return neighbors
        print(f'{kmer} not considered present ({count = })')
        return set()

    def navigate_from(self, kmer, k):
        KMER_MASK = (1 << (2 * k)) - 1

        visited_kmers = set()
        frontier = {kmer}
        unitig = kmer[:-1]
        while len(frontier) > 0:
            current_kmer = frontier.pop()
            unitig += current_kmer[-1]
            
            visited_kmers.add(current_kmer)
            neighbors = self.get_neighbors(current_kmer)
            if len(neighbors) == 1:
                frontier.update(neighbors)
            else:
                return (unitig, neighbors)

    def __repr__(self):
        return f"""
        {self.W = }, {self.D = }
        {self.hash_coefficients = }
        """

    def count_frequencies(self):
        # count_frequencies = {}
        # for i in range(1):
        #     for j in range(self.W):
        #         count = self.table[i][j] & self.counter_mask
        #         if count in count_frequencies:
        #             count_frequencies[count] += 1
        #         else:
        #             count_frequencies[count] = 1
        # print(sorted(count_frequencies))
        # print(sum(count_frequencies.keys()))
        plt.bar([i for i in range(self.W)], list(map(lambda x: x & self.counter_mask, self.table[0])))
        plt.show()
        plt.hist(list(map(lambda x: x & self.counter_mask, self.table[0])), bins=250)
        plt.show()
        # for i in range(0, self.W, 10):
        #     plt.bar([i for i in range(i, i+10)], list(map(lambda x: x & self.counter_mask, self.table[0][i:i+10])))
        #     plt.show()
        # for i in range(self.D):
        #     for j in range(self.W):
        #         print(self.table[i][j], end=' ')
        #     print()

if __name__ == '__main__':
    sketch = DeBruijnCountMin.from_file('./sketch.bin')
    print(sketch)
    # sketch.count_frequencies()
    to_navigate_from = { "CATGAATT" }
    navigated = set()
    unitigs = set()
    while len(to_navigate_from):
        source = to_navigate_from.pop()
        navigated.add(source)
        unitig, neighbors = sketch.navigate_from(source, 8)
        unitigs.add(unitig)
        # print(unitig, neighbors)
        to_navigate_from.update(neighbors - navigated)
    print(*unitigs, sep='\n')