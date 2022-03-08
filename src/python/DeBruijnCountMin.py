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
    counter_mask = (1 << 60) - 1
    large_prime = 18446744073709551557

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
        out_edges = reduce(lambda x, y: x & y, [self.table[i][hashes[i]] & ~self.counter_mask for i in range(self.D)], ~0) >> 60
        return count, out_edges

    def count_frequencies(self):
        count_frequencies = {}
        for i in range(1):
            for j in range(self.W):
                count = self.table[i][j] & self.counter_mask
                if count in count_frequencies:
                    count_frequencies[count] += 1
                else:
                    count_frequencies[count] = 1
        print(sorted(count_frequencies))
        print(sum(count_frequencies.keys()))
        plt.bar(count_frequencies.keys(), count_frequencies.values())
        plt.show()

if __name__ == '__main__':
    sketch = DeBruijnCountMin.from_file('./sketch.bin')
    sketch.count_frequencies()