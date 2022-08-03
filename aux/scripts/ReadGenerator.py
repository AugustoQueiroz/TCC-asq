from math import ceil
from random import randint, uniform

import argparse

from experiments.fastq_reader import FastQParser

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sequence-file', type=str)
parser.add_argument('-f', '--fastq', action='store_true')
parser.add_argument('-t', '--txt', action='store_true')
parser.add_argument('-l', '--read-length', type=int, default=100)
parser.add_argument('-c', '--coverage', type=int, default=80)
parser.add_argument('-rc', '--reverse-complement-probability', type=float, default=0)

class ReadGenerator:
    def __init__(self, sequence: str, reverse_complement_probability: float):
        self.sequence = sequence
        self.reverse_complement_probability = reverse_complement_probability

    def reverse_complement(self, read):
        complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        read_complement = ''.join(map(lambda base: complements[base], read))
        return read_complement[::-1]

    def generate_reads(self, *, length: int=100, coverage: int=50) -> list[str]:
        number_of_reads = ceil(len(self.sequence) * coverage / length)
        for _ in range(number_of_reads):
            start = randint(0, len(self.sequence) - length)
            if uniform(0, 1) < self.reverse_complement_probability:
                yield self.reverse_complement(self.sequence[start:start+length])
            else:
                yield self.sequence[start:start+length]

if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.sequence_file) as sequence_file:
        if args.fastq:
            sequence = ''.join(sequence_file.readlines()[1:]).replace('\n', '')
        else: 
            sequence = sequence_file.read().strip()
    print(len(sequence))
    generator = ReadGenerator(sequence, args.reverse_complement_probability)

    with open('./reads.txt', 'w') as output_file:
        for read in generator.generate_reads(length=args.read_length, coverage=args.coverage):
            output_file.write(read + '\n')