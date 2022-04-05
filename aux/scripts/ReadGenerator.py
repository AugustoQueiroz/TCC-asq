from math import ceil
from random import randint

class ReadGenerator:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def generate_reads(self, *, length: int=100, coverage: int=50) -> list[str]:
        number_of_reads = ceil(len(self.sequence) * coverage / length)
        for _ in range(number_of_reads):
            start = randint(0, len(self.sequence) - length)
            yield self.sequence[start:start+length]

if __name__ == '__main__':
    with open('results/exp1/ecoli.txt') as sequence_file:
        sequence = sequence_file.read().strip()
    print(len(sequence))
    generator = ReadGenerator(sequence)

    with open('results/exp1/reads.txt', 'w') as output_file:
        for read in generator.generate_reads(length=252, coverage=80):
            output_file.write(read + '\n')