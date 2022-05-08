from enum import Enum

import argparse

parser = argparse.ArgumentParser(description='Reads a fastq file and returns a list of tuples with the sequence and the quality score.')
parser.add_argument('fastq', help='The fastq file to read.')

class FastQParser:
    class ParserStep(Enum):
        ID = 1
        SEQUENCE = 2
        QUALITY_SCORES = 3

    def __init__(self):
        self.current_id = ""
        self.current_sequence = ""
        self.current_quality_scores = ""
        self.current_step = self.ParserStep.ID
        self.sequences = []

    def build_block(self):
        block = {
            'id': self.current_id,
            'sequence': self.current_sequence,
            'quality_scores': self.current_quality_scores
        }
        self.sequences.append(block)
        self.current_id = ""
        self.current_sequence = ""
        self.current_quality_scores = ""

    def parse_line(self, line):
        if self.current_step == self.ParserStep.ID:
            self.current_id = line.strip()
            self.current_step = self.ParserStep.SEQUENCE
        elif self.current_step == self.ParserStep.SEQUENCE:
            if line.startswith('+'):
                self.current_step = self.ParserStep.QUALITY_SCORES
            else:
                self.current_sequence += line.strip()
        elif self.current_step == self.ParserStep.QUALITY_SCORES:
            self.current_quality_scores += line.strip()
            if len(self.current_quality_scores) == len(self.current_sequence):
                self.current_step = self.ParserStep.ID
        
    def parse(self, fastq_file_path):
        self.sequences = []
        with open(fastq_file_path) as fastq_file:
            for i, line in enumerate(fastq_file):
                print(i, end='\r')
                self.parse_line(line)

        if self.current_step != self.ParserStep.ID:
            self.build_block()

        return self.sequences

    def process_block(self, block):
        ks = ['name', 'sequence', 'optional', 'quality_scores']
        return {k: v for k, v in zip(ks, block)}

    def naive_parse(self, fastq_file_path):
        self.sequences = []
        with open(fastq_file_path) as fastq_file:
            block = []
            for line in fastq_file:
                block.append(line)
                if len(block) == 4:
                    self.sequences.append(self.process_block(block))
                    block = []
        return self.sequences

if __name__ == '__main__':
    args = parser.parse_args()
    reader = FastQParser()
    records = reader.parse(args.fastq)
    print(records)