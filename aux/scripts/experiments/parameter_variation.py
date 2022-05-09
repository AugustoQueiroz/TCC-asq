from itertools import product

import os
import argparse
import multiprocessing as mp

from DeBruijnCountMin import DeBruijnCountMin
from fastq_reader import FastQParser

parser = argparse.ArgumentParser(description='Runs a parameter variation experiment.')

# Dataset Parameters
dataset_group = parser.add_mutually_exclusive_group(required=True)
dataset_group.add_argument('-r', '--reads-file', type=str, help='The reads file.')
dataset_group.add_argument('-f', '--fasta-file', type=str, help='The fasta file.')
dataset_group.add_argument('-s', '--synthetic-params', type=int, nargs=3, help='The parameters for the generation of a synthetic dataset.', metavar=('SEQUENCE_LENGTH', 'READ_LENGTH', 'COVERAGE'))
parser.add_argument('-sg', '--synthetic-dataset-generator', type=str, help='The path to the synthetic dataset generator. This must be used if the -s/--synthetic-params option is used. The generator must be an executable run as `./generator <SEQUENCE_LENGTH> <READ_LENGTH> <COVERAGE>.`')

# Experiment Parameters
parser.add_argument('-D', type=int, nargs='+', help='The different values of D to be used', default=[8])
w_group = parser.add_mutually_exclusive_group(required=True)
w_group.add_argument('-W', type=int, nargs='+', help='The different values of W to be used')
w_group.add_argument('-a', '--alpha', type=float, nargs='+', help='The different values of alpha to be used')
parser.add_argument('-k', type=int, nargs='+', help='The different values of k to be used')
parser.add_argument('-t', '--presence-threshold', type=int, nargs='+', help='The presence threshold to be used', default=[1])

# ..
parser.add_argument('-exe', '--executable', type=str, help='The path to the executable to be used.', required=True)
parser.add_argument('-nexe', '--navigation-executable', type=str, default=None, help='The path to the navigation executable to be used.')
parser.add_argument('-cexe', '--counting-executable', type=str, default=None, help='The path to the executable that perform k-mer counting.')

# Experiments to Run
parser.add_argument('--steps', type=int, nargs='+', default=[1, 2, 3], help='1 = dBCM construction, 2 = kmer counting, 3 = dBHT construction and traversal')

# Extra
parser.add_argument('--threads', type=int, help='The number of threads to be used.', default=mp.cpu_count())

# Output Parameters
parser.add_argument('-o', '--output-dir', type=str, help='The output directory.', default='./')

class DataSetHandler:
    def __init__(self, args):
        if args.reads_file:
            self.reads_file = args.reads_file
            with open(self.reads_file) as reads_file:
                self.read_length = len(reads_file.readline().strip())
        elif args.fasta_file:
            self.fasta_file = args.fasta_file
            parser = FastQParser()
            with open('reads.txt', 'w') as reads_file:
                for read in parser.naive_parse(self.fasta_file):
                    reads_file.write(read['sequence'].strip() + '\n')
                    self.read_length = len(read['sequence'])
            self.reads_file = 'reads.txt'
        else:
            self.sequence_length = args.synthetic_params[0]
            self.read_length = args.synthetic_params[1]
            self.coverage = args.synthetic_params[2]

            os.system(f'{args.synthetic_dataset_generator} {self.sequence_length} {self.read_length} {self.coverage}')
            self.reads_file = 'reads.txt'

    def get_reads(self):
        if self.reads_file:
            with open(self.reads_file) as reads_file:
                return list(map(lambda line: line.strip(), reads_file.readlines()))

    def get_k_mer_counts(self, k):
        reads = self.get_reads()
        k_mer_counts = {}
        for read in reads:
            for i in range(len(read) - k + 1):
                k_mer = read[i:i+k]
                if k_mer not in k_mer_counts:
                    k_mer_counts[k_mer] = 0
                k_mer_counts[k_mer] += 1
        return k_mer_counts
    
    def get_k_mers(self, k):
        return set(self.get_k_mer_counts(k).keys())

class ExperimentExecutables:
    def __init__(self, construction_executable: str, counting_executable: str = None, traversal_executable: str = None):
        self.construction_executable = construction_executable
        self.counting_executable = counting_executable
        self.traversal_executable = traversal_executable

class Experiment:
    valid_experiments = ['false_positive', 'false_negative']

    def __init__(self, D, W, k, presence_threshold, steps, dataset_handler, executables: ExperimentExecutables):
        self.D = D
        self.W = W
        self.k = k
        self.presence_threshold = presence_threshold
        self.experiment_name = f'K{k}W{W}D{D}T{presence_threshold}'
        self.steps = steps
        self.dataset_handler = dataset_handler
        self.results = {}

        self.executables = executables

    def run(self):
        if 1 in self.steps:
            self.run_construction()
        if self.executables.counting_executable and 2 in self.steps:
            self.run_counting()
        if self.executables.traversal_executable and 3 in self.steps:
            self.run_traversal()
    
    def run_construction(self):
        run_command = f'{self.executables.construction_executable} {self.k} {self.dataset_handler.read_length} {self.W} {self.D} {self.presence_threshold} < {self.dataset_handler.reads_file}'
        print(f'Executing: {run_command}')
        os.system(run_command)
        mv_command = f'mv sketch.bin {self.experiment_name}.sketch'
        print(f'Executing: {mv_command}')
        os.system(mv_command)
        mv_command = f'mv starting-kmers.txt {self.experiment_name}.starting'
        print(f'Executing: {mv_command}')
        os.system(mv_command)

    def run_counting(self):
        run_command = f'{self.executables.counting_executable} {self.k} {self.dataset_handler.reads_file} {self.experiment_name}.sketch {self.experiment_name}.counts'
        print(f'Executing: {run_command}')
        os.system(run_command)

    def run_traversal(self):
        run_command = f'{self.executables.traversal_executable} {self.k} {self.experiment_name}.sketch {self.experiment_name}.starting {self.experiment_name}.dbcm.results {self.experiment_name}.dbht.results {self.presence_threshold}'
        print(f'Executing: {run_command}')
        os.system(run_command)
    
class ExperimentSet:
    def __init__(self, args, dataset_handler):
        self.args = args
        self.experiments = []

        executables = ExperimentExecutables(self.args.executable, self.args.counting_executable, self.args.navigation_executable)

        variations = product(args.D, args.W, args.k, args.presence_threshold)
        for variation in variations:
            experiment = Experiment(*variation, self.args.steps, dataset_handler, executables)
            self.experiments.append(experiment)
    
    def run(self):
        with mp.Pool(self.args.threads) as pool:
            pool.map(Experiment.run, self.experiments)

if __name__ == '__main__':
    args = parser.parse_args()
    
    # Make sure all the paths are absolute
    if args.reads_file:
        args.reads_file = os.path.abspath(args.reads_file)
    if args.fasta_file:
        args.fasta_file = os.path.abspath(args.fasta_file)
    if args.synthetic_dataset_generator:
        args.synthetic_dataset_generator = os.path.abspath(args.synthetic_dataset_generator)
    args.executable = os.path.abspath(args.executable)
    if args.counting_executable:
        args.counting_executable = os.path.abspath(args.counting_executable)
    if args.navigation_executable:
        args.navigation_executable = os.path.abspath(args.navigation_executable)
    args.output_dir = os.path.abspath(args.output_dir)

    #
    if args.threads != 1:
        raise NotImplementedError("Multithreaded execution is not implemented yet")

    # Change working directory to output directory
    os.chdir(args.output_dir)

    dataset_handler = DataSetHandler(args)

    experiment_set = ExperimentSet(args, dataset_handler)
    experiment_set.run()