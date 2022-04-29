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
parser.add_argument('-nexe', '--navigation-executable', type=str, help='The path to the navigation executable to be used.')

# Experiments to Run
experiments_group = parser.add_argument_group('Experiments to run')
parser.add_argument('-fp', '--false-positive', action='store_true', help='Run the false positive experiment.')
parser.add_argument('-fn', '--false-negative', action='store_true', help='Run the false negative experiment.')

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

class Experiment:
    valid_experiments = ['false_positive', 'false_negative']

    def __init__(self, D, W, k, presence_threshold, experiments_to_run, dataset_handler, path_to_exe, path_to_navigation_exe):
        self.D = D
        self.W = W
        self.k = k
        self.presence_threshold = presence_threshold
        self.experiments_to_run = experiments_to_run
        self.dataset_handler = dataset_handler
        self.results = {}

        self.path_to_exe = path_to_exe
        self.path_to_navigation_exe = path_to_navigation_exe

    def run(self):
        run_command = f'{self.path_to_exe} {self.k} {self.dataset_handler.read_length} {self.W} {self.D} {self.presence_threshold} < {self.dataset_handler.reads_file}'
        print(f'Executing: {run_command}')
        os.system(run_command)
        mv_command = f'mv sketch.bin K{self.k}W{self.W}D{self.D}T{self.presence_threshold}.sketch'
        print(f'Executing: {mv_command}')
        os.system(mv_command)
        mv_command = f'mv starting-kmers.txt K{self.k}W{self.W}D{self.D}T{self.presence_threshold}.starting'
        print(f'Executing: {mv_command}')
        os.system(mv_command)

        self.run_navigation()

    def run_navigation(self):
        run_command = f'{self.path_to_navigation_exe} {self.k} K{self.k}W{self.W}D{self.D}T{self.presence_threshold}.sketch K{self.k}W{self.W}D{self.D}T{self.presence_threshold}.starting K{self.k}W{self.W}D{self.D}T{self.presence_threshold}.results {self.presence_threshold}'
        print(f'Executing: {run_command}')
        os.system(run_command)

    def something(self):
        actual_kmer_counts = self.dataset_handler.get_k_mer_counts(self.k)
        navigatable_kmers = self.results['sketch'].navigatable_k_mers(set(read[:self.k] for read in self.dataset_handler.get_reads()), self.presence_threshold)
        present_kmers = set(filter(lambda kmer: navigatable_kmers[kmer][0] >= self.presence_threshold, navigatable_kmers))

        self.results['actual_kmer_counts'] = actual_kmer_counts
        self.results['navigatable_kmers'] = navigatable_kmers

        print(f'Total possible distinct k-mers: {4**self.k}')
        print(f'Actual distinct k-mer count: {len(actual_kmer_counts)}')
        print(f'Navigatable distinct k-mer count: {len(navigatable_kmers)}')
        if len(actual_kmer_counts) > len(navigatable_kmers):
            print(f'Unnavigated k-mers from reads: {set(actual_kmer_counts) - set(navigatable_kmers)}')
        print(f'Present distinct k-mer count: {len(present_kmers)}')

        false_positive_kmers = self.get_false_positives()
        print(f'False Positives: {len(false_positive_kmers)}')
        print(f'\tRate (fp / navigated): {100*(len(false_positive_kmers) / len(present_kmers))}%')

        false_negative_kmers = self.get_false_negatives()
        print(f'False Negatives: {len(false_negative_kmers)}')
        print(f'\tRate (fn / actual count): {100*(len(false_negative_kmers) / len(actual_kmer_counts))}%')
    
    def get_false_positives(self):
        kmers_found = set(filter(lambda kmer: self.results['navigatable_kmers'][kmer][0] >= self.presence_threshold, self.results['navigatable_kmers']))
        kmers_present = set(self.results['actual_kmer_counts'])

        return kmers_found - kmers_present
    
    def get_false_negatives(self):
        kmers_found = set(filter(lambda kmer: self.results['navigatable_kmers'][kmer][0] >= self.presence_threshold, self.results['navigatable_kmers']))
        kmers_present = set(self.results['actual_kmer_counts'])

        return kmers_present - kmers_found

    def save_results(self, output_dir):
        with open(os.path.join(output_dir, f'K{self.k}W{self.W}D{self.D}T{self.presence_threshold}'), 'w') as results_file:
            results_file.write(f'{self.results["navigatable_kmers"]}\n')
            results_file.write(f'{len(self.get_false_positives())}\n')
            results_file.write(f'{len(self.get_false_negatives())}\n')
    
class ExperimentSet:
    def __init__(self, args, dataset_handler):
        self.args = args
        self.experiments = []

        experiments_to_run = list(map(lambda experiment: vars(args)[experiment], Experiment.valid_experiments))

        variations = product(args.D, args.W, args.k, args.presence_threshold)
        for variation in variations:
            experiment = Experiment(*variation, experiments_to_run, dataset_handler, self.args.executable, self.args.navigation_executable)
            self.experiments.append(experiment)
    
    def run(self):
        with mp.Pool(self.args.threads) as pool:
            pool.map(Experiment.run, self.experiments)

    def save_results(self):
        for experiment in self.experiments:
            experiment.save_results(self.args.output_dir)

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
    # experiment_set.save_results()