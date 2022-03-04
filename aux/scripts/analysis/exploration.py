import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Explore the results of from the representation model. Count total explored k-mers, the number of possible k-mers that weren\'t explored, the number of false positives, and the number of edges each k-mer had in the graph.')
parser.add_argument('-a', '--all-kmers-file', type=str, required=False, default=None, help='File containing all possible k-mers and their query results on the representation structure.')
parser.add_argument('-n', '--navigated-kmers-file', type=str, required=True, help='File containing all k-mers that were found by navigating the graph, and their query results on the representation structure.')
parser.add_argument('-r', '--reads-file', type=str, required=True, help='File containing all the reads.')
parser.add_argument('-k', '--kmer-length', type=int, required=True, help='Length of the k-mers.')

def get_data(all_kmers_file_path, navigated_kmers_file_path, read_kmers_file_path, kmer_length):
    with open(read_kmers_file_path, 'r') as read_kmers_file:
        kmers_in_sequence = set()
        for line in read_kmers_file:
            for i in range(len(line) - kmer_length):
                kmers_in_sequence.add(line[i:i+kmer_length])

    with open(navigated_kmers_file_path, 'r') as output_file:
        navigated_kmers_query_results = {}
        for line in output_file:
            try:
                kmer, query_result = line.split(':')
            except:
                print(line)
            if kmer in navigated_kmers_query_results:
                print("Duplicate kmer found: {}".format(kmer))
            navigated_kmers_query_results[kmer] = int(query_result)

    all_kmers_query_results = None
    if all_kmers_file_path is not None:
        with open(all_kmers_file_path, 'r') as all_kmers_file:
            all_kmers_query_results = {}
            for line in all_kmers_file:
                kmer, query_result = line.split(':')
                if kmer in all_kmers_query_results:
                    print("Duplicate kmer found: {}".format(kmer))
                all_kmers_query_results[kmer] = int(query_result)
    
    return kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results

def print_data_summary(kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results, k):
    print("Total possible k-mers: {}".format(4**k))
    print("Total k-mers in sequence: {}".format(len(kmers_in_sequence)))
    print("Total visited k-mers: {}".format(len(navigated_kmers_query_results)))
    print("Total unvisited k-mers: {}".format(4**k - len(navigated_kmers_query_results)))
    print("Total unvisited k-mers (from sequence): {}".format(len(kmers_in_sequence - set(navigated_kmers_query_results))))
    if len(kmers_in_sequence - set(navigated_kmers_query_results)) > 0:
        print(kmers_in_sequence - set(navigated_kmers_query_results))

def print_false_positive_results(kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results):
    if all_kmers_query_results is not None:
        all_false_positives = set()
        for kmer, result in all_kmers_query_results.items():
            if kmer not in kmers_in_sequence and result != 255:
                all_false_positives.add(kmer)
        print('False positives (all k-mers): {}'.format(len(all_false_positives)))
    
    false_positives_from_navigation = set()
    for kmer, result in navigated_kmers_query_results.items():
        if kmer not in kmers_in_sequence and result != 255:
            false_positives_from_navigation.add(kmer)

    print('False positives (only visited in traversal): {}'.format(len(false_positives_from_navigation)))

def print_edge_count_results(kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results):
    edge_count_frequencies = {}
    for kmer, result in navigated_kmers_query_results.items():
        if result == 255:
            continue
        edge_count = 0
        for _ in range(4):
            edge_count += result & 1
            result >>= 1
        if edge_count not in edge_count_frequencies:
            edge_count_frequencies[edge_count] = 0
        edge_count_frequencies[edge_count] += 1

    for edge_count, frequency in edge_count_frequencies.items():
        print('k-mers with {} edges: {}'.format(edge_count, frequency))

def main(all_kmers_file_path, navigated_kmers_file_path, read_kmers_file_path, kmer_length):
    kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results = get_data(all_kmers_file_path, navigated_kmers_file_path, read_kmers_file_path, kmer_length)

    print_data_summary(kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results, kmer_length)
    print_false_positive_results(kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results)
    print_edge_count_results(kmers_in_sequence, navigated_kmers_query_results, all_kmers_query_results)
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args.all_kmers_file, args.navigated_kmers_file, args.reads_file, args.kmer_length)