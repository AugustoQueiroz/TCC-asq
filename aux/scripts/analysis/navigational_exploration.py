import argparse
from functools import reduce

parser = argparse.ArgumentParser(description='This is a test')
parser.add_argument('-k', '--kmer-length', required=True, type=int, default=31, help='kmer length')
parser.add_argument('-r', '--reads-file', required=True, type=str, default='../../../src/log/reads.log', help='reads file')
parser.add_argument('-a', '--all-kmers-file', required=True, type=str, default='../../../src/log/all-kmers.log', help='all kmers file')
parser.add_argument('--navigate_without_edges', action='store_true', help='navigate without edges')

args = parser.parse_args()
K = args.kmer_length

print("Reading reads file...")
with open(args.reads_file, 'r') as reads_file:
    initial_kmers = set()
    all_kmers = {}
    for i, line in enumerate(reads_file):
        print("Line {}".format(i), end='\r')
        initial_kmers.add(line[0:K])
        for i in range(len(line) - K):
            if line[i:i + K] not in all_kmers: all_kmers[line[i:i + K]] = 0
            all_kmers[line[i:i + K]] += 1

print("Reading all kmers file...")
with open(args.all_kmers_file, 'r') as all_kmer_results_file:
    query_results = {}
    for i, line in enumerate(all_kmer_results_file):
        print("Line {}".format(i), end='\r')
        kmer, query_result = line.split(':')
        query_results[kmer] = int(query_result)

print("Navigating on the graph...")
kmers_to_visit = set(initial_kmers)
visited_kmers = set()
while kmers_to_visit:
    current_kmer = kmers_to_visit.pop()
    visited_kmers.add(current_kmer)

    # print("Visiting: '{}', Query result: {}".format(current_kmer, query_results[current_kmer]))

    if query_results[current_kmer] == 255:
        continue

    if not args.navigate_without_edges:
        if query_results[current_kmer] & 0b1000:
            next_kmer = current_kmer[1:] + 'A'
            if next_kmer not in visited_kmers:
                kmers_to_visit.add(next_kmer)
        if query_results[current_kmer] & 0b0100:
            next_kmer = current_kmer[1:] + 'C'
            if next_kmer not in visited_kmers:
                kmers_to_visit.add(next_kmer)
        if query_results[current_kmer] & 0b0010:
            next_kmer = current_kmer[1:] + 'G'
            if next_kmer not in visited_kmers:
                kmers_to_visit.add(next_kmer)
        if query_results[current_kmer] & 0b0001:
            next_kmer = current_kmer[1:] + 'T'
            if next_kmer not in visited_kmers:
                kmers_to_visit.add(next_kmer)
    else:
        for base in ['A', 'C', 'G', 'T']:
            next_kmer = current_kmer[1:] + base
            if next_kmer not in visited_kmers:
                kmers_to_visit.add(next_kmer)

print("== Results ==============================")
print("Total possible k-mers: {}".format(4**K))
print("Total k-mers in reads: {}".format(len(all_kmers)))
print("Total visited k-mers: {}".format(len(visited_kmers)))
positive_results = set(filter(lambda kmer: query_results[kmer] != 255, visited_kmers))
print("Total found k-mers: {}".format(len(positive_results)))
false_negatives = set(all_kmers) - visited_kmers
print("Total false negatives: {}".format(len(false_negatives)))
if false_negatives:
    print("False negatives: {}".format(false_negatives))
    max_count = reduce(lambda max_count, curr: max(max_count, all_kmers[curr]), false_negatives, 0)
    print("All false negatives appeared at most {} times in reads".format(max_count))
false_positives = positive_results - set(all_kmers)
print("Total false positives: {}".format(len(false_positives)))

edge_count_frequencies = {}
for kmer in query_results:
    if query_results[kmer] != 255:
        edge_count = 0
        for _ in range(4):
            if query_results[kmer] & 1:
                edge_count += 1
            query_results[kmer] >>= 1
        if edge_count not in edge_count_frequencies:
            edge_count_frequencies[edge_count] = 0
        edge_count_frequencies[edge_count] += 1
print(edge_count_frequencies)