import argparse

parser = argparse.ArgumentParser(description='This is a test')
parser.add_argument('-k', '--kmer-length', required=True, type=int, default=31, help='kmer length')
parser.add_argument('-r', '--reads-file', required=True, type=str, default='../../../src/log/reads.log', help='reads file')
parser.add_argument('-a', '--all-kmers-file', required=True, type=str, default='../../../src/log/all-kmers.log', help='all kmers file')

args = parser.parse_args()
K = args.kmer_length

print("Reading reads file...")
with open(args.reads_file, 'r') as reads_file:
    initial_kmers = set()
    all_kmers = set()
    for i, line in enumerate(reads_file):
        print("Line {}".format(i), end='\r')
        initial_kmers.add(line[0:K])
        for i in range(len(line) - K):
            all_kmers.add(line[i:i + K])

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

print("== Results ==============================")
print("Total possible k-mers: {}".format(4**K))
print("Total k-mers in reads: {}".format(len(all_kmers)))
print("Total visited k-mers: {}".format(len(visited_kmers)))
positive_results = set(filter(lambda kmer: query_results[kmer] != 255, visited_kmers))
print("Total found k-mers: {}".format(len(positive_results)))
false_negatives = all_kmers - visited_kmers
print("Total false negatives: {}".format(len(false_negatives)))
if false_negatives:
    print("False negatives: {}".format(false_negatives))
false_positives = positive_results - all_kmers
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