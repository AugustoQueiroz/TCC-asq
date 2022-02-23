import argparse
import os
import re
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Read analysis')
parser.add_argument('-r', '--results-folder', type=str, help='Folder containing the results file for each value of alpha + one file containing all the k-mers in the read that was used (named "kmers_in_read.log"). The results files should be named <alpha>.log (e.g.: 0.5.log), and the format inside it should be, one per line, <k-mer>: <k-mer code> - <hash result> - <fingerprint> - <query result>.')

def main(results_folder):
    files = os.listdir(results_folder)
    result_files = list(filter(lambda filename: re.compile("^[0-9]+(\.[0-9]+)?\.log$").match(filename), files))
    read_file = "kmers_in_read.log"

    kmers_in_read = set()
    with open(results_folder + '/' + read_file) as read_file:
        for line in read_file:
            kmer = line.strip()
            kmers_in_read.add(kmer)
    
    results = {}
    for result_filename in result_files:
        with open(results_folder + '/' + result_filename) as result_file:
            alpha = float(result_filename[:-4])
            results[alpha] = {}

            for line in result_file:
                kmer, result = line.strip().split(':')
                kmer_code, hash_result, fingerprint, query_result = result.strip().split(' - ')
                results[alpha][kmer] = {
                    'kmer_code': int(kmer_code),
                    'hash_result': int(hash_result),
                    'fingerprint': int(fingerprint),
                    'query_result': int(query_result),
                    'in_read': kmer in kmers_in_read
                }

    false_positive_counts = {}
    for alpha in sorted(results):
        false_positive_counts[alpha] = 0
        for kmer in results[alpha]:
            if results[alpha][kmer]['in_read'] and results[alpha][kmer]['query_result'] != 255:
                false_positive_counts[alpha] += 1
        
        print('False positive count (alpha = {}): {}'.format(alpha, false_positive_counts[alpha]))

    x = sorted(false_positive_counts)
    y = list(map(lambda alpha: false_positive_counts[alpha], x))
    plt.plot(x, y)
    plt.xlabel('Alpha')
    plt.ylabel('False positive count')
    plt.title('False positive count vs. Alpha\n(|R| = 1024, k = 8)')
    plt.xlim(right=max(x))
    plt.show()

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.results_folder)