#!/usr/bin/env python3

""" This will take in a bam and will check barcodes that have been added by
    umi_tools and will check them against a whitelist barcode and correct
    those barcodes that are not on the list
"""

import argparse
import math
import itertools
import concurrent.futures
import pygtrie

QUAL_OFFSET = 33

# Global shared across workers
whitelist_trie = None
bc_probabilities = None
_max_edit_dist = None
_min_post_prob = None


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", required=True, help="Input TSV")
    parser.add_argument("-o", "--outfile", required=True, help="Output TSV")
    parser.add_argument("-w", "--whitelist", required=True, help="Whitelist barcodes")
    parser.add_argument("-b", "--barcode_count", required=True, help="Barcode count file")
    parser.add_argument("--max_edit_dist", type=int, default=2, help="Max edit distance")
    parser.add_argument("--min_post_prob", type=float, default=0.975, help="Min posterior prob")
    parser.add_argument("--skip_header", action="store_false", dest="print_header", help="Skip header")

    return parser.parse_args()


def read_whitelist(path):
    trie = pygtrie.CharTrie()
    with open(path) as f:
        for line in f:
            trie[line.strip()] = True
    return trie


def calculate_bc_ratios(path):
    ratios = {}
    total = 0
    with open(path) as f:
        for line in f:
            bc, count = line.strip().split(",")
            count = int(count)
            ratios[bc] = count + 1  # avoid zero
            total += count
    return ratios


def get_edit_probability(q_score):
    return math.pow(10, -1 * (ord(q_score) - QUAL_OFFSET) / 10)


def get_mismatch_locs(bc1, bc2):
    return [i for i, (a, b) in enumerate(zip(bc1, bc2)) if a != b]


def get_bc_probability(query_bc, query_qual, potential_bc):
    mismatch_idxs = get_mismatch_locs(query_bc, potential_bc)
    prob = 1
    for idx in mismatch_idxs:
        prob *= get_edit_probability(query_qual[idx])
    prob *= bc_probabilities.get(potential_bc, 0)
    return prob


def get_mutated_bcs(bc, max_dist):
    for dist in range(1, max_dist + 1):
        for positions in itertools.combinations(range(len(bc)), dist):
            bc_list = [[base] for base in bc]
            for pos in positions:
                orig = bc_list[pos][0]
                bc_list[pos] = [b for b in "ACGT" if b != orig]
            for variant in itertools.product(*bc_list):
                yield "".join(variant)


def get_similar_bcs(query_bc, query_qual):
    candidates = []
    for alt_bc in get_mutated_bcs(query_bc, _max_edit_dist):
        if whitelist_trie.has_key(alt_bc):
            prob = get_bc_probability(query_bc, query_qual, alt_bc)
            if prob > 0:
                candidates.append((alt_bc, prob))
    return candidates


def find_correct_barcode(candidates):
    total_prob = sum(p for _, p in candidates)
    best_bc = ""
    best_prob = 0
    for bc, prob in candidates:
        likelihood = prob / total_prob
        if likelihood > best_prob and likelihood > _min_post_prob:
            best_bc = bc
            best_prob = likelihood
    return best_bc


def correct_line(line):
    line = line.strip()
    if not line or line.startswith("read_id"):
        return None
    try:
        read_id, bc, bc_qual, *_ = line.split('\t')
    except ValueError:
        return None
    if not bc:
        return None
    if whitelist_trie.has_key(bc):
        return line + '\t' + bc
    candidates = get_similar_bcs(bc, bc_qual)
    if candidates:
        corrected = find_correct_barcode(candidates)
        if corrected:
            return line + '\t' + corrected
    return None


def main():
    global whitelist_trie, bc_probabilities, _max_edit_dist, _min_post_prob

    args = parse_args()

    whitelist_trie = read_whitelist(args.whitelist)
    bc_probabilities = calculate_bc_ratios(args.barcode_count)
    _max_edit_dist = args.max_edit_dist
    _min_post_prob = args.min_post_prob

    with open(args.infile) as f:
        lines = f.readlines()

    header = lines[0].strip()
    data_lines = lines[1:]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(correct_line, data_lines))

    with open(args.outfile, 'w') as out:
        if args.print_header:
            out.write(header + '\tcorrected_bc\n')
        for r in results:
            if r:
                out.write(r + '\n')


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
