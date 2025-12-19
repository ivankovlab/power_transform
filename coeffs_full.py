import numpy as np
import scipy as sp
import argparse
import sys

# this script is for full-landscape computations only. do not use it for hypercube-wise.

np.printoptions(threshold=sys.maxsize, linewidth=np.inf)

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--hypercubes')
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--order')
parser.add_argument('-f', '--output')
args = parser.parse_args()

land = dict()
hcubes = list()

with open(args.hypercubes) as f: # Parsing HypercubeME-2 output (hash table)
    for line in f:
        hcubes.append(list())
        hcubes[-1].append(line.split('\t')[1].strip())
        muts = line.split('\t')[0].split(':')
        for mut in muts:
            temp = list()
            for seq in hcubes[-1]:
                temp.append(seq[:int(mut[1:-1])]+mut[0]+seq[int(mut[1:-1])+1:])
            hcubes[-1].extend(temp)

with open(args.landscape) as f: # Parsing fitness landscape (hash table)
    for line in f:
        spl = line.split('\t')
        land[spl[0]] = spl[1]

order = int(args.order)
H = sp.linalg.hadamard(2 ** order)  # as in Poelwijk (2016)

diag = np.array([1.0])
for _ in range(order):
    diag = np.concatenate([0.5 * diag, -1.0 * diag])
V = np.diag(diag)  # as in Poelwijk (2016)
VH = np.matmul(V, H)

P = np.zeros(2 ** order)
i = 0

with open(args.output, 'w') as f:
    for hcube in hcubes:
        for i, genotype in enumerate(hcube):
            P[i] = land[genotype]
        K = np.matmul(VH, P)
        print(K, file=f)
