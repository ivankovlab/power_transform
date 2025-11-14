import numpy as np
import scipy as sp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--hypercubes')
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--order')
parser.add_argument('-f', '--output')
args = parser.parse_args()

land = dict()
hcubes = list()

with open(args.hypercubes) as f:
    for line in f:
        hcubes.append(list())
        hcubes[-1].append(line.split('\t')[1].strip())
        muts = line.split('\t')[0].split(':')
        for mut in muts:
            temp = list()
            for seq in hcubes[-1]:
                temp.append(seq[:int(mut[1:-1])]+mut[0]+seq[int(mut[1:-1])+1:])
            hcubes[-1].extend(temp)

with open(args.landscape) as f:
    for line in f:
        spl = line.split('\t')
        land[spl[0]] = spl[1]

order = int(args.order)
H = sp.linalg.hadamard(2 ** order)
X = sp.linalg.inv(H)

P = np.zeros(2 ** order)
i = 0

f = open(args.output, 'a')

for hcube in hcubes:
    for i, genotype in enumerate(hcube):
        P[i] = land[genotype]
        K = np.matmul(X, P)
        f.write(K)

f.close()
