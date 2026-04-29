'''
This script is for calculating epistasis coeffients.
Use it only in one particular case: after you have used either bc_hcubewise.py or yj_hcubewise.py.
It is designed for taking the output of these scripts as input.
'''

import numpy as np
import scipy as sp
import argparse
import sys

np.printoptions(threshold=sys.maxsize, linewidth=np.inf)

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--order')
parser.add_argument('-f', '--output')
args = parser.parse_args()

land = dict()
hcubes = list()

order = int(args.order)
H = sp.linalg.hadamard(2 ** order) # as in Poelwijk (2016)

diag = np.array([1.0])
for _ in range(order):
    diag = np.concatenate([0.5 * diag, -1.0 * diag])
V = np.diag(diag) # as in Poelwijk (2016)
VH = np.matmul(V, H)

with open(args.landscape) as fin, open(args.output, 'w') as fout:
    for line in fin:
        if line != 'FAILED TO LINEARIZE THIS HYPERCUBE\n':
            P = np.array([float(part.split(',')[1].strip()) for part in line.strip().split()]) # Phenotypes in one hypercube
            K = np.matmul(VH, P)
            print(K, file=fout)
        else:
            fout.write(line) # for consistency
