import numpy as np
import scipy as sp
import argparse
import sys

#this script is for a single combinatorially complete landscape

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
    P.append(line.strip()[1].split())
  P = np.array(P)
  K = np.matmul(VH, P)
  print(K, file=fout)
