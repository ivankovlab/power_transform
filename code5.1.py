import numpy as np
import scipy as sp
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('order')
args = parser.parse_args()

H = sp.linalg.hadamard(2 ** int(args.order))
X = sp.linalg.inv(H)
P = np.zeros(2 ** int(args.order))
i = 0
   
with open(args.filename) as f:
    for line in f:
        P[i] = float(line.split(',')[1].strip())
        K = np.matmul(X, P)
        i += 1

with open(args.filename[:-4] + '_coeffs.txt', 'w') as fw:
    fw.write(str(K))
