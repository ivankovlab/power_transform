import numpy as np
import scipy as sp
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('order')
args = parser.parse_args()

H = sp.linalg.hadamard(2 ** int(args.order))
X = sp.linalg.inv(H)

for filename in os.listdir('data' + args.order):
    f = os.path.join(args.directory, filename)
    P = np.zeros(2 ** args.order)
    i = 0
    with open(f) as f_:
        for line in f_:
            P[i] = int(line.split('\t')[1].strip())
    K = np.matmul(X, P)
    with open(f + '_coeffs.txt', 'w') as fw:
        fw.write(K)
