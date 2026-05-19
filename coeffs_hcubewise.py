import numpy as np
import argparse
import sys

np.set_printoptions(threshold=sys.maxsize, linewidth=np.inf)

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--order', type=int)
parser.add_argument('-f', '--output')
args = parser.parse_args()

# ---------------- FAST HADAMARD ----------------
def fast_hadamard_transform(x):
    n = x.shape[0]
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                a = x[j]
                b = x[j + h]
                x[j] = a + b
                x[j + h] = a - b
        h *= 2
    return x

# ---------------- PREPARE ----------------
order = args.order
n = 2 ** order

# diagonal V
diag = np.array([1.0])
for _ in range(order):
    diag = np.concatenate([diag, -diag])
diag *= 2.0 ** (-order)

# ---------------- MAIN ----------------
with open(args.landscape) as fin, open(args.output, 'w') as fout:
    for line in fin:
        # parse a hypercube
        P = np.array(
            [float(part.split(',')[1].strip()) for part in line.strip().split()],
            dtype=np.float64
        )

        # FWHT
        K = fast_hadamard_transform(P.copy())

        # apply V
        K *= diag

        print(K, file=fout)
