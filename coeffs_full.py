import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--hypercubes')
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

# ---------------- READ LANDSCAPE ----------------
land = dict()
with open(args.landscape) as f:
    for line in f:
        spl = line.split('\t')
        land[spl[0]] = float(spl[1])

# ---------------- PREPARE ----------------
order = args.order
n = 2 ** order

# диагональ V
diag = np.array([1.0])
for _ in range(order):
    diag = np.concatenate([diag, -diag])
diag *= 2.0 ** (-order)

# один буфер на всё время
P = np.zeros(n, dtype=np.float64)

# ---------------- STREAM PROCESSING ----------------
with open(args.hypercubes) as fin, open(args.output, 'w') as fout:
    for line in fin:
        # --- строим один гиперкуб ---
        cube = []
        base = line.split('\t')[1].strip()
        cube.append(base)

        muts = line.split('\t')[0].split(':')

        for mut in muts:
            pos = int(mut[1:-1])
            new_seqs = []
            for seq in cube:
                new_seq = seq[:pos] + mut[0] + seq[pos+1:]
                new_seqs.append(new_seq)
            cube.extend(new_seqs)

        # --- заполняем P ---
        for i, genotype in enumerate(cube):
            P[i] = land[genotype]

        # --- FWHT ---
        K = fast_hadamard_transform(P.copy())

        # --- применяем V ---
        K *= diag

        print(K, file=fout)
