import numpy as np
import pwlf
import matplotlib.pyplot as plt
import math
import sklearn
import argparse
from pulp import *

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

genotypes = list()
phenotypes = list()

with open(args.filename) as f:
    for line in f:
        genotypes.append(line.split(',')[0])
        phenotypes.append(float(line.split(',')[1].strip()))

A = [(0,)]
nummut = int(math.log2(len(phenotypes)))

for i in range(1, nummut + 1):
    A.append((i,))
    for elem in A[1:-1]:
        A.append(elem+(i,))

Z = list(zip(phenotypes, A))
A_sorted = [i[1] for i in sorted(Z, key=lambda x: x[0])]

var = np.concatenate([np.array([0]), np.array(list(map(LpVariable, np.full(nummut, 'a', dtype=str) + np.array(list(map(str, np.arange(1, nummut+1)))))))], dtype=object)

prob = LpProblem('FP', LpMinimize)

for i in range(len(A_sorted)-1):
    prob += sum([var[j] for j in A_sorted[i]]) <= sum([var[k] for k in A_sorted[i+1]]) - 1

prob.solve()

coefs = list()

for i in range(1, len(var)):
    coefs.append(value(var[i]))

print(prob)
print(A_sorted)
print(coefs)

wildtype = genotypes[0]
fps = list()

for genotype in genotypes:
    fp = 0
    for i in range(len(genotype)):
        if genotype[i] != wildtype[i]:
            fp += coefs[i]
    fps.append(fp)

print(fps)

pwlf_ = pwlf.PiecewiseLinFit(fps, phenotypes)
pwlf_.fit_with_breaks(fps)
x = np.linspace(min(fps), max(fps), 10000)
y = pwlf_.predict(x)

plt.figure()
plt.plot(fps, phenotypes, 'ko')
plt.plot(x, y, 'g-')
plt.savefig(args.filename[:-4]+'_potential.png')
