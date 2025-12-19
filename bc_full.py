import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--output')
parser.add_argument('-m', '--model')
args = parser.parse_args()

# this is the box-cox linearization script for full-landscape computation

genotypes = list()
Pobs = list()

with open(args.landscape) as f:     # parsing the landscape file (hash table format)
    for line in f:
        spl = line.split()
        genotypes.append(spl[0])
        Pobs.append(float(spl[1]))

import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import PowerTransformer

def model(Padd, lmbda, A, B):       # defining the box-cox based transform
    return ((Padd + A) ** lmbda - 1) / (lmbda * sp.stats.gmean(Padd + A) ** (lmbda - 1)) + B if lmbda !=0 else sp.stats.gmean(Padd + A)*np.log(Padd + A) + B

from collections import defaultdict

def calculate_additive_phenotypes(genotypes, phenotypes):    # calculating additive phenotypes
    n = len(genotypes)
    if n == 0:
        return []
    L = len(genotypes[0])
    wt_geno = genotypes[0]
    wt_phen = phenotypes[0]
    effects_lists = defaultdict(list)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            g1 = genotypes[i]
            g2 = genotypes[j]
            diff_pos = []
            for p in range(L):
                if g1[p] != g2[p]:
                    diff_pos.append(p)
            if len(diff_pos) == 1:
                p = diff_pos[0]
                from_aa = g1[p]
                to_aa = g2[p]
                delta = phenotypes[j] - phenotypes[i]
                key = (p, from_aa, to_aa)
                effects_lists[key].append(delta)
    effects = {}
    for key, lst in effects_lists.items():
        if lst:
            effects[key] = sum(lst) / len(lst)
    print(effects)
    additive_phenotypes = []
    for g in genotypes:
        add = wt_phen
        for p in range(L):
            if g[p] != wt_geno[p]:
                from_aa = wt_geno[p]
                to_aa = g[p]
                key = (p, from_aa, to_aa)
                add += effects[key]
        additive_phenotypes.append(add)

    return additive_phenotypes


Padd = np.array(calculate_additive_phenotypes(genotypes, Pobs))

pt = PowerTransformer()
pt.fit(Padd.reshape(-1, 1))
lambdas = pt.lambdas_        # initial guess for lambda parameter
try:
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],0,0],
        bounds=([0, -min(Padd), -np.inf], [2, np.inf, min(Pobs)]), max_nfev=1e6)     # if the initial lambda guess is in [0, 2] interval, we try to avoid very large or very little lambda values
except ValueError:
    try:
        popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],0,0],
            bounds=([-np.inf, -min(Padd), -np.inf], [np.inf, np.inf, min(Pobs)]), max_nfev=1e6) # if it is not possible, we fit with no restrictions other than mathematical
    except ValueError:
        print('FAILED TO LINEARIZE THE LANDSCAPE\n')
        sys.exit(1) # sometimes (very rarely) the landscape can not be linearized using this method
Pobs_linear = list()
for p in Pobs:
    Pobs_linear.append((popt[0] * sp.stats.gmean(Padd + popt[1]) ** (popt[0] - 1) * (p - popt[2]) + 1) ** (1 / popt[0]) - popt[1]) # applying reverse transform to observed phenotypes to get linearized values

with open(args.output, 'w') as fout: # making the output file
    for j in range(len(genotypes)):
        fout.write(genotypes[j] + '\t' + str(Pobs_linear[j]) + '\n')

x = np.linspace(min(Padd), max(Padd), 1000)
y = np.zeros(1000)
y = model(x, popt[0], popt[1], popt[2])

with open(args.model, 'w') as fmod: # making the model file for plotting
    print('x\ty', file=fmod)
    for xi, yi in zip(x, y):
        print(f'{xi},{yi}')
