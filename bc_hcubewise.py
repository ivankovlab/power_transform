import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--hypercubes')
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--output')
args = parser.parse_args()

# this is the box-cox linearization script for hypercube-wise computation

land = dict()

with open(args.landscape) as f:  # parsing the landscape file (hash table format)
    for line in f:
        spl = line.split()
        land[spl[0]] = spl[1]

import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import PowerTransformer
import pwlf

def model(Padd, lmbda, A, B):   # defining the box-cox based transform
    return ((Padd + A) ** lmbda - 1) / (lmbda * sp.stats.gmean(Padd + A) ** (lmbda - 1)) + B if lmbda !=0 else sp.stats.gmean(Padd + A)*np.log(Padd + A) + B

genotypes = list()
Pobs = list()

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


with open(args.hypercubes) as fin, open(args.output, 'w') as fout:
    for line in fin:
        genotypes = list()
        genotypes.append(line.split('\t')[1].strip())
        muts = line.split('\t')[0].split(':')
        for mut in muts:
            temp = list()
            for seq in genotypes:
                temp.append(seq[:int(mut[1:-1])]+mut[0]+seq[int(mut[1:-1])+1:])
            genotypes.extend(temp) # defining each hypercube in situ

        Pobs = list(map(float, [land[genotype] for genotype in genotypes])) # finding phenotypes for needed genotypes
        Padd = np.array(calculate_additive_phenotypes(genotypes, Pobs))

        pt = PowerTransformer()
        pt.fit(Padd.reshape(-1, 1))
        lambdas = pt.lambdas_ # initial guess for lambda parameter
        print(lambdas)
        try:
            popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],0,0],   
                bounds=([0, -min(Padd), -np.inf], [2, np.inf, min(Pobs)]), max_nfev=1e6)  # if the initial lambda guess is in [0, 2] interval, we try to avoid very large or very little lambda values
        except ValueError:
            try:
                popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],0,0],
                    bounds=([-np.inf, -min(Padd), -np.inf], [np.inf, np.inf, min(Pobs)]), max_nfev=1e6) # if it is not possible, we fit with no restrictions other than mathematical
            except ValueError:
                print('FAILED TO LINEARIZE THE HYPERCUBE\n', file=fout)
                continue   # sometimes (very rarely) the hypercube can not be linearized using this method. in this case we give up with this hypercube and continue with other hypercubes
        print(popt)
        print(pcov)
        Pobs_linear = list()
        for p in Pobs:
            Pobs_linear.append((popt[0] * sp.stats.gmean(Padd + popt[1]) ** (popt[0] - 1) * (p - popt[2]) + 1) ** (1 / popt[0]) - popt[1]) # applying reverse transform to observed phenotypes to get linearized values

        for j in range(len(genotypes)):
            fout.write(genotypes[j] + ',' + str(Pobs_linear[j]) + '\t') # output
        fout.write('\n')
