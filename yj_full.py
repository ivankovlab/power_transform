import argparse
import sys

# this is the yeo-johnston linearization script for full-landscape computation

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--output')
parser.add_argument('-m', '--model')
args = parser.parse_args()

genotypes = list()
Pobs = list()

with open(args.landscape) as f:  # parsing the landscape file (hash table format)
    for line in f:
        spl = line.split()
        genotypes.append(spl[0])
        Pobs.append(float(spl[1]))

import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import PowerTransformer
import pwlf

def model(P_add, lam, A, B):   # defining the yeo-johnston based transform
    try:
        shifted_P = P_add + A
    except ValueError:
        return np.full_like(P_add, np.inf)

    if np.any(shifted_P <= 0):
        return np.full_like(P_add, np.inf)

    gm = sp.stats.gmean(shifted_P)

    if np.isclose(lam, 0):
        psi = np.log(shifted_P + 1)
    else:
        psi = (np.power(shifted_P + 1, lam) - 1) / lam

    if np.isclose(lam, 1.0):
        norm_factor = 1.0
    else:
        norm_factor = np.power(gm, lam - 1)

    if np.isclose(norm_factor, 0):
        return np.full_like(P_add, np.inf)

    P_obs = (psi / norm_factor) + B

    return P_obs

def model_inverse(P_obs, lmbda, A, B, GM):
    P_obs = np.array(P_obs, dtype=float)
    
    scaling_factor = GM ** (lmbda - 1)
    W = (P_obs - B) * scaling_factor
    
    Y_shifted = np.zeros_like(W)
    
    pos_mask = W >= 0
    neg_mask = ~pos_mask

    if np.isclose(lmbda, 0):
        Y_shifted[pos_mask] = np.exp(W[pos_mask]) - 1
    else:
        base = np.maximum(lmbda * W[pos_mask] + 1, 0)
        Y_shifted[pos_mask] = np.power(base, 1.0 / lmbda) - 1

    if np.isclose(lmbda, 2):
        Y_shifted[neg_mask] = 1 - np.exp(-W[neg_mask])
    else:
        term = (2 - lmbda)
        base = np.maximum(1 - term * W[neg_mask], 0)
        Y_shifted[neg_mask] = 1 - np.power(base, 1.0 / term)

    return Y_shifted - A

from collections import defaultdict

def calculate_additive_phenotypes(genotypes, phenotypes):  # calculating additive phenotypes
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


Padd = np.array(calculate_additive_phenotypes(genotypes, Pobs))

pt = PowerTransformer()
pt.fit(Padd.reshape(-1, 1))
lambdas = pt.lambdas_       # initial guess for lambda parameter
try:
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],-min(Padd)+5*max(Padd),0],
        bounds=([0, -min(Padd), -np.inf], [2, np.inf, np.inf]), max_nfev=1e6)   # if the initial lambda guess is in [0, 2] interval, we try to avoid very large or very little lambda values
except ValueError:
    try:
        popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],-min(Padd)+5*max(Padd),0],
            bounds=([-np.inf, -min(Padd), -np.inf], [np.inf, np.inf, np.inf]), max_nfev=1e6)  # if it is not possible, we fit with no restrictions other than mathematical
        except ValueError:
            print('FAILED TO LINEARIZE THE LANDSCAPE\n')
            sys.exit(1) # sometimes (very rarely) the landscape can not be linearized using this method
Pobs_linear = list()
for p in Pobs:
    Pobs_linear.append(model_reverse(p,popt[0],popt[1],popt[2],sp.stats.gmean(Padd+popt[1])))   # applying reverse transform to observed phenotypes to get linearized values

with open(args.output, 'w') as fout:  # making the output file
    for j in range(len(genotypes)):
        fout.write(genotypes[j] + '\t' + str(Pobs_linear[j]) + '\n')

x = np.linspace(min(Pobs), max(Pobs), 1000)
y = np.zeros(1000)
y = model(x, popt[0], popt[1], popt[2])

with open(args.model, 'w') as fmod:  # making the model file for plotting
    print('x\ty', file=fmod)
    for xi, yi in zip(x, y):
        print(f'{xi},{yi}',file=fmod)
