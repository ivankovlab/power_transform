import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--output')
parser.add_argument('-m', '--model')
args = parser.parse_args()

genotypes = list()
Pobs = list()

with open(args.landscape) as f:
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

def model(P_add, lam, A, B):
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

from collections import defaultdict

def calculate_additive_phenotypes(genotypes, phenotypes):
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
lambdas = pt.lambdas_
try:
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],-min(Padd)+5*max(Padd),0],
        bounds=([0, -min(Padd), -np.inf], [2, np.inf, np.inf]), max_nfev=1e6)
except ValueError:
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],-min(Padd)+5*max(Padd),0],
        bounds=([-np.inf, -min(Padd), -np.inf], [np.inf, np.inf, np.inf]), max_nfev=1e6)
Pobs_linear = list()
for p in Pobs:
    Pobs_linear.append((popt[0] * sp.stats.gmean(Padd + popt[1]) ** (popt[0] - 1) * (p - popt[2]) + 1) ** (1 / popt[0]) - popt[1])

with open(args.output, 'w') as fout:
    for j in range(len(genotypes)):
        fout.write(genotypes[j] + ',' + str(Pobs_linear[j]) + '\n')

x = np.linspace(min(Padd), max(Padd), 1000)
y = np.zeros(1000)
y = model(x, popt[0], popt[1], popt[2])

with open(args.model, 'w') as fmod:
    print('x,y', file=fmod)
    for xi, yi in zip(x, y):
        print(f'{xi},{yi}',file=fmod)
