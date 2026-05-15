'''
This is the Box-Cox linearization script for full-landscape computation
'''

import numpy as np
import scipy as sp
from scipy.optimize import least_squares
import argparse
from collections import defaultdict

def model(Padd, lmbda, A, B):
    """
    Модель зависимости Pobs от Padd.
    """
    x = Padd + A

    # защита от log/степени отрицательных чисел
    if np.any(x <= 0):
        return np.full_like(Padd, np.nan)

    gmean = sp.stats.gmean(x)

    # случай lambda ~= 0
    if np.abs(lmbda) < 1e-8:
        return gmean * np.log(x) + B

    return ((x ** lmbda - 1) /
            (lmbda * gmean ** (lmbda - 1))) + B


def residuals(params, Padd, Pobs):
    """
    Вектор остатков для оптимизации.
    """
    lmbda, A, B = params

    pred = model(Padd, lmbda, A, B)

    # штраф за невалидные значения
    if np.any(np.isnan(pred)) or np.any(np.isinf(pred)):
        return np.ones_like(Pobs) * 1e10

    return pred - Pobs


def fit_model(Padd, Pobs,
              lambda0=1.0,
              A0=0.1,
              B0=0.0):
    """
    Подбор оптимальных lambda, A, B.
    """

    Padd = np.asarray(Padd, dtype=float)
    Pobs = np.asarray(Pobs, dtype=float)

    # начальные параметры
    x0 = np.array([lambda0, A0, B0])

    # ограничение:
    # Padd + A должно быть > 0
    Amin = -np.min(Padd) + 1e-9

    bounds_lower = [-10, Amin, -np.inf]
    bounds_upper = [10, np.inf, np.inf]

    result = least_squares(
        residuals,
        x0=x0,
        bounds=(bounds_lower, bounds_upper),
        args=(Padd, Pobs)
    )

    lmbda_opt, A_opt, B_opt = result.x

    return lmbda_opt, A_opt, B_opt

def calculate_additive_phenotypes(genotypes: list[str], phenotypes: list[float]):
    '''
    Calculate additive phenotypes

    Parameters
    ----------
    genotypes : list[str]
        GENOTYPES IN SEQUENCE FORMAT.
    phenotypes : list[float]
        OBSERVED PHENOTYPES.

    Returns
    -------
    list[float]
        ADDITIVE PHENOTYPES.

    '''
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

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--output')
parser.add_argument('-m', '--model')
args = parser.parse_args()

genotypes = list()
Pobs = list()
land = dict()

with open(args.landscape) as f:     # Parse the landscape file (genotypes in sequence format)
    for line in f:
        spl = line.split()
        genotypes.append(spl[0])
        Pobs.append(float(spl[1]))
        
Padd = np.array(calculate_additive_phenotypes(genotypes, Pobs))

# Fit the model
lmbda_opt, A_opt, B_opt = fit_model(Padd, Pobs, lambda0=1.0, A0=-np.min(Padd)+1, B0=0.0)
print(lmbda_opt, A_opt, B_opt)

# Apply reverse transform to observed phenotypes to get linearized values
Pobs_linear = list()
for p in Pobs:
    Pobs_linear.append((lmbda_opt * sp.stats.gmean(Padd + A_opt) ** (lmbda_opt - 1) * (p - B_opt) + 1) ** (1 / lmbda_opt) - A_opt)

# Write the output to file
with open(args.output, 'w') as fout:
    for j in range(len(genotypes)):
        fout.write(genotypes[j] + '\t' + str(Pobs_linear[j]) + '\n')

x = np.linspace(min(Padd), max(Padd), 1000)
y = np.zeros(1000)
y = model(x, lmbda_opt, A_opt, B_opt)

with open(args.model, 'w') as fmod: # Make the model file for plotting
    print('x\ty', file=fmod)
    for xi, yi in zip(x, y):
        print(f'{xi},{yi}')
