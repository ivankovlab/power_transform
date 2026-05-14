'''
This is the Box-Cox linearization script for a single combinatorially complete landscape
'''

import numpy as np
import scipy as sp
from scipy.optimize import least_squares
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns

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
parser.add_argument('filename')
args = parser.parse_args()

genotypes = list()
Pobs = list()

# Parse the landscape file (genotypes in sequence format)
with open(args.filename) as f:
    for line in f:
        genotypes.append(line.split()[0])
        phenotype = float(line.split()[1].strip())
        Pobs.append(phenotype)

# Calculate additive phenotypes
Padd = np.array(calculate_additive_phenotypes(genotypes, Pobs))

with open(args.filename[:-4]+'_additive.tsv', 'w') as f:
    for gt, pt in zip(genotypes, Padd):
        print(gt+'\t'+str(pt), file=f)

# Fit the model
lmbda_opt, A_opt, B_opt = fit_model(Padd, Pobs, lambda0=1.0, A0=-np.min(Padd)+1, B0=0.0)

print(lmbda_opt, A_opt, B_opt)

# Apply reverse transform to observed phenotypes to get linearized values
Pobs_linear = list()
for p in Pobs:
    Pobs_linear.append((lmbda_opt * sp.stats.gmean(Padd + A_opt) ** (lmbda_opt - 1) * (p - B_opt) + 1) ** (1 / lmbda_opt) - A_opt)

# Write the output to file
with open(args.filename[:-4] + '_linearized.tsv', 'w') as f: # Make the output file
    for j in range(len(genotypes)):
        f.write(genotypes[j] + '\t' + str(Pobs_linear[j]) + '\n')

# Plot the observed, additive and linearized phenotypes
x = np.linspace(min(Padd), max(Padd), 1000)
y = np.zeros(1000)

y = model(x, lmbda_opt, A_opt, B_opt)

print(model(Padd, lmbda_opt, A_opt, B_opt))

labels = ['wt', 'single mutant', 'single mutant', 'double mutant']
# labels = ['wt', 'single mutant', 'single mutant', 'single mutant', 'double mutant', 'double mutant', 'double mutant', 'triple mutant']

fig, ax = plt.subplots()
ax.set_xlabel('Additive / Linearized phenotype')
ax.set_ylabel('Observed phenotype')
plt.plot(Padd, Pobs, 'ko')
plt.plot(x, y, 'r')

for i, label in enumerate(labels):
    plt.annotate(
            label,
            (Padd[i], Pobs[i]),
            textcoords='offset points',
            xytext=(0, 10),
            ha='center',
            fontsize=9
            )

plt.savefig(args.filename[:-4] + '_linearized.png')'''
This is the Box-Cox linearization script for a single combinatorially complete landscape
'''

import numpy as np
import scipy as sp
import argparse
import matplotlib.pyplot as plt
from sklearn.preprocessing import PowerTransformer
from collections import defaultdict

def model(Padd: np.ndarray, lmbda: float, A: float, B: float):
    '''
    Define the Box-Cox based transform

    Parameters
    ----------
    Padd : np.ndarray
        ADDITIVE PHENOTYPES.
    lmbda : float
        LMBDA PARAMETER.
    A : float
        A PARAMETER.
    B : float
        B PARAMETER.

    Returns
    -------
    np.ndarray[float]
        OBSERVED PHENOTYPES.

    '''
    if lmbda !=0:
        return ((Padd + A) ** lmbda - 1) / (lmbda * sp.stats.gmean(Padd + A) ** (lmbda - 1)) + B
    else:
        return sp.stats.gmean(Padd + A) * np.log(Padd + A) + B

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

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

genotypes = list()
Pobs = list()

# Parse the landscape file (genotypes in sequence format)
with open(args.filename) as f:
    for line in f:
        genotypes.append(line.split()[0])
        phenotype = float(line.split()[1].strip())
        Pobs.append(phenotype)

# Calculate additive phenotypes
Padd = np.array(calculate_additive_phenotypes(genotypes, Pobs))

# Initial guess for lambda parameter
pt = PowerTransformer()
pt.fit(Padd.reshape(-1, 1))
lambdas = pt.lambdas_

lambdas[0] = 4.0

print(lambdas)
try:
     # If the initial lambda guess is in [0, 2] interval, we try to avoid very large or very little lambda values
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=np.full(len(Pobs),0.01), p0=[lambdas[0],0,0],
        bounds=([0, -min(Padd), -np.inf], [2, np.inf, min(Pobs)]), max_nfev=1e6)
except ValueError:
    # If it is not possible, we fit with no restrictions other than mathematical
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=np.full(len(Pobs),0.01), p0=[lambdas[0],-min(Padd), min(Pobs)],
        bounds=([-np.inf, -min(Padd), -np.inf], [np.inf, np.inf, min(Pobs)]), max_nfev=1e6)
print(popt)
print(pcov)

# Apply reverse transform to observed phenotypes to get linearized values
Pobs_linear = list()
for p in Pobs:
    Pobs_linear.append((popt[0] * sp.stats.gmean(Padd + popt[1]) ** (popt[0] - 1) * (p - popt[2]) + 1) ** (1 / popt[0]) - popt[1])

# Write the output to file
with open(args.filename[:-4] + '_linearized.csv', 'w') as f: # Make the output file
    for j in range(len(genotypes)):
        f.write(genotypes[j] + '\t' + str(Pobs_linear[j]) + '\n')

# Plot the observed, additive and linearized phenotypes
x = np.linspace(min(Padd), max(Padd), 1000)
y = np.zeros(1000)

y = model(x, popt[0], popt[1], popt[2])

print(Padd)

print(model(Padd, popt[0], popt[1], popt[2]))

fig, ax = plt.subplots()
ax.set_xlabel('Additive / Linearized phenotype')
ax.set_ylabel('Observed phenotype')
plt.plot(Padd, Pobs, 'ko')
plt.plot(x, y, 'r')
plt.savefig(args.filename[:-4] + '_linearized.png')
