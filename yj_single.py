'''
This is the Yeo-Johnston linearization script for a single combinatorially complete landscape
'''

import numpy as np
import scipy as sp
import argparse
import matplotlib.pyplot as plt
from sklearn.preprocessing import PowerTransformer
from scipy.stats import gmean
from collections import defaultdict

def model(P_add: np.ndarray[float], lam: float, A: float, B: floa):
    '''
    Define the Yeo-Johnson based transform

    Parameters
    ----------
    Padd : np.ndarray[float]
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
    
    # --- Step 1: Apply the shift parameter 'A' ---
    try:
        shifted_P = P_add + A
    except ValueError:
        return np.full_like(P_add, np.inf)

    # --- Step 2: Critical Constraint Check (for GM) ---
    # The geometric mean (GM) in your denominator REQUIRES positive inputs.
    # If 'A' makes any value non-positive, this parameter set is invalid.
    if np.any(shifted_P <= 0):
        return np.full_like(P_add, np.inf)

    # --- Step 3: Calculate the Geometric Mean (GM) ---
    gm = gmean(shifted_P)

    # --- Step 4: Calculate the Core Yeo-Johnson Transform (Psi) ---
    # We only need the case for positive inputs (Y >= 0)
    # where Y = shifted_P.
    
    if np.isclose(lam, 0):
        # Case for lam = 0
        psi = np.log(shifted_P + 1)
    else:
        # Case for lam != 0
        psi = (np.power(shifted_P + 1, lam) - 1) / lam

    # --- Step 5: Calculate the Normalization Factor ---
    # This is the (GM**(lam - 1)) term
    
    # Handle the lam = 1 edge case (GM**0 = 1)
    if np.isclose(lam, 1.0):
        norm_factor = 1.0
    else:
        norm_factor = np.power(gm, lam - 1)
    
    # Avoid division by zero if norm_factor is somehow zero
    if np.isclose(norm_factor, 0):
        return np.full_like(P_add, np.inf)

    # --- Step 6: Apply scaling and final shift 'B' ---
    P_obs = (psi / norm_factor) + B
    
    return P_obs

def model_inverse(P_obs: np.ndarray[float], lmbda: float, A: float, B: float, GM: float):
    '''
    Define the inverse Yeo-Johnson based transform
    
    Parameters
    ----------
    P_obs : np.ndarray[float]
        OBSERVED PHENOTYPES.
    lmbda : float
        LMBDA PARAMETER.
    A : float
        A PARAMETER.
    B : float
        B PARAMETER.
    GM : float
        GEOMETRIC MEAN OF ADDITIVE PHENOTYPES + A.

    Returns
    -------
    np.ndarray[float]
        LINEARIZED PHENOTYPES.

    '''
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

def calculate_additive_phenotypes(genotypes: list[str], phenotypes: list[float]):    # Calculate additive phenotypes
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
print(lambdas)

try:
    # If the initial lambda guess is in [0, 2] interval, we try to avoid very large or very little lambda values
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],-min(Padd)+5*abs(max(Padd)),0],
        bounds=([0, -min(Padd), -np.inf], [2, np.inf, np.inf]), max_nfev=1e6)  
except ValueError:
    # If it is not possible, we fit with no restrictions other than mathematical
    popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],-min(Padd)+5*abs(max(Padd)),0],
        bounds=([-np.inf, -min(Padd), -np.inf], [np.inf, np.inf, np.inf]), max_nfev=1e6) 
print(popt)
print(pcov)

# Apply inverse transform to observed phenotypes to get linearized values
Pobs_linear = list()
for p in Pobs:
    Pobs_linear.append(model_inverse(p,popt[0],popt[1],popt[2],sp.stats.gmean(Padd+popt[1]))) 

# Write the output to file
with open(args.filename[:-4] + '_linearized_47.csv', 'w') as f:
    for j in range(len(genotypes)):
        f.write(genotypes[j] + '\t' + str(Pobs_linear[j]) + '\n')

# Make the model file for plotting the power transform function
x = np.linspace(min(Padd), max(Padd), 1000)
y = np.zeros(1000)

y = model(x, popt[0], popt[1], popt[2])

fig, ax = plt.subplots()
ax.set_xlabel('Observed phenotype')
ax.set_ylabel('Additive / Linearized phenotype')
plt.plot(Padd, Pobs, 'ko')
plt.plot(x, y, 'r')
plt.savefig(args.filename[:-4] + '_linearized_47.png')
