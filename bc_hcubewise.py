import argparse
import numpy as np
import scipy as sp
from sklearn.preprocessing import PowerTransformer
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--hypercubes')
parser.add_argument('-l', '--landscape')
parser.add_argument('-o', '--output')
args = parser.parse_args()

# This is the Box-Cox linearization script for hypercube-wise computation

land = dict()

with open(args.landscape) as f:  # Parse the landscape file (genotypes in sequence format)
    for line in f:
        spl = line.split()
        land[spl[0]] = spl[1]

def model(Padd: np.array[float], lmbda: float, A: float, B: float):  # Define the Box-Cox based transform
    """
    

    Parameters
    ----------
    Padd : np.array[float]
        ADDITIVE PHENOTYPES.
    lmbda : float
        LMBDA PARAMETER.
    A : float
        A PARAMETER.
    B : float
        B PARAMETER.

    Returns
    -------
    np.array[float]
        OBSERVED PHENOTYPES.

    """
    if lmbda !=0:
        return ((Padd + A) ** lmbda - 1) / (lmbda * sp.stats.gmean(Padd + A) ** (lmbda - 1)) + B
    else:
        return sp.stats.gmean(Padd + A) * np.log(Padd + A) + B

genotypes = list()
Pobs = list()

def calculate_additive_phenotypes(genotypes: list[str], phenotypes: list[float]):    # Calculate additive phenotypes
    """
    

    Parameters
    ----------
    genotypes : list[str]
        GENOTYPES.
    phenotypes : list[float]
        OBSERVED PHENOTYPES.

    Returns
    -------
    list[float]
        ADDITIVE PHENOTYPES.

    """
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


with open(args.hypercubes) as fin, open(args.output, 'w') as fout:
    for line in fin:
        genotypes = list()
        genotypes.append(line.split('\t')[1].strip())
        muts = line.split('\t')[0].split(':')
        for mut in muts:
            temp = list()
            for seq in genotypes:
                temp.append(seq[:int(mut[1:-1])]+mut[0]+seq[int(mut[1:-1])+1:])
            genotypes.extend(temp) # Define each hypercube in situ

        Pobs = list(map(float, [land[genotype] for genotype in genotypes])) # Find phenotypes for needed genotypes
        Padd = np.array(calculate_additive_phenotypes(genotypes, Pobs))

        pt = PowerTransformer()
        pt.fit(Padd.reshape(-1, 1))
        lambdas = pt.lambdas_ # Initial guess for lambda parameter
        print(lambdas)
        try:
            popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],0,0],   
                bounds=([0, -min(Padd), -np.inf], [2, np.inf, min(Pobs)]), max_nfev=1e6)  
            # If the initial lambda guess is in [0, 2] interval, we try to avoid very large or very little lambda values
        except ValueError:
            try:
                popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],0,0],
                    bounds=([-np.inf, -min(Padd), -np.inf], [np.inf, np.inf, min(Pobs)]), max_nfev=1e6)
                # If it is not possible, we fit with no restrictions other than mathematical
            except ValueError:
                print('FAILED TO LINEARIZE THE HYPERCUBE\n', file=fout)
                continue   
            # Sometimes (very rarely) the hypercube can not be linearized using this method.
            # in this case we give up with this hypercube and continue with other hypercubes
        print(popt)
        print(pcov)
        Pobs_linear = list()
        for p in Pobs:
            Pobs_linear.append((popt[0] * sp.stats.gmean(Padd + popt[1]) ** (popt[0] - 1) * (p - popt[2]) + 1) ** (1 / popt[0]) - popt[1])
            # Apply reverse transform to observed phenotypes to get linearized values

        for j in range(len(genotypes)):
            fout.write(genotypes[j] + ',' + str(Pobs_linear[j]) + '\t') # Output
        fout.write('\n')
