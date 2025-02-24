import numpy as np
import scipy as sp
import math
import argparse
import matplotlib.pyplot as plt
from sklearn.preprocessing import PowerTransformer
import pwlf

def model(Padd, lmbda, A, B):
    return ((Padd + A) ** lmbda - 1) / (lmbda * sp.stats.gmean(Padd + A) ** (lmbda - 1)) + B if lmbda !=0 else sp.stats.gmean(Padd + A)*np.log(Padd + A) + B

def two_pass_variance(data):
    n = len(data)
    mean = sum(data) / n
    variance = sum((x - mean) ** 2 for x in data) / (n - 1)
    return variance

def two_pass_covariance(data1, data2):
    n = len(data1)
    mean1 = sum(data1) / n
    mean2 = sum(data2) / n
    covariance = 0
    for i1, i2 in zip(data1, data2):
        a = i1 - mean1
        b = i2 - mean2
        covariance += a * b / n
    return covariance

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

genotypes = list()
Pobs = list()

with open(args.filename) as f:
    for line in f:
        genotypes.append(line.split(',')[0])
        phenotype = float(line.split(',')[1].strip())
        Pobs.append(phenotype)

Padd = np.zeros(len(Pobs))
Padd[0] = Pobs[0]
for i in range(1, len(Pobs)):
    if math.ceil(math.log2(i)) == math.floor(math.log2(i)):
        indices = [j for j in range(i, len(Pobs)) if j % (2 * i) >= i]
        temp = list()
        for j in indices:
            temp.append(Pobs[j] - Pobs[j-i])
        Padd[i] = np.mean(temp) + Padd[0]
    else:
        indices = [2 ** index for index, bit in enumerate(bin(i)[:1:-1]) if bit == '1']
        Padd[i] = np.sum(Padd[indices] - Padd[0]) + Padd[0]
print(Padd)

pt = PowerTransformer()
pt.fit(Padd.reshape(-1, 1))
lambdas = pt.lambdas_
print(lambdas)

popt, pcov = sp.optimize.curve_fit(f=model, xdata=Padd, ydata=Pobs, sigma=0.01, p0=[lambdas[0],0,0],
        bounds=([0, -min(Padd), -np.inf], [2, np.inf, min(Pobs)]), max_nfev=1e6)
print(popt)
print(pcov)

Pobs_linear = list()

for p in Pobs:
    Pobs_linear.append((popt[0] * sp.stats.gmean(Padd + popt[1]) ** (popt[0] - 1) * (p - popt[2]) + 1) ** (1 / popt[0]) - popt[1])

with open(args.filename[:-4] + '_linearized.csv', 'w') as f:
    for j in range(len(genotypes)):
        f.write(genotypes[j] + ',' + str(Pobs_linear[j]) + '\n')

Padd_ = list(Padd)
x = np.linspace(min(Padd_+Pobs), max(Padd_+Pobs), 1000)
y = np.zeros(1000)

y = model(x, popt[0], popt[1], popt[2])

plt.plot(Padd, Pobs, 'ko')
plt.plot(x, y, 'r')
plt.savefig(args.filename[:-4] + '_linearized.png')

pwlf_ = pwlf.PiecewiseLinFit(Padd, Pobs)
pwlf_.fit_with_breaks(Padd)
y_pwlf = pwlf_.predict(x)

cov_pwlf = two_pass_covariance(x, y_pwlf)
cov_model = two_pass_covariance(x, y)
cov_default = two_pass_covariance(Padd_, Pobs)
print('PWLF ', cov_pwlf, ' MY MODEL ', cov_model, ' DEFAULT ', cov_default)

var_default = two_pass_variance(Pobs)
var_linear = two_pass_variance(Pobs_linear)
print('DEFAULT VARIANCE ', var_default, ' LINEARIZED VARIANCE ', var_linear)
