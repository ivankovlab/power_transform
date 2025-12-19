import math
import os
import argparse

# this is the script for calculating the quantity of the 2nd order epistasis

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

coeffs = list()

with open(args.filename) as f:
    for line in f:
        if line.startswith('['):
            if ' ]' not in line:
                coeffs.append(line.split()[-1][:-1])
            else:
                coeffs.append(line.split()[-2])

print(f'Total amount of epistasis is {math.sqrt((sum((float(i) ** 2 for i in coeffs))) / len(coeffs))}')
