import math
import os
import argparse

# this is the script for calculating the quantity of the 3rd order epistasis. it also calculates 2nd order epistasis in 3d hypercubes.
# do not use it for 2d hypercubes.

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

coeffs2 = list()
coeffs3 = list()

with open(args.filename) as f:
    for line in f:
        if line.startswith(('[', 'E')):
            spl = line.split()
            spl = [i.strip('[]') for i in spl]
            while '' in spl:
                spl.remove('')
            coeffs3.append(spl[-1])
            coeffs2.extend(spl[-4:-1])

print(f'Total amount of 2nd order epistasis is {math.sqrt((sum((float(i) ** 2 for i in coeffs2))) / len(coeffs2))}')
print(f'Total amount of 3rd order epistasis is {math.sqrt((sum((float(i) ** 2 for i in coeffs3))) / len(coeffs3))}')
