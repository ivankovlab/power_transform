#!/bin/python3

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--hypercubes')
parser.add_argument('-l', '--landscape')
args = parser.parse_args()

land = dict()
hcubes = list()

with open(args.hypercubes) as f:
    for line in f:
        hcubes.append(list())
        hcubes[-1].append(line.split('\t')[1].strip())
        muts = line.split('\t')[0].split(':')
        for mut in muts:
            temp = list()
            for seq in hcubes[-1]:
                temp.append(seq[:int(mut[1:-1])]+mut[0]+seq[int(mut[1:-1])+1:])
            hcubes[-1].extend(temp)

with open(args.landscape) as f:
    for line in f:
        spl = line.split('\t')
        land[spl[0]] = spl[1]

dr = f'data{args.hypercubes[-5]}'
os.mkdir(dr)
os.chdir(dr)

for i in range(len(hcubes)):
    with open(args.hypercubes[-5]+'D_hcube'+str(i), 'w') as f:
        for j in range(len(hcubes[i])):
            f.write(hcubes[i][j]+'\t'+land[hcubes[i][j]])
