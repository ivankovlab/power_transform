import argparse
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--hypercubes')
parser.add_argument('-l', '--landscape')
parser.add_argument('-i', '--indices')
parser.add_argument('-o', '--output')
args = parser.parse_args()

land = dict()

with open(args.indices) as f:
    lines = f.readlines()
    line = lines[-1]
    line = line.split('[')[2]
    line = line.split(']')[0]
    indices = line.split(', ')

with open(args.landscape) as f:
    for line in f:
        spl = line.split()
        land[spl[0]] = spl[1]

hcubes = list()
genotypes = list()

with open(args.hypercubes) as fin, open(args.output, 'w') as fout:
    for i, line in enumerate(fin):
        if str(i) in indices:
            hcubes.append(list())
            hcubes.append(line.split()[1].strip())
            muts = line.split()[0].split(':')
            for mut in muts:
                temp = list()
                for seq in hcubes[-1]:
                    temp.append(seq[:int(mut[1:-1])]+mut[0]+seq[int(mut[1:-1])+1:])
                hcubes[-1].extend(temp)

    flat_iter = itertools.chain.from_iterable(hcubes)
    genotypes = list(set(list(flat_iter)))

    for genotype in genotypes:
        print(f'{genotype}\t{land[genotype]}\n', file=fout)
