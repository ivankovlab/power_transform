#!/bin/python3

#this script is used for formatting the fitness landscapes into the hash table format

import argparse
import csv

#Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename')
parser.add_argument('-s', '--sequence')
args = parser.parse_args()

with open(args.filename, newline = '') as fin:
    grand = list()
    st = set()
    reader = csv.reader(fin, delimiter='\t')
    next(reader)
    for row in reader:
        if '*' in row:
            continue
        grand.append(list())
        grand[-1].append(row[1])
        split = row[0].split(':')
        grand[-1].append(list())
        for mut in split:
            grand[-1][-1].append(list())
            grand[-1][-1][-1].append(mut[1:-1])
            grand[-1][-1][-1].append(mut[-1])
            st.add(mut[1:-1])

if '' in st:
    st.remove('')
st_ = set(map(int, st))

seqlst = list(args.sequence)
seqtrunc = list()
for i in range(len(seqlst)):
    if (i+1) in st_:
        seqtrunc.append(seqlst[i])

num = sorted(list(st_))
enum = dict(enumerate(num))
enum_ = dict()
for k, v in enum.items():
    enum_[v] = k

with open(args.filename.split('.')[0] + '_hashtable.tsv', 'w') as fout:
    for i in range(len(grand)):
        seqtrunc_copy = seqtrunc.copy()
        for j in range(len(grand[i][-1])):
            if grand[i][-1][j][1] == 't':
                continue
            seqtrunc_copy[enum_[int(grand[i][-1][j][0])]] = grand[i][-1][j][1]
        fout.write(''.join(seqtrunc_copy) + '\t' + grand[i][0] + '\n')
