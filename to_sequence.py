#!/bin/python3

''' 
This program changes format of the fitness landscapes from mutational view 
to the sequence view.
Mutational view file is as follows (with sequence ABCDE):
  
Genotype	Phenotype
A1F	1.0
C3G	2.0
A1H:D4I	3.0
B2J:D4K:E5L	4.0

Expected sequence view output file is as follows:

Sequence	Phenotype
FBCDE	1.0
ABGDE	2.0
HBCIE	3.0
AJCKL	4.0
'''

import argparse
import csv
import sys
import inspect

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename')
parser.add_argument('-s', '--sequence')
args = parser.parse_args()

def apply_mutations(wt_seq: str, mutations: str) -> str:
    '''
    Apply list of mutations to the reference/wildtype sequence.
    Format of string with mutations is 'A123C:D45E:S125M'
    '''
    if mutations.lower() == 'wt':
        return wt_seq

    seq = list(wt_seq)

    for mut in mutations.split(':'):
        orig_aa = mut[0]
        mut_aa = mut[-1]
        pos = int(mut[1:-1])  # 1-based
        if orig_aa != seq[pos - 1]:
            print(inspect.cleandoc(f'''
                ERROR: amino acid in pos {pos} is {seq[pos - 1]} while \
                    {orig_aa} in mutational list. 
                    The output file was not produced. 
                    Correct the input sequence and run again.'''),
                  file=sys.stderr)
            sys.exit(0)
        seq[pos - 1] = mut_aa

    return ''.join(seq)


def convert_tsv(
    input_tsv: str,
    output_tsv: str,
    wt_sequence: str
):
    
    with open(input_tsv, newline='') as fin, \
         open(output_tsv, 'w', newline='') as fout:

        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        row = next(reader)
        print(' '.join(row), file=fout)
#        print('Sequence', 'Phenotype', sep='\t', newline='')

        for row in reader:
            mutations = row[0].strip()
            phenotype = ''.join(row[1:])

            mutant_seq = apply_mutations(wt_sequence, mutations)

            writer.writerow([mutant_seq, phenotype])


convert_tsv(args.filename, args.filename.split('.')[0] + '_seq.tsv', args.sequence)
