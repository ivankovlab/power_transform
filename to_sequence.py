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

Sequence can be given either as a string consisting of capital letters only or as a FASTA file.
If a FASTA file is given, only the first sequence from it is considered.
Your FASTA file name must end with either '.fasta' or '.FASTA'.
'''

import argparse
import csv
import sys
import inspect

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename')
parser.add_argument('-s', '--sequence')
args = parser.parse_args()

import re

def parse_fasta(fasta_path):
   """
   Return the first sequence from a FASTA file as a single string.

   Args:
       fasta_path (str): Path to the FASTA file.
       check_invalid_chars (bool): If True, raise error on invalid characters.
       allowed_chars (set): Set of allowed characters (e.g., {'A','T','C','G','N'}).
                             If None and check_invalid_chars is True, default to
                             standard amino acid / nucleotide characters.

   Returns:
       str: The first sequence (concatenated).

   Raises:
       FileNotFoundError: If the file does not exist.
       ValueError: For broken FASTA (no header, empty sequence, malformed, or invalid chars).
   """
   # Common allowed symbols for DNA/RNA/AA (including ambiguity)
   allowed_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*')

   seq_parts = []
   header_found = False
   line_count = 0

   with open(fasta_path, 'r') as f:
       for line in f:
           line = line.strip()
           if not line:
               continue
           line_count += 1

           if line.startswith('>'):
               if header_found:
                   # Second header reached header line
           if header_found:
               # Check for invalid characters 
               invalid = set(line) - allowed_chars
               if invalid:
                   print(
                       f"Invalid character(s) in sequence line {line_count}: {invalid}\n",
                       f"Line: {line}",
                       file=sys.stderr
                   )
                   sys.exit(0)
           seq_parts.append(line)

   # Error conditions
   if not header_found:
       print(f"Broken FASTA: No header line starting with '>' found in '{fasta_path}'", file=sys.stderr)
       sys.exit(0)

   sequence = ''.join(seq_parts)
   if not sequence:
       print(
           f"Broken FASTA: Header found but sequence is empty (no sequence lines) in '{fasta_path}'",
           file=sys.stderr
       )
       sys.exit(0)

   return sequence

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
       if pos - 1 >= len(seq):
           print(inspect.cleandoc(f'''
           ERROR: the sequence is too short.
           The output file was not produced.
           Correct the input sequence and run again.'''),
           file=sys.stderr)
           sys.exit(0)
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

       for row in reader:
           mutations = row[0].strip()
           phenotype = ''.join(row[1:])

           mutant_seq = apply_mutations(wt_sequence, mutations)

           writer.writerow([mutant_seq, phenotype])

if args.sequence.isalpha() and args.sequence.isupper():
   seq = args.sequence
elif args.sequence.endswith('.fasta') or args.sequence.endswith('.FASTA'):
   seq = parse_fasta(args.sequence)
else:
   print(inspect.cleandoc(f'''
   ERROR: this sequence format is not supported.
   The output file was not produced.
   Correct the input sequence and run again.'''),
   file=sys.stderr)
   sys.exit(0)

convert_tsv(args.filename, args.filename.split('.')[0] + '_seq.tsv', seq)
