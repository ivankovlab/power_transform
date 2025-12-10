Here are the scripts for application of Power Transform to the GFP random mutagenesis landscape.

The landscapes must be given to the scripts in hash table TSV format.
Format example:
AA  1.0
AC  1.2
CA  1.5
CC  2.4

1) Scripts for Box-Cox method

   bc_single.py - the script for application of Box-Cox method to a single combinatorially complete landscape
   Usage: python3 bc_single.py filename.tsv

   bc_hcubewise.py - the script for application of Box-Cox method to a random mutagenesis landscape hypercube-wisely.
   Usage: python3 bc_hcubewise.py -c hypercubes.txt -l landscape.txt -o output.tsv
   hypercubes.txt means the output of application of HypercubeME program to the landscape

   bc_full.py - the script for application of Box-Cox method to a whole simply connected random mutagenesis landscape.
   It will not work for GFP random mutagenesis landscape.
   Usage: python3 bc_full.py -l landscape.tsv -o output.tsv -m model.tsv

2) Scripts for Yeo-Johnston method

   yj_single.py - the script for application of Box-Cox method to a single combinatorially complete landscape
   Usage: python3 yj_single.py filename.tsv

   yj_hcubewise.py - the script for application of Box-Cox method to a random mutagenesis landscape hypercube-wisely.
   Usage: python3 yj_hcubewise.py -c hypercubes.txt -l landscape.txt -o output.tsv
   hypercubes.txt means the output of application of HypercubeME program to the landscape

   yj_full.py - the script for application of Box-Cox method to a whole simply connected random mutagenesis landscape.
   Usage: python3 yj_full.py -l landscape.tsv -o output.tsv -m model.tsv
