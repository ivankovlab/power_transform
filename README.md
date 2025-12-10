Here are the scripts for application of Power Transform to the GFP random mutagenesis landscape.

The landscapes must be given to the scripts in hash table TSV format.

Format example:

AA  1.0

AC  1.2

CA  1.5

CC  2.4

1) Scripts for Box-Cox method
   

   bc_single.py - the script for application of Box-Cox method to a single combinatorially complete landscape.
   
   Usage: python3 bc_single.py filename.tsv
   

   bc_hcubewise.py - the script for application of Box-Cox method to a random mutagenesis landscape hypercube-wisely.
   
   Usage: python3 bc_hcubewise.py -c hypercubes.txt -l landscape.txt -o output.tsv
   
   hypercubes.txt means the output of application of HypercubeME program to the landscape.
   

   bc_full.py - the script for application of Box-Cox method to a whole simply connected random mutagenesis landscape.
   
   It will not work for GFP random mutagenesis landscape.
   
   Usage: python3 bc_full.py -l landscape.tsv -o output.tsv -m model.tsv
   

2) Scripts for Yeo-Johnston method


   yj_single.py - the script for application of Box-Cox method to a single combinatorially complete landscape.
   
   Usage: python3 yj_single.py filename.tsv
   

   yj_hcubewise.py - the script for application of Box-Cox method to a random mutagenesis landscape hypercube-wisely.
   
   Usage: python3 yj_hcubewise.py -c hypercubes.txt -l landscape.txt -o output.tsv
   
   hypercubes.txt means the output of application of HypercubeME program to the landscape.
   

   yj_full.py - the script for application of Box-Cox method to a whole simply connected random mutagenesis landscape.
   
   Usage: python3 yj_full.py -l landscape.tsv -o output.tsv -m model.tsv
   

3) Scripts for coefficients calculation
   

   coeffs_single.py - the script for calculating the coefficients in a single combinatorially complete landscape.
   
   Usage: python3 coeffs_single.py -l landscape.tsv -o <order> -f output.txt
   
   <order> means an integer number: coefficients of which order need to be calculated.
   

   coeffs_hcubewise.py - the script for calculating the coefficients after running bc_hcubewise.py or yj_hcubewise.py.
   
   Usage: python3 coeffs_hcubewise.py -l landscape.tsv -o <order> -f output.txt
   
   <order> means an integer number: coefficients of which order need to be calculated.
   
   Only landscapes obtained via bc_hcubewise.py of yj_hcubewise.py are supported.
   

   coeffs_full.py - the script for calculating the coefficients in any arbitrary landscape.
   
   Usage: python3 coeffs_full.py -c hypercubes.txt -l landscape.tsv -o <order> -f output.txt
   
   hypercubes.txt means the output of application of HypercubeME program to the landscape.
   
   <order> means an integer number: coefficients of which order need to be calculated.
   

4) Script for quantifying the epistasis - quant.py
   
   Usage: python3 quant.py landscape.tsv
