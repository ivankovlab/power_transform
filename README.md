Here are the scripts for application of Power Transform to the GFP random mutagenesis landscape.

The approach considered by Sailer and Harms to reduce nonlinearities in fitness landscapes (Sailer, Harms, 2017) uses the Power Transform method with a modification of the Box-Cox power transformation. For each single mutation present in the landscape, the averaged effect across all genetic contexts is computed, and then, based on these effects, an additive phenotype is calculated for each genotype.

Note that the Box-Cox-based Power Transform cannot be applied when either the observed or additive phenotypes take non-positive values. Although measured phenotypes in the GFP landscape are strictly positive, computed additive phenotypes can occasionally become zero or negative during analysis, which makes the Box-Cox transformation inapplicable. To address this limitation, we adopted an alternative formulation of the Power Transform based on the Yeo-Johnson transformation, which can accommodate non-positive values.

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


   yj_single.py - the script for application of Yeo-Johnson method to a single combinatorially complete landscape.
   
   Usage: python3 yj_single.py filename.tsv
   

   yj_hcubewise.py - the script for application of Yeo-Johnson method to a random mutagenesis landscape hypercube-wisely.
   
   Usage: python3 yj_hcubewise.py -c hypercubes.txt -l landscape.txt -o output.tsv
   
   hypercubes.txt means the output of application of HypercubeME program to the landscape.
   

   yj_full.py - the script for application of Yeo-Johnson method to a whole simply connected random mutagenesis landscape.
   
   Usage: python3 yj_full.py -l landscape.tsv -o output.tsv -m model.tsv
   

3) Scripts for coefficients calculation
   

   coeffs_single.py - the script for calculating the coefficients in a single combinatorially complete landscape.
   
   Usage: python3 coeffs_single.py -l landscape.tsv -o order -f output.txt
   
   order means an integer number: coefficients of which order need to be calculated.
   

   coeffs_hcubewise.py - the script for calculating the coefficients after running bc_hcubewise.py or yj_hcubewise.py.
   
   Usage: python3 coeffs_hcubewise.py -l landscape.tsv -o order -f output.txt
   
   order means an integer number: coefficients of which order need to be calculated.
   
   Only landscapes obtained via bc_hcubewise.py of yj_hcubewise.py are supported.
   

   coeffs_full.py - the script for calculating the coefficients in any arbitrary landscape.
   
   Usage: python3 coeffs_full.py -c hypercubes.txt -l landscape.tsv -o order -f output.txt
   
   hypercubes.txt means the output of application of HypercubeME program to the landscape.

   order means an integer number: coefficients of which order need to be calculated.
   

5) Script for quantifying the epistasis - quant.py
   
   Usage: python3 quant.py landscape.tsv
   

6) Script for finding simply connected components in the random mutagenesis landscape - find_max_subset.py

   Usage: python3 find_max_subset.py landscape.tsv
   

7) Script for formatting the fitness landscape file to the hash table format - to_hash.py

   Usage: python3 to_hash.py -f filename.tsv -s sequence

   sequence means the wildtype amino acid or nucleotide sequence


The GFP random mutagenesis landscape (Sarkisyan et al., 2016) was used with all these scripts to obtain all the results published in our paper.
To obtain the GFP random mutagenesis landscape file, you need to download the file named amino_acid_genotypes_to_brightness.tsv from http://dx.doi.org/10.6084/m9.figshare.3102154, and then use the Perl script process_data.pl
Usage: ./process_data.pl amino_acid_genotypes_to_brightness.tsv
Then the file GFP.txt will be created.
