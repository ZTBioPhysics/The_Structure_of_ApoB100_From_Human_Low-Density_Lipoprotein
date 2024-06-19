This repository contains all the scripts, functions, and data files necessary to reproduce all custom image analysis and crosslink analysis found in the paper

The MatLab script LDL_diamterANDeccentricity.m will perform morphological analysis of cryo-EM 2D class averages of LDL particles

1) download and unzip the class average mrc file
2) download the ReadMRC script from: https://www.mathworks.com/matlabcentral/fileexchange/27021-imagic-mrc-dm-and-star-file-i-o
3) run script

The Python script xlink_analysis_script will reproduce the crosslinking analysis plots and tables

- the module file xlink_analysis_functions is imported and supplies all the functions needed in the script
- the 4 CSV files contain the crosslinking data organized by residue pairs with Ca Distances for both apoB100 models and the average of the two. They also contain the spectral count, sequence distance, and domain associations for each unique crosslinked. 
- the file all_commoni_xlinks_small_and_large.csv contains all crosslinks that were found in common between the two independent datasets and the spectral count column is the average spectral count for that common crosslink between the two datasets
