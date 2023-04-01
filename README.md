[![DOI:10.1038/s41586-023-05844-9](http://img.shields.io/badge/DOI-10.1038/s41586-023-05844-9-B31B1B.svg)](https://www.nature.com/articles/s41586-023-05844-9)

[![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.omicspred.org/)

# An atlas of genetic scores to predict multi-omic traits
This repository houses and documents the codes used to train genetic scores of omic traits using INTERVAL data and internally validate these scores in the study: Xu Y et al. An atlas of genetic scores to predict multi-omic traits. Nature (2023) doi: 10.1038/s41586-023-05844-9 ([https://www.nature.com/articles/s41586-023-05844-9](https://www.nature.com/articles/s41586-023-05844-9)).

All genetic score models trained in the study, and their internal as well as external validation results were all deposited in a cloud service (boxing.com) and are publicly accessible through our online portal (www.omicspred.org).


## The following  softwares and versions were used to perform the analyses:

- Scientific Linux release 7.7 (Nitrogen) (HPC operating system)
- slurm version 19.05.5 (HPC queue manager and job submission system)
- GNU bash version 4.2.46(2) (shell environment used to run bash scripts)
- PLINK v1.90b6.11 64-bit (24 Oct 2019) (www.cog-genomics.org/plink/1.9/)
- PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   (www.cog-genomics.org/plink/2.0/)
- STAR v2.7.3.a (https://github.com/alexdobin/STAR)
- featureCounts v2.0.0 (http://subread.sourceforge.net/)
- QTLtools v1.3.1 (https://qtltools.github.io/qtltools/)

- Python version 3.6.8 with the following Python packages:
   - numpy version 1.19.5 
   - pandas version 1.1.5
   - scikit-learn version 0.21.2 
   - scipy version 1.5.4
   - statsmodels version 0.12.2
   - lifelines version 0.26.0    

- R version 3.6.1 with the following R packages: 
  - cowplot version 1.0.0
  - data.table version 1.13.6
  - dplyr version 1.0.8
  - foreach version 1.5.1
  - ggplot2 version 3.3.5
  - ggpubr version 0.2.5
  - grid version 3.6.1
  - plyr version 1.8.6
  - reshape2 version 1.4.4
  - RcolorBrewer version 1.1-2
  - stringr version 1.4.0
  - tibble version 3.1.0
  - bigsnpr version 1.10.8

## Description of scripts in each sub-folder:

- Genetic score development for multi-omic traits: 
   - **01_convert_bgen**: convert genotype data from bgen to Plink pgen format and remove duplacate variants;
   - **02_ldthin**: remove multi-allelic, ambiguous (A/T, G/G) variants and variants with a MAF < 0.5%, and ld-thin variants with r2=0.8 (i.e. indep-pairwise 1000kb 0.8 in plink2);
   - **03_collate_QTLs**: curate the list of QTLs infomration needed for variant selection from GWAS summary statistics;
   - **04_extract_QTLs**: select the list QTLs with given p-value thresholds, and extract their dosages of the effect alleles as input data of Bayesian Ridge;
   - **05_genetic_score_training**: training genetic score models using Bayesian ridge;
   

- Others: 
   - **06_all_omics_UKB_phecode_assoc_test**: perform PheWAS with the genetic scores of omics traits in UK Biobank
   - **GCTB**: scripts attempted to use SbayesS to estimate heritability of omics traits
   - **LDpred2**: Scripts used to develop genetic scores of omic traits using LDpred2-auto.
