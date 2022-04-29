# --------------------------------------------------------------------------------------
# Load libraries/dependencies
# --------------------------------------------------------------------------------------
library(data.table)
library(openxlsx)
library(foreach)
library(doMC)

# --------------------------------------------------------------------------------------
# Set global script options
# --------------------------------------------------------------------------------------

out_dir = "geno_files"
ldthinned = "geno_files/genotype_data/ldthinned" # variant set to consider

# Data comes from other projects:
soma_GWAS = "/rds/project/jmmh2/rds-jmmh2-projects/somalogic_proteomics/interval/gwas/BAKEOFF151001/gwas_output/imputed/somalogic/meta"
nmr_GWAS = "/rds/project/jmmh2/rds-jmmh2-projects/nightingale_metabolomics/interval/gwas/nmr/results/HPC_results"
metabo_GWAS = "/rds/project/jmmh2/rds-jmmh2-results/private/metabolomics/metabolon_hd4/interval_gwas/raw_results"
olink_GWAS = "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/jp549/olink-merged-output"
olink_neu_GWAS = "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/interval_gwas_discovery/neu/interval_subset_olink/neuro/full_set/output/formatted_assoc_results"

## ======================================================================================
## First, we want to load the summary statistics for all platforms, and filter to the 
## ld-thinned variant set, and add a basic filter of P < 0.01 - we want to find a 
## reasonable P-value threshold across all platforms but need to load the summary stats
## for all measurements
## ======================================================================================

varset = foreach(chr_id = 1:22, .combine=rbind) %do% {
  fread(sprintf("%s/impute_%s_interval_dedup_unambig_SNPs_maf0.005_ldthin0.8.pvar", ldthinned, chr_id))
}
setnames(varset, c("chr", "pos", "id", "ref", "alt"))

##' # --------------------------------------------------------------------------------------
##' # Load in the NMR GWAS summary stats
##' # ---------------------------------------------------------------------------------------
##' 
##' nmr_files = list.files(path=nmr_GWAS, pattern="*.gz$")
##' nmr_ss = foreach(ff = nmr_files, .combine=rbind) %do% {
##'   # Load from the summary stats variants with P < 0.05 (filtered using awk). 
##'   # The header row is discarded because it is malformed.
##'   ss = fread(cmd = sprintf('zcat %s/%s | tail -n +2 | awk \'BEGIN { OFS="\t" } { if ( $11 < 0.01 ) { print $2,$3,$5,$6,$9,$11 } }\'', nmr_GWAS, ff))
##'   setnames(ss, c("chr", "pos", "effect_allele", "other_allele", "beta", "pval"))
##'   ss = ss[varset[, .(chr, pos)], on = .(chr, pos), nomatch=0]
##'   ss[, phenotype := gsub(".tar.gz", "", ff)]
##'   return(ss)
##' }
##' fwrite(nmr_ss, file=sprintf("%s/nightingale_p_less_0.1.txt", out_dir), sep="\t", quote=F)
##' rm(nmr_ss)
##' gc()
##' 
##' # --------------------------------------------------------------------------------------
##' # Load in the Metabolon HD4 GWAS summary stats
##' # ---------------------------------------------------------------------------------------
##' 
##' metabo_files = list.files(path=metabo_GWAS, pattern="*.gz$")
##' metabo_ss = foreach(ff = metabo_files, .combine=rbind) %do% {
##'   ss = fread(cmd = sprintf('zcat %s/%s | tail -n +2 | awk \'BEGIN { OFS="\t" } { if ( $9 < 0.01 ) { print $1,$2,$3,$7,$9 } }\'', metabo_GWAS, ff))
##'   setnames(ss, c("markername", "effect_allele", "other_allele", "beta", "pval"))
##'   ss[, chr := as.integer(gsub("chr", "", gsub(":.*", "", markername)))]
##'   ss[, pos := as.integer(gsub(":.*", "", gsub("chr.*?:", "", markername)))]
##'   ss = ss[varset[, .(chr, pos)], on = .(chr, pos), nomatch=0]
##'   ss = ss[, .(chr, pos, effect_allele, other_allele, beta, pval)]
##'   ss[, phenotype := gsub("_.*", "", gsub("INTERVAL_", "", ff))]
##'   return(ss)
##' }
##' fwrite(metabo_ss, file=sprintf("%s/metabolon_p_less_0.1.txt", out_dir), sep="\t", quote=F)
##' rm(metabo_ss)
##' gc()

# --------------------------------------------------------------------------------------
# Load in the olink data
# ---------------------------------------------------------------------------------------

olink_files = list.files(path=olink_GWAS, pattern="*.gz$")
olink_ss = foreach(ff = olink_files, .combine=rbind) %do% {
  ss = fread(cmd = sprintf('zcat %s/%s | tail -n +2 | awk \'BEGIN { OFS="\t" } { if ( $22 < 0.01 ) { print $3,$4,$5,$6,$22 } }\'', olink_GWAS, ff))
  setnames(ss, c("chr", "pos", "effect_allele", "other_allele", "pval"))
  ss[, chr := as.integer(chr)]
  ss = ss[varset[, .(chr, pos)], on = .(chr, pos), nomatch=0]
  ss[, phenotype := gsub("_chr_merged.gz", "", gsub("INTERVAL_", "", ff))]
  return(ss)
}

olink_neu_files = list.files(path=olink_neu_GWAS, pattern="*.gz$")
olink_neu_ss = foreach(ff = olink_neu_files, .combine=rbind) %do% {
  ss = fread(cmd = sprintf('zcat %s/%s | tail -n +2 | awk \'BEGIN { OFS="\t" } { if ( $22 < 0.01 ) { print $3,$4,$5,$6,$22 } }\'', olink_neu_GWAS, ff))
  setnames(ss, c("chr", "pos", "effect_allele", "other_allele", "pval"))
  ss[, chr := as.integer(chr)]
  ss = ss[varset[, .(chr, pos)], on = .(chr, pos), nomatch=0]
  ss[, phenotype := paste0("neu_", gsub("_olink.*", "", ff))]
  return(ss)
}

olink_ss = rbind(olink_ss, olink_neu_ss)
fwrite(olink_ss, file=sprintf("%s/olink_p_less_0.1.txt", out_dir), sep="\t", quote=F)
rm(olink_ss, olink_neu_ss)
gc()

# --------------------------------------------------------------------------------------
# Load in the somalogic data
# ---------------------------------------------------------------------------------------

soma_dirs = list.files(path=soma_GWAS)
soma_ss = foreach(dd = soma_dirs, .combine=rbind) %do% {
  soma_files = list.files(path=sprintf("%s/%s", soma_GWAS, dd), pattern="*.gz$")
  foreach(ff = soma_files, .combine=rbind) %do% {
    ss = fread(cmd = sprintf('zcat %s/%s/%s | tail -n +2 | awk \'BEGIN { OFS="\t" } { if ( $8 < -2 ) { print $1,$2,$4,$5,$8 } }\'', soma_GWAS, dd, ff))
    setnames(ss, c("chr", "pos", "effect_allele", "other_allele", "pval"))
    ss[, pval := 10^pval]
    ss = ss[varset[, .(chr, pos)], on = .(chr, pos), nomatch=0]
    ss[, phenotype := dd]
    return(ss)
  }
}

fwrite(soma_ss, file=sprintf("%s/somalogic_p_less_0.1.txt", out_dir), sep="\t", quote=F)



