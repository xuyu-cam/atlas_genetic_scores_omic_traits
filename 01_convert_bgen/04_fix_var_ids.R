library(data.table)

# Remove extra row identifier tacked onto the variant IDs
for (chr_id in 1:22) {
  pvar = fread(sprintf("geno_files/genotype_data/impute_%s_interval_dedup.pvar", chr_id))
  pvar[, ID := gsub(":[0-9]+?$", "", ID)]
  fwrite(pvar, sep="\t", quote=FALSE, file=sprintf("geno_files/genotype_data/impute_%s_interval_dedup.pvar", chr_id))
}

