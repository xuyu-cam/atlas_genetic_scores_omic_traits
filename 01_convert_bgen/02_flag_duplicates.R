library(data.table)

# Fix missing rsIDs and flag duplicate variants for removal.
for (chr_id in 1:22) {
  # Load variant information and add common identifier to help identify duplicates
  # and fix missing rsIDs.
  pvar = fread(sprintf("geno_files/genotype_data/impute_%s_interval.pvar", chr_id))
  pvar[, row := .I]
  pvar[, sorted_alleles := paste(sort(c(REF, ALT)), collapse=":"), by=row]
  pvar[, var_id := paste(`#CHROM`, POS, sorted_alleles, sep=":")]
  pvar[ID == ".", ID := var_id]

  # Load in variant statistics so we can use INFO scores to flag 
  # which of each pair of duplicates to remove
  snpstats = fread(sprintf("/rds/project/jmmh2/rds-jmmh2-projects/polygenic/internal/interval_grs_scan/data/INTERVAL/reference_files/imputed_genotypes/impute_%s_interval.snpstats", chr_id))
  pvar[snpstats, on = .(`#CHROM`=chromosome, POS=position, REF=A_allele, ALT=B_allele), INFO := i.information]
  pvar[snpstats, on = .(`#CHROM`=chromosome, POS=position, REF=B_allele, ALT=A_allele), INFO := i.information]

  # Identify and flag duplicates for removal, keeping the entry with
  # the highest INFO score in each case. There are two cases to deal with:
  #
  # (1) variants that are duplicates by position and alleles
  # (2) variants that are duplicates by rsid (and alleles), but which may
  #     have different positions.
  #
  # The reason we match by alleles as well is that it appears that multi-allelic
  # variants are split into multiple entries even in the BGEN files, so we don't
  # want to incorrectly remove these.

  # First, identify any variant that is duplicate by position, or duplicate by
  # id, then filter the pvar table to all remaining variants ('ok').
  dups_by_pos = pvar[,.N, by=.(var_id)][N > 1]
  dups_by_id = pvar[,.N,by=.(ID, sorted_alleles)][N > 1]
  ok = pvar[!dups_by_pos, on = .(var_id)][!dups_by_id, on=.(ID, sorted_alleles)]

  # Extract the remaining variants, which are all duplicates in either sense
  dups = pvar[!ok, on = .(row)]

  # First, considering all duplicates by position, take the variant with the max
  # INFO score (first one if there are multiple with the same INFO).
  max_info_by_pos = dups[,.SD[which.max(INFO)], by=.(var_id)] 
 
  # Then, from these remaining variants, take the max INFO score by rsID to handle
  # cases where > 1 variant may have the same rsID, but differeing positions.
  max_info_by_id = max_info_by_pos[,.SD[which.max(INFO)], by=.(ID, sorted_alleles)]

  # Add these back to the "ok" table
  ok = rbind(ok, max_info_by_id)

  # Flag in the pvar table the variants to remove - important to preserve order
  # of variants in the output table as the row number corresponds to row in the
  # genotype data.
  pvar[, remove := FALSE]
  pvar[!ok, on = .(row), remove := TRUE]
  
  # Make sure every variant has a unique identifier so we can accurately flag 
  # variants for removal with plink.
  pvar[, ID := paste(ID, row, sep=":")]

  # Flag variants for removal in the identifier
  pvar[(remove), ID := paste(ID, "remove", sep=":")]

  # Overwrite pvar file:
  fwrite(pvar[, .(`#CHROM`, POS, ID, REF, ALT)], sep="\t", quote=FALSE,
         file=sprintf("geno_files/genotype_data/impute_%s_interval.pvar", chr_id))

  # Remove objects and garbage collect before going to next loop
  rm(list=ls())
  gc()
}

