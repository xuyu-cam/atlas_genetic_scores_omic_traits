library(data.table)

chr_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

ref_dir = "geno_files/genotype_data"
out_dir = sprintf("%s/ldthinned", ref_dir)

pvar = fread(sprintf("%s/impute_%s_interval_dedup.pvar", ref_dir, chr_id))

# Remove multi-allelic sites:
multi = pvar[grepl("^rs", ID), .N, by=ID][N > 1]
pvar = pvar[!multi, on = .(ID)]
multi = pvar[!grepl("^rs", ID), .N, by=.(`#CHROM`, POS)][N > 1]
pvar = pvar[!multi, on = .(`#CHROM`, POS)]

# Function for flipping the strand of an allele.
# Uses a series of gsub calls to replace A's with T's,
# G's with C's, and vice-versa. Also works for alleles
# with more than one nucleotide (e.g. indels).
flip_strand <- function(x) {
  # Swap each letter for a dummy, we need this intermediate
  # step so we can distinguish between alleles when swapping.
  # E.g if we did A -> T then T -> A we'd end up with all A's
  # and no T's. instead we do A -> V -> T and T -> X -> A.
  x <- gsub("A", "V", x)
  x <- gsub("T", "X", x)
  x <- gsub("C", "Y", x)
  x <- gsub("G", "Z", x)
  x <- gsub("V", "T", x)
  x <- gsub("X", "A", x)
  x <- gsub("Y", "G", x)
  x <- gsub("Z", "C", x)
  return(x)
}

# Remove strand ambiguous alleles:
pvar = pvar[REF != flip_strand(ALT)]

# Filter to SNPs
pvar = pvar[nchar(REF) == 1 & nchar(ALT) == 1]

# Write out list of variants to extract prior to LD-thinning
fwrite(pvar[,.(ID)], col.names=FALSE, quote=FALSE, file=sprintf("%s/chr%s_keep.txt", out_dir, chr_id))

