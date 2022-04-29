# --------------------------------------------------------------------------------------
# Load libraries/dependencies
# --------------------------------------------------------------------------------------
library(data.table)
library(foreach)
library(doMC)
library(AnnotationHub) # Bioconductor package
library(annotables) # remotes::install_github("stephenturner/annotables")

args = commandArgs(trailingOnly=TRUE)
trans_pthresh = as.numeric(args[1])
cis_pthresh = as.numeric(args[2])

if (is.na(trans_pthresh) || trans_pthresh > 0.001 || trans_pthresh < 0) {
  stop("Trans/genome-wide P-value threshold must <= 0.001")
}

if (is.na(cis_pthresh) || cis_pthresh > 1 || cis_pthresh < 0) {
  stop("Cis P-value threshold must be between 0 and 1")
}


ncores = 25

parallelise_fread = function() {
  setDTthreads(ncores)
  registerDoMC(1)
}

parallelise_foreach = function() {
  setDTthreads(1)
  registerDoMC(ncores)
}

# --------------------------------------------------------------------------------------
# Set paths
# --------------------------------------------------------------------------------------

out_dir = "geno_files/ml_inputs"
geno_dir = "geno_files/genotype_data/ldthinned" # variant set to consider
trait_dir = "/rds/project/jmmh2/rds-jmmh2-projects/polygenic/internal/interval_grs_scan/analyses/processed_traits"
trait_dir2 = "/rds/project/jmmh2/rds-jmmh2-projects/polygenic/internal/INTERVAL_gwasqc_technical_covariates_only/qced"

olink_GWAS = "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/jp549/olink-merged-output"
olink_neu_GWAS = "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/interval_gwas_discovery/neu/interval_subset_olink/neuro/full_set/output/formatted_assoc_results"
soma_GWAS = "/rds/project/jmmh2/rds-jmmh2-projects/somalogic_proteomics/interval/gwas/BAKEOFF151001/gwas_output/imputed/somalogic/meta"

processed_GWAS = "geno_files"

tmpdir = sprintf("%s/tmpdir", out_dir)
dir.create(tmpdir, showWarnings=FALSE)

# Define high complexity regions to extend 1MB cis windows when handling protein cis-QTLs
# See flashpca exclusion regions: https://github.com/gabraham/flashpca
# Coordinates are HG19
complex_ld <- data.table(
  region_chr=c(5, 6, 8, 11),
  region_start=c(44000000, 25000000, 8000000, 45000000),
  region_end=c(51500000, 33500000, 12000000, 57000000),
  region_name=c("r1", "MHC", "r3", "r4")
)


## ======================================================================================
## First, we want to load the summary statistics for all platforms, and filter to the
## ld-thinned variant set, and add a basic filter of P < 0.01 - we want to find a
## reasonable P-value threshold across all platforms but need to load the summary stats
## for all measurements
## ======================================================================================

parallelise_fread()
varset = foreach(chr_id = 1:22, .combine=rbind) %do% {
  fread(sprintf("%s/impute_%s_interval_dedup_unambig_SNPs_maf0.005_ldthin0.8.pvar", geno_dir, chr_id))
}
setnames(varset, c("chr", "pos", "id", "ref", "alt"))
setkey(varset, chr, pos)

# --------------------------------------------------------------------------------------
# Load in NMR GWAS results, filter to P < trans_pthresh and output variant effect files
# --------------------------------------------------------------------------------------

if (!file.exists(sprintf("%s/Nightingale_phenotype_info.txt", out_dir))) {
  nmr_SS = fread(sprintf("%s/nightingale_p_less_0.1.txt", processed_GWAS))
  nmr_SS = nmr_SS[pval < trans_pthresh]

  # Make sure we've accurately filtered to varset snps
  nmr_SS = rbind(
    nmr_SS[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=alt, other_allele=ref), nomatch=0],
    nmr_SS[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=ref, other_allele=alt), nomatch=0]
  )
  nmr_SS = nmr_SS[order(pos)][order(chr)][order(phenotype)]

  nmr_SS[varset, on = .(chr, pos), rsid := id] # add rsid
  nmr_SS = nmr_SS[, .SD[which.min(pval)], by=.(chr, pos, phenotype)] # some duplicate results (duplicate SNPs with different INFO). Take best estimate.

  # Load phenotype data and info
  nmr_info = fread(sprintf("%s/nmr_metabolomics/trait_info.tsv", trait_dir))
  nmr_pheno = fread(sprintf("%s/nmr_metabolomics/traits.tsv", trait_dir))

  # Some of the trait names have "_" at the end (measurements with %). Since we use "_" as
  # a file name separator, we'll just replace these 
  nmr_info[, variable := gsub("_", ".pct", variable)]
  nmr_pheno[, variable := gsub("_", ".pct", variable)]
  nmr_SS[, phenotype := gsub("_", ".pct", phenotype)]

  # Filter phenotype data 
  nmr_pheno = nmr_pheno[!is.na(value)]

  # Fix column names:
  setnames(nmr_info, "variable", "PhenotypeCompName")
  setnames(nmr_pheno, "variable", "PhenotypeCompName")
  setnames(nmr_SS, "phenotype", "PhenotypeCompName")

  # Filter info sheet and phenotype data to measurements 
  # with at least 1 variant passing the P-value threshol
  nmr_info = nmr_info[PhenotypeCompName %chin% nmr_SS$PhenotypeCompName]
  nmr_pheno = nmr_pheno[PhenotypeCompName %chin% nmr_SS$PhenotypeCompName]

  # write out variant effects file for each phenotype
  foreach(phen_id = unique(nmr_SS$PhenotypeCompName)) %do% {
    fwrite(nmr_SS[PhenotypeCompName == phen_id, .(rsid, chr, pos, effect_allele, other_allele, effect=beta, pval)],
           sep="\t", quote=FALSE, file=sprintf("%s/%s_variant_effects.txt", out_dir, phen_id))
  }

  # Write out info sheet for NMR variables:
  if (nrow(nmr_info) > 0) {
    fwrite(nmr_info[, .(PhenotypeCompName, Name, Description, Units, Group, Sub.Group)],
           sep="\t", quote=FALSE, file=sprintf("%s/Nightingale_phenotype_info.txt", out_dir))
  }

  # Free objects to free memory
  rm(nmr_info, nmr_SS)
  gc()
} else {
  cat("Nightingale NMR GWAS summary statistics filtered by previous run. Loading and filtering phenotype data...\n")
  nmr_info = fread(sprintf("%s/Nightingale_phenotype_info.txt", out_dir))
  nmr_pheno = fread(sprintf("%s/nmr_metabolomics/traits.tsv", trait_dir))
  nmr_pheno[, variable := gsub("_", ".pct", variable)]
  nmr_pheno = nmr_pheno[!is.na(value)]
  setnames(nmr_pheno, "variable", "PhenotypeCompName")
  nmr_pheno = nmr_pheno[PhenotypeCompName %chin% nmr_info$PhenotypeCompName]
  rm(nmr_info); gc()
} 

# --------------------------------------------------------------------------------------
# Do the same for the metabolon data
# --------------------------------------------------------------------------------------

if (!file.exists(sprintf("%s/Metabolon_phenotype_info.txt", out_dir))) {
  metabo_SS = fread(sprintf("%s/metabolon_p_less_0.1.txt", processed_GWAS))
  metabo_SS = metabo_SS[pval < trans_pthresh]

  # Make sure we've accurately filtered to varset snps
  metabo_SS = rbind(
    metabo_SS[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=alt, other_allele=ref), nomatch=0],
    metabo_SS[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=ref, other_allele=alt), nomatch=0]
  )
  metabo_SS = metabo_SS[order(pos)][order(chr)][order(phenotype)]

  metabo_SS[varset, on = .(chr, pos), rsid := id] # add rsid
  metabo_SS = metabo_SS[, .SD[which.min(pval)], by=.(chr, pos, phenotype)] # some duplicate results (duplicate SNPs with different INFO). Take best estimate.

  # Load phenotype data and info
  metabo_info = fread(sprintf("%s/metabolon_metabolomics/trait_info.tsv", trait_dir))
  metabo_pheno = fread(sprintf("%s/metabolon_metabolomics/traits.tsv", trait_dir))

  # Make it so variable nmae in phenotype data table matches GWAS
  metabo_pheno[, variable := gsub("^m", "M", variable)]

  # Filter phenotype data
  metabo_pheno = metabo_pheno[!is.na(value)]

  # Filter info table to analysed metabolites:
  metabo_info = metabo_info[comp_id %in% metabo_pheno$variable]

  # Fix column names:
  setnames(metabo_info, "comp_id", "PhenotypeCompName")
  setnames(metabo_pheno, "variable", "PhenotypeCompName")
  setnames(metabo_SS, "phenotype", "PhenotypeCompName")

  # Filter info sheet and phenotype data to measurements 
  # with at least 1 variant passing the P-value threshol
  metabo_info = metabo_info[PhenotypeCompName %chin% metabo_SS$PhenotypeCompName]
  metabo_pheno = metabo_pheno[PhenotypeCompName %chin% metabo_SS$PhenotypeCompName]

  # write out variant effects file for each phenotype
  foreach(phen_id = unique(metabo_SS$PhenotypeCompName)) %do% {
    fwrite(metabo_SS[PhenotypeCompName == phen_id, .(rsid, chr, pos, effect_allele, other_allele, effect=beta, pval)],
           sep="\t", quote=FALSE, file=sprintf("%s/%s_variant_effects.txt", out_dir, phen_id))
  }

  # Write out info sheet
  if (nrow(metabo_info) > 0) {
    fwrite(metabo_info[, .(PhenotypeCompName, MASS, RI, BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY)],
           sep="\t", quote=FALSE, file=sprintf("%s/Metabolon_phenotype_info.txt", out_dir))
  }

  # Free objects to free memory
  rm(metabo_info, metabo_SS)
  gc()
} else {
  cat("Metabolon HD4 GWAS summary statistics filtered by previous run. Loading and filtering phenotype data...\n")
  metabo_info = fread(sprintf("%s/Metabolon_phenotype_info.txt", out_dir))
  metabo_pheno = fread(sprintf("%s/metabolon_metabolomics/traits.tsv", trait_dir))
  metabo_pheno[, variable := gsub("^m", "M", variable)]
  metabo_pheno = metabo_pheno[!is.na(value)]
  setnames(metabo_pheno, "variable", "PhenotypeCompName")
  metabo_pheno = metabo_pheno[PhenotypeCompName %chin% metabo_info$PhenotypeCompName]
  rm(metabo_info); gc()
}

# --------------------------------------------------------------------------------------
# Olink data is slightly more complicated: we need to also load the cis-region for each
# protein, and to average proteins measured on multiple platforms
# --------------------------------------------------------------------------------------

if (!file.exists(sprintf("%s/Olink_phenotype_info.txt", out_dir))) {
  olink_SS = fread(sprintf("%s/olink_p_less_0.1.txt", processed_GWAS))

  # Load phenotype data and info
  olink_info = fread(sprintf("%s/olink_proteins/trait_info.tsv", trait_dir))
  olink_pheno = fread(sprintf("%s/olink_proteins/traits.tsv", trait_dir))

  # Make variable names match across data.tables
  olink_SS[, varmatch := tolower(phenotype)]
  olink_SS[, varmatch := gsub("\\.", "", varmatch)]
  olink_SS[, varmatch := gsub("--", "", varmatch)]
  olink_SS[, varmatch := gsub("^inf1", "inf", varmatch)]
  olink_SS[varmatch %like% "inf_dner___q8nft8", varmatch := "inf_dner___q8nft8"] # trailing whitespace
  olink_SS[varmatch == "inf_4ebp1___q13541", varmatch := "inf_ebp1___q13541"]

  # Olink proteins are unique by UniProt identifier, so we will use this as the unique
  # phenotype ID. some fixes are needed.
  olink_info[, PhenotypeCompName := UniProt]
  olink_info[, PhenotypeCompName := gsub(";", ".", PhenotypeCompName)]
  olink_info[PhenotypeCompName == "", PhenotypeCompName := Olink_id]

  # Add unique identifier to olink_SS table so we have full list of genome-wide P < trans_pthresh SNPs for
  # each protein
  olink_SS[olink_info, on = .(varmatch=variable), PhenotypeCompName := i.PhenotypeCompName]

  # Obtain the genomic location of each protein
  loc = olink_info[, .(UniProt = strsplit(UniProt, "\\.")[[1]]), by=.(PhenotypeCompName)]

  system("mkdir -p $HOME/.cache/AnnotationHub", wait=TRUE) # so we dont get interactive prompt below
  ah = AnnotationHub()
  orgdb = query(ah, c("OrgDb", "org.Hs.eg.db"))[[1]]
  txdb <- query(ah, c("TxDB", "TxDb.Hsapiens.UCSC.hg19.knownGene"))[[1]]

  up2gene = select(orgdb, unique(loc$UniProt), c("UNIPROT", "GENENAME", "SYMBOL", "ENTREZID"), "UNIPROT")
  setDT(up2gene)
  up2gene[, ENTREZID := as.integer(ENTREZID)]
  up2gene = up2gene[!is.na(ENTREZID)]
  gene2loc = as.data.table(grch37)
  gene2loc = gene2loc[, .(entrez, chr, start)]
  up2loc = up2gene[gene2loc, on = .(ENTREZID=entrez), nomatch=0]
  up2loc = up2loc[!grepl("_", chr)]
  loc[up2loc, on = .(UniProt=UNIPROT), c("chr", "start") := .(chr, start)]

  # manually annotate a few missing ones by manually looking up UniProt entry
  # and cross referencing with NCBI gene
  loc[UniProt == "Q8WWJ7", c("chr", "start") := .(11, 60739113)]
  loc[UniProt == "Q8NF90", c("chr", "start") := .(4, 81187742)]
  loc[UniProt == "P16284", c("chr", "start") := .(17, 62396775)]
  loc = rbind(loc, data.table(PhenotypeCompName = "OID00195", UniProt=NA, chr=1, start=11917521))

  # Add to olink info table
  olink_info[loc, on = .(PhenotypeCompName), c("chr", "TSS") := .(paste(i.chr, collapse="|"), paste(i.start, collapse="|")), by=.EACHI]

  # Now load summary stats for all SNPs passing P < trans_pthresh genome wide and 
  # cis-snps < cis_pthresh for each protein
  parallelise_foreach()
  foreach(phen_id = unique(olink_SS$PhenotypeCompName), .combine=c) %dopar% {
    if (file.exists(sprintf("%s/%s_variant_effects.txt", out_dir, phen_id))) {
      cat("Olink protein", phen_id, "already filtered at P <", trans_pthresh, "(trans),", cis_pthresh, "(cis) by previous run, skipping.\n")
      return(NULL) 
    }
    gc()
    phen_ss = foreach(var_id = olink_SS[PhenotypeCompName == phen_id, unique(phenotype)], .combine=rbind) %do% {
      # Load full summary stats
      varmatch = olink_SS[PhenotypeCompName == phen_id & phenotype == var_id, unique(varmatch)]
      panel = olink_info[variable == varmatch, panel]
      if (panel == "neu") {
        this_ss = fread(sprintf("%s/%s_olink_neuro_full_set_autosomal_imputed_all_chrs_combined.snptest.out.gz", olink_neu_GWAS, gsub("^neu_", "", var_id)), tmpdir=tmpdir)
      } else {
        this_ss = fread(sprintf("%s/INTERVAL_%s_chr_merged.gz", olink_GWAS, var_id), tmpdir=tmpdir)
      }
      this_ss = this_ss[, .(chr=chromosome, pos=position, effect_allele=alleleB, other_allele=alleleA, effect=frequentist_add_beta_1, pval=frequentist_add_pvalue)] 
      this_ss = this_ss[pval < pmax(trans_pthresh, cis_pthresh)]

      # Filter to varset snps
      if (nrow(this_ss) > 0) {
        this_ss = rbind(
          this_ss[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=alt, other_allele=ref), nomatch=0],
          this_ss[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=ref, other_allele=alt), nomatch=0]
        )
        this_ss = this_ss[!is.na(effect)]
      }
      return(this_ss)
    }

    # R, do you even garbage collect? WTF.
    if(exists("this_ss")) { rm(this_ss) }
    gc()

    # Remove variants that did not pass the p-value threshold for all panels for each protein
    npan = olink_info[PhenotypeCompName == phen_id, length(unique(phenotype))]
    varn = phen_ss[,.N, by=.(chr, pos, effect_allele, other_allele)]
    pass = varn[N == npan]
    pass[,N := NULL]
    phen_ss = phen_ss[pass, on = .(chr, pos, effect_allele, other_allele)]

    # average across measurements
    phen_ss = phen_ss[, .(effect=mean(effect), pval=mean(pval)), by = .(chr, pos, effect_allele, other_allele)]

    # Get SNPs at genome-wide P < trans_pthresh 
    gw = phen_ss[pval < trans_pthresh]

    # Get cis-SNPs with P < cis_pthresh
    chrs = strsplit(olink_info[PhenotypeCompName == phen_id, unique(chr)], "\\|")[[1]]
    starts = strsplit(olink_info[PhenotypeCompName == phen_id, unique(TSS)], "\\|")[[1]]
    cis = foreach(idx = seq_along(chrs), .combine=rbind) %do% {
      window = data.table(chr=as.integer(chrs[idx]), TSS=as.integer(starts[idx]))
      window[, start := pmax(0, TSS - 1e6)]
      window[, end := TSS + 1e6]
      window[complex_ld, on = .(chr=region_chr, start<=region_end, start>=region_start), start := region_start]
      window[complex_ld, on = .(chr=region_chr, end<=region_end, end>=region_start), start := region_start]
      phen_ss[window, on = .(chr, pos>=start, pos<=end), .(chr, pos=x.pos, effect_allele, other_allele, effect, pval)]
    }
    cis = cis[pval < cis_pthresh]

    # remove full phen_ss object to free up memory
    rm(phen_ss)
    gc()

    phen_ss = unique(rbind(gw, cis))

    # If any QTLs, proceed
    if (nrow(phen_ss) > 0) {

      phen_ss = phen_ss[, .SD[which.min(pval)], by=.(chr, pos)] # some duplicate results (duplicate SNPs with different INFO). Take best estimate.
      phen_ss[varset, on = .(chr, pos), rsid := id] # add rsid
      phen_ss = phen_ss[, .(rsid, chr, pos, effect_allele, other_allele, effect, pval)][order(pos)][order(chr)]

      # write out
      fwrite(phen_ss, sep="\t", quote=FALSE, file=sprintf("%s/%s_variant_effects.txt", out_dir, phen_id))

      # Free up memory and garbage collect
      rm(phen_ss)
      gc()

      return(phen_id)
    }
  }
  parallelise_fread()
  gc()

  # Filter phenotype data 
  olink_pheno = olink_pheno[!is.na(value)]

  # Average phenotype data across different platform measures:
  olink_pheno[olink_info, on = .(variable), PhenotypeCompName := PhenotypeCompName]
  olink_pheno = olink_pheno[, .(value=mean(value)), by = .(PhenotypeCompName, IID)]

  # Filter info sheet and phenotype data to measurements 
  # with at least 1 variant passing the P-value threshol
  olink_info = olink_info[PhenotypeCompName %chin% unique(olink_SS$PhenotypeCompName)]
  olink_pheno = olink_pheno[PhenotypeCompName %chin% unique(olink_SS$PhenotypeCompName)]

  # Make sure info table has one entry per phenotype and write out
  if (nrow(olink_info) > 0) {
    olink_info = olink_info[, .(panels = paste(panel, collapse=","), protein = paste(unique(protein), collapse="/")),
                            by = .(PhenotypeCompName, UniProt, chr, TSS)]
    fwrite(olink_info[, .(PhenotypeCompName, UniProt, protein, chr, TSS, panels)],
           sep="\t", quote=FALSE, file=sprintf("%s/Olink_phenotype_info.txt", out_dir))
  }

  # Free objects to free memory
  rm(olink_info, olink_SS)
  gc()
} else {
  cat("Olink protein GWAS summary statistics filtered by previous run. Loading and filtering phenotype data...\n")
  olink_info = fread(sprintf("%s/Olink_phenotype_info.txt", out_dir))
  olink_info_full = fread(sprintf("%s/olink_proteins/trait_info.tsv", trait_dir))
  olink_pheno = fread(sprintf("%s/olink_proteins/traits.tsv", trait_dir))
  olink_info = olink_info[, .(protein=strsplit(protein, "/")[[1]], panel = strsplit(panels, ",")[[1]]), by=.(PhenotypeCompName, UniProt)]
  olink_info[olink_info_full, on = .(UniProt, protein), variable := i.variable]
  olink_pheno = olink_pheno[!is.na(value)]
  olink_pheno[olink_info, on = .(variable), PhenotypeCompName := PhenotypeCompName]
  olink_pheno = olink_pheno[, .(value=mean(value)), by = .(PhenotypeCompName, IID)]
  olink_pheno = olink_pheno[PhenotypeCompName %chin% olink_info$PhenotypeCompName]
  rm(olink_info); gc()
}

# --------------------------------------------------------------------------------------
# Do the same for SomaLogic data. Only this time, we don't have a list of P < 0.01 
# GWAS summary stats to work from.
# --------------------------------------------------------------------------------------

if (!file.exists(sprintf("%s/Somalogic_phenotype_info.txt", out_dir))) {
	# Load phenotype data and info
	soma_info = fread(sprintf("%s/somalogic_proteins/trait_info.tsv", trait_dir))
	soma_pheno = fread(sprintf("%s/somalogic_proteins/traits.tsv", trait_dir))

	# Remove bad aptamers and make unique protien names
	soma_info = soma_info[Type == "Protein"]
	soma_info = soma_info[, .(SeqId, SOMAMER_ID, variable, Target, TargetFullName, UniProt=UniProt.Id.Current.at.Uniprot, Gene=Gene.Name, chr, TSS=start)]

	soma_info[, PhenotypeCompName := gsub("\\..*", "", SOMAMER_ID)]

	# Some are not unique with respect to different Full Target Names, so we need to fix these:
	soma_info[PhenotypeCompName == "VEGFA", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "PLG", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "C3", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "C4A", PhenotypeCompName := Target]
	soma_info[Target == "C5b, 6 Complex", PhenotypeCompName := "C5b"]
	soma_info[Target == "MPIF-1", PhenotypeCompName := "MPIF.1"]
	soma_info[Target == "Ck-b-8-1", PhenotypeCompName := "Ck.b.8.1"]
	soma_info[PhenotypeCompName == "CGA", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "Luteinizing hormone", PhenotypeCompName := "CGA.LHB"]
	soma_info[PhenotypeCompName == "Glycoprotein hormones a-chain", PhenotypeCompName := "CGA"]
	soma_info[Target == "SCGF-beta", PhenotypeCompName := "SCGFb"]
	soma_info[Target == "SCGF-alpha", PhenotypeCompName := "SCGFa"]
	soma_info[Target == "Coagulation Factor Xa", PhenotypeCompName := "F10a"]
	soma_info[Target == "FN1.3", PhenotypeCompName := "FN1.3"]
	soma_info[Target == "Haptoglobin, Mixed Type", PhenotypeCompName := "HPm"]
	soma_info[PhenotypeCompName == "LRP1", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "LYN", PhenotypeCompName := Target]
	soma_info[Target == "PILRA isoform FDF03-M14", PhenotypeCompName := "PILRA.M14"]
	soma_info[Target == "PILRA isoform FDF03-deltaTM", PhenotypeCompName := "PILRA.dTM"]
	soma_info[Target == "Ubiquitin+1", PhenotypeCompName := "RPS27Aplus1"]
	soma_info[Target == "alpha-1-antichymotrypsin complex", PhenotypeCompName := "SERPINA3cmplx"]
	soma_info[Target == "14-3-3 protein beta/alpha", PhenotypeCompName := "14.3.3.pba"]
	soma_info[Target == "14-3-3", PhenotypeCompName := "14.3.3"]
	soma_info[PhenotypeCompName == "C5", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "F2", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "EGFR", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "FN1", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "NRXN1", PhenotypeCompName := Target]
	soma_info[PhenotypeCompName == "ADCYAP1", PhenotypeCompName := gsub("-", ".", Target)]
	soma_info[PhenotypeCompName == "CKB", PhenotypeCompName := gsub("-", ".", Target)]
	soma_info[PhenotypeCompName == "EGFR", PhenotypeCompName := gsub("-", ".", Target)]
	soma_info[PhenotypeCompName == "FGA", PhenotypeCompName := gsub("-", ".", Target)]
	soma_info[PhenotypeCompName == "FGF8", PhenotypeCompName := gsub("-", ".", Target)]
	soma_info[PhenotypeCompName == "PPBP", PhenotypeCompName := gsub("-", ".", Target)]
	soma_info[Target == "CLF-1/CLC Complex", PhenotypeCompName := "CLF1.CLC.complex"]
	soma_info[Target == "CK2-A1:B", PhenotypeCompName := "CK2.A1.B"]
	soma_info[Target == "Coagulation Factor IX", PhenotypeCompName := "CF.IX"]
	soma_info[Target == "Coagulation Factor IXab", PhenotypeCompName := "CF.IXab"]
	soma_info[Target == "GDF-11/8", PhenotypeCompName := "GDF11.8"]
	soma_info[Target == "IgG2, Kappa", PhenotypeCompName := "IgG2"]
	soma_info[Target == "IgG4, Kappa", PhenotypeCompName := "IgG4"]
	soma_info[Target == "N-terminal pro-BNP", PhenotypeCompName := "NPPB.Nt"]
	soma_info[Target == "Activated Protein C", PhenotypeCompName := "PROC.activated"]
	soma_info[Target == "TLR4:MD-2 complex", PhenotypeCompName := "TLR4.MD2.complex"]
	soma_info[Target == "Activin A", PhenotypeCompName := "INHBA.A"]
	soma_info[Target == "Activin AB", PhenotypeCompName := "INHBA.AB"]
	soma_info[Target == "Lymphotoxin a1/b2", PhenotypeCompName := "LTA.A1.B2"]
	soma_info[Target == "Lymphotoxin a2/b1", PhenotypeCompName := "LTA.A2.B1"]
	soma_info[PhenotypeCompName == "POMC", PhenotypeCompName := gsub("-", ".", TargetFullName)]
	soma_info[Target == "SEM6C", Target := "SEMA6C"] 
	soma_info[PhenotypeCompName == "14.3.3", UniProt := "P61981|Q04917"]

  # Actually we want to train models for each SeqId apparently...
  setnames(soma_info, "PhenotypeCompName", "UniqueShortName")
  soma_info[, PhenotypeCompName := paste0("SeqId_", gsub("-", "_", SeqId))]

	# Now load summary stats for all SNPs passing P < trans_pthresh genome wide and 
	# all cis-snps for each protein. Note this code handles averaging across multiple
  # aptamers if we want to return to predicting protein levels.
	parallelise_foreach()
	has_qtls = foreach(phen_id = unique(soma_info$PhenotypeCompName), .combine=c) %dopar% {
    if (file.exists(sprintf("%s/%s_variant_effects.txt", out_dir, phen_id))) {
      cat("Somalogic protein", phen_id, "already filtered at P <", trans_pthresh, "(trans),", cis_pthresh, "(cis) by previous run, skipping.\n")
      return(phen_id) 
    }
		gc()
		# Load all summary stats for all aptamers targetting this protein, filtering to SNPs
		# in the LD-thinned variant set.
		phen_ss = foreach(var_id = soma_info[PhenotypeCompName == phen_id, unique(SOMAMER_ID)], .combine=rbind) %do% {
			this_ss = foreach(chr_id = 1:22, .combine=rbind) %do% {
				chr_ss = fread(sprintf("%s/%s/%s_chrom_%s_meta_1.tbl.gz", soma_GWAS, var_id, var_id, chr_id), tmpdir=tmpdir)
				chr_ss = chr_ss[, .(PhenotypeCompName=phen_id, SOMAMER_ID=var_id, chr=chromosome, pos=position, 
														effect_allele=toupper(Allele1), other_allele=toupper(Allele2), effect=Effect, 
														pval=10^`log(P)`)]
        chr_ss = chr_ss[pval < pmax(trans_pthresh, cis_pthresh)]
        if (nrow(chr_ss) > 0) {
          chr_ss = rbind(
            chr_ss[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=alt, other_allele=ref), nomatch=0],
            chr_ss[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=ref, other_allele=alt), nomatch=0]
          )
        }
				return(chr_ss)
			}
			this_ss[!is.na(effect)]
		}

    # Remove variants that did not pass the p-value threshold for all measurements per protein
    napt = soma_info[PhenotypeCompName == phen_id, length(unique(SOMAMER_ID))]
    varn = phen_ss[,.N, by=.(chr, pos, effect_allele,other_allele)]
    pass = varn[N == napt]
    pass[,N := NULL]
    phen_ss = phen_ss[pass, on = .(chr, pos, effect_allele,other_allele)]

		# average across measurements
		phen_ss = phen_ss[, .(effect=mean(effect), pval=mean(pval)), by = .(chr, pos, effect_allele, other_allele)]

		# R, do you even garbage collect? WTF.
		if(exists("this_ss")) { rm(this_ss) }
		if(exists("chr_ss")) { rm(chr_ss) }
		gc()

		# Identify and extract all SNPs with P < trans_pthresh for all aptamers
		gw_snps = unique(phen_ss[pval < trans_pthresh, .(chr, pos)])
		gw = phen_ss[gw_snps, on = .(chr, pos)]

		# Identify and extract all SNPS in cis with any gene with P < cis_pthresh
		chrs = strsplit(soma_info[PhenotypeCompName == phen_id, unique(chr)], "\\|")[[1]]
		starts = strsplit(soma_info[PhenotypeCompName == phen_id, unique(TSS)], "\\|")[[1]]
		cis = foreach(idx = seq_along(chrs), .combine=rbind) %do% {
			window = data.table(chr=as.integer(chrs[idx]), TSS=as.integer(starts[idx]))
			window[, start := pmax(0, TSS - 1e6)]
			window[, end := TSS + 1e6]
			window[complex_ld, on = .(chr=region_chr, start<=region_end, start>=region_start), start := region_start]
			window[complex_ld, on = .(chr=region_chr, end<=region_end, end>=region_start), start := region_start]
			chr_ss = phen_ss[window, on = .(chr, pos>=start, pos<=end), .(chr, pos=x.pos, effect_allele, other_allele, effect, pval)]
			chr_ss[pval < cis_pthresh]
		}
		rm(phen_ss)
		gc()

		phen_ss = unique(rbind(gw, cis))

		if(nrow(phen_ss) > 0) {
			phen_ss = phen_ss[, .SD[which.min(pval)], by=.(chr, pos)] # some duplicate results (duplicate SNPs with different INFO). Take best estimate.
			phen_ss[varset, on = .(chr, pos), rsid := id] # add rsid
			phen_ss = phen_ss[, .(rsid, chr, pos, effect_allele, other_allele, effect, pval)][order(pos)][order(chr)]

			# write out
			fwrite(phen_ss, sep="\t", quote=FALSE, file=sprintf("%s/%s_variant_effects.txt", out_dir, phen_id))

			# Free up some memory
			rm(phen_ss)
			gc()

			return(phen_id)
		}
	}
	parallelise_fread()

	# Filter phenotype data
	soma_pheno = soma_pheno[variable %chin% unique(soma_info$variable)]
	soma_pheno = soma_pheno[!is.na(value)]

	# Average phenotype data across different platform measures:
	soma_pheno[soma_info, on = .(variable), PhenotypeCompName := PhenotypeCompName]
	soma_pheno = soma_pheno[, .(value=mean(value)), by = .(PhenotypeCompName, IID)]

	# Filter info sheet and phenotype data to measurements 
	# with at least 1 variant passing the P-value threshol
	soma_info = soma_info[PhenotypeCompName %chin% has_qtls]
	soma_pheno = soma_pheno[PhenotypeCompName %chin% has_qtls]

	if (nrow(soma_info) > 0) {
		soma_info = soma_info[, .(SeqId=paste(SeqId, collapse=","), SOMAMER_ID=paste(SOMAMER_ID, collapse=",")),
													 by = .(PhenotypeCompName, UniqueShortName, Target, TargetFullName, UniProt, Gene, chr, TSS)]
		fwrite(soma_info, sep="\t", quote=FALSE, file=sprintf("%s/Somalogic_phenotype_info.txt", out_dir))
	}

	# Adjust phenotype levels for measurement batch:
	if (nrow(soma_pheno) > 0) {
		batch = fread(sprintf("%s/somalogic_proteins/covariates.tsv", trait_dir))
		soma_pheno = soma_pheno[batch, on = .(IID), nomatch=0]
		soma_pheno[, value := lm(value ~ factor(batch))$residuals, by = .(PhenotypeCompName)]
		soma_pheno[, batch := NULL]
	}
} else {
  cat("Somalogic protein GWAS summary statistics filtered by previous run. Loading and filtering phenotype data...\n")
  soma_info = fread(sprintf("%s/Somalogic_phenotype_info.txt", out_dir))
  soma_pheno = fread(sprintf("%s/somalogic_proteins/traits.tsv", trait_dir))
	batch = fread(sprintf("%s/somalogic_proteins/covariates.tsv", trait_dir))
  soma_pheno = soma_pheno[variable %chin% unique(soma_info$variable)]
  soma_pheno = soma_pheno[!is.na(value)]
  soma_pheno = soma_pheno[soma_info[,.(variable, PhneotypeCompName)], on = .(variable), nomatch=0]
  soma_pheno = soma_pheno[, .(value=mean(value)), by = .(PhenotypeCompName, IID)]
  rm(soma_info); gc()
	soma_pheno = soma_pheno[batch, on = .(IID), nomatch=0]
	soma_pheno[, value := lm(value ~ factor(batch))$residuals, by = .(PhenotypeCompName)]
	soma_pheno[, batch := NULL]
}

# --------------------------------------------------------------------------------------
# Combine phenotype data into one table and write out
# --------------------------------------------------------------------------------------

pcs = fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files/annot_INT_50PCs_pcs.txt")

pheno = rbind(
  Nightingale = nmr_pheno,
  Metabolon = metabo_pheno,
  Olink = olink_pheno,
  Somalogic = soma_pheno,
  idcol = "platform")

# Adjust for 10 genotype PCs
pheno = pheno[pcs, on = .(IID=ID), nomatch=0]
pheno[,value := lm(value ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + 
                           PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, 
      by = .(PhenotypeCompName)]

if (nrow(pheno) > 0) {
  fwrite(pheno[, .(PhenotypeCompName, platform, IID, value)],
         sep="\t", quote=FALSE, file=sprintf("%s/phenotypes.txt", out_dir))
}

# Phenotype data adjusted for technical covariates only.
soma_pheno_tech = fread(sprintf("%s/soma4000_gwasQC_adj_technical.txt", trait_dir2))
soma_pheno_tech = melt(soma_pheno_tech, id.vars="aliquot_id", variable.name="SOMAMER_ID")
apt2prot = soma_info[,.(SeqID=strsplit(SeqID, ",")[[1]], SOMAMER_ID=strsplit(SOMAMER_ID, ",")[[1]]), by=PhenotypeCompName]
soma_pheno_tech = soma_pheno_tech[apt2prot, on = .(SOMAMER_ID)]
idmap = fread("/rds/project/jmmh2/rds-jmmh2-projects/polygenic/general/INTERVAL_data/1074/omicsMap.csv")
soma_pheno_tech[idmap, on = .(aliquot_id=soma4000_gwasQC_bl), IID := Affymetrix_gwasQC_bl]
soma_pheno_tech[, platform := "Somalogic"]
soma_pheno_tech = soma_pheno_tech[, .(value=mean(value)), by = .(PhenotypeCompName, platform, IID)]
soma_pheno_tech = soma_pheno_tech[batch, on = .(IID), nomatch=0]
soma_pheno_tech[, value := lm(value ~ factor(batch))$residuals, by = .(PhenotypeCompName)]
soma_pheno_tech = soma_pheno_tech[pcs, on = .(IID=ID), nomatch=0]
soma_pheno_tech[, value := lm(value ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10)$residuals, by = .(PhenotypeCompName)]
soma_pheno_tech = soma_pheno_tech[PhenotypeCompName %chin% pheno[platform == "Somalogic", PhenotypeCompName]]

if (nrow(soma_pheno_tech) > 0) {
  fwrite(soma_pheno_tech[,.(PhenotypeCompName, platform, IID, value)],, sep="\t", quote=FALSE, file=sprintf("%s/phenotypes_no_agesex.txt", out_dir))

  # write out age and sex information
  agesex = fread(sprintf("%s/phenotypes.tsv", trait_dir))
  agesex = agesex[,.(IID, sex=sexPulse, age=agePulse)]
  agesex = agesex[IID %in% unique(soma_pheno_tech$IID)]
  fwrite(agesex, sep="\t", quote=FALSE, file=sprintf("%s/agesex.txt", out_dir))
}


# remove leftover temporary directory
system(sprintf("rm -rf %s", tmpdir), wait=TRUE)

