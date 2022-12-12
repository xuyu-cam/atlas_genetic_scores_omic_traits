library(bigsnpr)
library(data.table)
library(foreach)
library(tictoc)
library(ggplot2)
library(bit64)


setwd("/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/OmicsPred_LDpred_benchmarks")
args = commandArgs(trailingOnly=TRUE)

# the trait name from the gwas
gwas <- args[1]

# Make sure we're running on the appropriate compute node to make use of ramdisks partition
# File back shared memory obejcts can only effectively be used when using the /ramdisks/ partition
# as a temporary working directory, which is only available here on the skylake partitions.
# Trying to use the lustre filesystem means LDpred2 grinds to a halt unless run on a single core
# (and only one instance across the whole cluster) due to the consistency checks made by the lustre
# filesystem on shared memory objects.

if (!(Sys.getenv("SLURM_JOB_PARTITION") %like% "skylake")) {
	stop("Script must be run on compute node on skylake or skylake-himem partitions")
}


## create ramdir
starttime = Sys.time()
ramdir <- "/ramdisks/ldpred2/"
if (dir.exists(ramdir)) system(sprintf("rm -rf %s", ramdir), wait=TRUE)
system(sprintf("mkdir -p %s", ramdir), wait=TRUE)


# Copy genotype data to ramdisks
system(sprintf("cp data/INTERVAL_genotypes/filtered_interval_chr*{rds,bk} %s", ramdir), wait=TRUE)


gwas_file <- sprintf("data/GWAS_SumStats/gwas_filteredCleaned_%s.txt",gwas)
stopifnot(file.exists(gwas_file))

# Setup output directory
outdir <- sprintf("output/ldpred2/train/%s", gwas)
system(sprintf("mkdir -p %s", outdir), wait=TRUE)

# Set up temporary directories - clean up if already exists
tmpdir <- sprintf("tmp/ldpred2/%s", gwas)
if (dir.exists(tmpdir)) system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
system(sprintf("mkdir -p %s", tmpdir), wait=TRUE)

### Do per-SNP QC of GWAS summary statistics
# Load the gwas_ss, match to the genotype data, and obtain per-SNP standard deviations and 
# allele frequencies for downstream SNP QC.
gwas_ss <- fread(gwas_file)
endtime = Sys.time()
data_prep_time = as.numeric(difftime(endtime, starttime, units="secs"))


# This loop needs at least 80GB of memory and takes ~16min to run
tic("#1")
starttime = Sys.time()
gwas_ss <- foreach(this_chr =  1:22, .combine=rbind) %dopar% { 
	cat("\nChromosome:", this_chr, "\n")
	# Attach file-backed genotype data
	geno <- snp_attach(sprintf("%s/filtered_interval_chr%s.rds", ramdir, this_chr))

	## Match summary stats to genotype data. A few notes:
	# - Summary stats for all GWAS have already been filtered to HapMap3 variant
	#   set that intersects with the variants in INTERVAL
	# - Strand orientation of alleles has also already been harmonized to INTERVAL for
	#   all GWAS

	cat("- snp_match \n")
	map <- geno$map[-2]
	names(map) <- c("chr", "rsid", "pos", "a0", "a1")
	
	# snp_match not needed as it's been done before
	matched_snps <- snp_match(as.data.frame(gwas_ss[chr == this_chr]), map, strand_flip=FALSE)	
	setDT(matched_snps)	
	# matched_snps <- gwas_ss[chr == this_chr]

	cat("- allele frequency \n")
	# Obtain allele frequencies in the training data - dosages count 'a0'
	# this is the time consuming step
	matched_snps[, a1freq := Matrix::colSums(geno$genotype[, `_NUM_ID_`], na.rm=TRUE) / (Matrix::colSums(!is.na(geno$genotype[, `_NUM_ID_`]))*2)]
	matched_snps[, test := Matrix::colSums(geno$genotype[, `_NUM_ID_`], na.rm=TRUE) ]

	cat("- sd of allele frequency \n")
	# Compute the standard deviation of the allele frequency
	# See https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html
	matched_snps[, sd_val := sqrt(2 * a1freq * (1 - a1freq))]

	# Return
	return(matched_snps)
}
toc()
endtime = Sys.time()
geno_load_time = as.numeric(difftime(endtime, starttime, units="secs"))
# about 20 mins


# Do Per SNP QC of the summary stats.
# See https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html and
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8016455/

gwas_ss[sd_val > 0, sd_y_est := median(sd_val * beta_se * sqrt(n_eff))]
gwas_ss[sd_val > 0, sd_ss := sd_y_est / (beta_se * sqrt(n_eff))]
gwas_ss[, fail_qc := sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05]

gwas_qc <- gwas_ss[, .(rsid=rsid.ss, chr, pos, effect_allele=a1, other_allele=a0, gwas_beta=beta, gwas_se=beta_se, 
													 n_eff, trainingset_EAF=a1freq, sd_val, sd_y_est, sd_ss, fail_qc)]

fwrite(gwas_qc, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/ldpred2_gwasqc.txt.gz", outdir))


# Diagnostic plot
g <- ggplot(gwas_qc, aes(x=sd_val, y=sd_ss, color=fail_qc)) +
  theme_bigstatsr() +
  geom_point(shape=19, size=0.5, alpha=0.5) +
  geom_abline(intercept=0, slope=1, linetype=2, colour="red") +
  scale_colour_manual(name="SNP failed ldpred2 QC", values=c("TRUE"="purple", "FALSE"="yellow")) +
  xlab("SNP dosage SD in training dataset") + 
  ylab("SNP dosage SD in GWAS") +
  theme(legend.position="bottom")
ggsave(g, width=7.2, height=6, file=sprintf("%s/ldpred2_gwasqc_%s.png", outdir, gwas))


tic("#2") # ~10min (or 40min sometimes...)
starttime = Sys.time()
for (this_chr in 1:22) {
  cat("\nLoading correlation matrix for chromosome ", this_chr, "\n")
  
	# Load correlation matrix precomputed for all candidate variants
	corr_file <- sprintf("data/ldpred2/filtered_interval_ldcorr_chr%s.rds", this_chr)
	stopifnot(file.exists(corr_file))
	corr0 <- readRDS(corr_file)

	# Filter to those passing QC for this GWAS
	cat("- Filter out failed QC \n")
	corr0 <- gwas_ss[chr == this_chr & !(fail_qc), corr0[`_NUM_ID_`, `_NUM_ID_`]]

	# Compute LD score
	cat("- Compute LD score \n")
	gwas_ss[chr == this_chr & !(fail_qc), LDsum := Matrix::colSums(corr0^2)]

	foo <- sum(is.na(gwas_ss[chr == this_chr & !(fail_qc),]$LDsum))
	if (foo > 0) {
		cat("(!) LDsum missing:", foo)
	}

	# Aggregate into a single sparse big matrix
	if (this_chr == 1) { 
    cat("Initialized SFBM\n")
		genocorr <- as_SFBM(corr0, backingfile=sprintf("%s/ldcorr_passqc", tmpdir), compact = TRUE)
	} else {
    cat("Adding matrix to SFBM\n")
		genocorr$add_columns(corr0, nrow(genocorr))
	}
}
toc()
endtime= Sys.time()
corr_load_time = as.numeric(difftime(endtime, starttime, units="secs"))


cat("Moving SFBM backing file to /ramdisks\n")
system(sprintf("cp %s/ldcorr_passqc.sbk %s/", tmpdir, ramdir), wait=TRUE)
system(sprintf("rm %s/ldcorr_passqc.sbk", tmpdir), wait=TRUE)
system(sprintf("ln -s %s/ldcorr_passqc.sbk %s/", ramdir, tmpdir), wait=TRUE)

## Calculate LDSC results (~3min)
tic("LDSC")
starttime = Sys.time()
ldsc <- gwas_ss[!(fail_qc), snp_ldsc(
		ld_score = LDsum, ld_size = .N, 
		chi2 = (beta / beta_se)^2,
		sample_size = n_eff,
		ncores = 1
	)]
toc()
endtime= Sys.time()
ldsc_cal_time = as.numeric(difftime(endtime, starttime, units="secs"))

#save ldsc estimates
write.csv(ldsc,sprintf("%s/ldsc_results.csv", outdir),quote=FALSE)

# Extract estimated heritability
h2_est <- ldsc[["h2"]]

# assgin a small heritability estimate when a negative value is returned in ldsc
if (h2_est<0){
  h2_est=0.001
}


### Run auto model
cat("Running auto model\n")

tic("Auto model")
starttime = Sys.time()
multi_auto <- snp_ldpred2_auto(
	genocorr, gwas_ss[!(fail_qc)], h2_init = h2_est, allow_jump_sign=FALSE,
	vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
	ncores = nb_cores()
)
toc()
endtime= Sys.time()
ldpred_train_time = as.numeric(difftime(endtime, starttime, units="secs"))




# check for "chain" convergence
auto_params <- rbindlist(lapply(multi_auto, function(x) {
	data.table(p_init = x$p_init, h2_init = x$h2_init, p_est = x$p_est, h2_est = x$h2_est)
}))
auto_params[, paramset := .I]

auto_path <- foreach(pIdx = seq_along(multi_auto), .combine=rbind) %do% {
	auto = multi_auto[[pIdx]]
	data.table(paramset = pIdx, path_iter = seq_along(auto$path_p_est), 
						 p_est = auto$path_p_est, h2_est = auto$path_h2_est)
}

g1 <- ggplot(auto_path) + aes(x = path_iter, y=p_est) +
	theme_bigstatsr() + 
	geom_hline(data = auto_params, aes(yintercept=p_est), col="blue") +
	geom_point(shape=19, size=0.5) +
	scale_y_log10(name="p") + xlab("") +
	facet_wrap(~ paramset, ncol=10, labeller = label_both) + 
	theme(strip.background=element_blank(), strip.text=element_text(size=6), 
				axis.text=element_text(size=6), axis.title=element_text(size=10))

g2 <- ggplot(auto_path) + aes(x = path_iter, y=h2_est) +
	theme_bigstatsr() + 
	geom_hline(data = auto_params, aes(yintercept=h2_est), col="blue") +
	geom_point(shape=19, size=0.5) +
	ylab("h2") + xlab("") +
	facet_wrap(~ paramset, ncol=10, labeller = label_both) +
	theme(strip.background=element_blank(), strip.text=element_text(size=6), 
				axis.text=element_text(size=6), axis.title=element_text(size=10))

g <- plot_grid(g1, g2, nrow=2) 
ggsave(g, width=20, height=10, units="in", file=sprintf("%s/ldpred2_auto_chain_convergence.png", outdir))


## select genetic score models to keep
# and use the mean of betas of these selected models as the beta of the final genetic score model with LDpred2-auto
# see https://privefl.github.io/bigsnpr/articles/LDpred2.html
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- (range > (0.95 * quantile(range, 0.95))))
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))


pgs_auto <- foreach(this_chr = 1:22, .combine=`+`) %do% {
  geno <- snp_attach(sprintf("%s/filtered_interval_chr%s.rds", ramdir, this_chr))
  geno <- snp_fastImputeSimple(geno$genotypes, ncores=nb_cores()) # doesn't work with missing genotypes, so need to impute as median
  big_prodVec(
    X = geno, 
    y.col = beta_auto[gwas_ss[!(fail_qc), which(chr == this_chr)]], 
    ind.col = gwas_ss[!(fail_qc) & chr == this_chr, `_NUM_ID_`],
    ncores = nb_cores()
  )
}

# get all sample IDs of the INTERVAL data
geno <- snp_attach(sprintf("%s/filtered_interval_chr%s.rds", ramdir, 1))
sample_IDs = geno$fam$sample.ID

# save calculated genetic scores of the LDpred2-auto model for INTERVAL individuals
auto_pgs_data = data.table(sample_IDs,pgs_auto)
write.table(auto_pgs_data,sprintf("%s/ldpred2_auto_sample_pgs.csv", outdir),quote=FALSE,sep="\t",row.names=FALSE)

# save the genetic score model developed using LDpred2-auto
gwas_ss_sub = gwas_ss[!(fail_qc)]
gwas_ss_sub$auto_beta = beta_auto
auto_pgs_model = gwas_ss_sub[,c('rsid.ss','chr','pos','a1','a0',"auto_beta")]
colnames(auto_pgs_model) =c("rsid",'chr','pos','effect_allele','other_allele','effect')
write_file = sprintf("%s/ldpred2_auto_pgs_model.txt", outdir)
write.table(auto_pgs_model,write_file,quote=FALSE,sep="\t",row.names=FALSE)

# save running times at each stage of the genetic score development with LDpred2-auto
values = c(data_prep_time,geno_load_time,corr_load_time,ldsc_cal_time,ldpred_train_time)
time_name = c('data_prep_time','geno_load_time','corr_load_time','ldsc_cal_time','ldpred_train_time')
df_time <- data.frame(time_name, values)
write_file = sprintf("%s/running_time.txt", outdir)
write.table(df_time,write_file,quote=FALSE,sep="\t",row.names=FALSE)

# save the estimated heritability of all the trained genetic score models with LDpred2-auto above
Kept = keep
Heritability = sapply(multi_auto, get, x="h2_est")
df_heri <- data.frame(Kept, Heritability)
write_file = sprintf("%s/auto_model_heritability.txt", outdir)
write.table(df_heri,write_file,quote=FALSE,sep="\t",row.names=FALSE)



### Clean up

system(sprintf("rm -rf %s", tmpdir), wait=TRUE)
system(sprintf("rm -rf %s", ramdir), wait=TRUE)