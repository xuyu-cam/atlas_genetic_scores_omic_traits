library(data.table)

args = commandArgs(trailingOnly=TRUE)

dt = fread(sprintf("%s/%s_dosages.txt", args[1], args[2]))
varID = dt[,.(varID=IID)] # should take first IID column encountered
# drop these columns, will occur 1 time per each chromosome file pasted
while ("FID" %in% names(dt)) {
  dt[, c('FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE') := NULL];
}
# add back the varID column as the first column
dt = cbind(varID, dt)
# drop extracted allele from variant name in header
setnames(dt, gsub('_.*', '', names(dt)));
# write out
fwrite(dt, sep='\t', quote=FALSE, compress="gzip", file=sprintf("%s/%s_dosages.txt.gz", args[1], args[2]))
