simes = function (q, n) {
	q = sort(q)
	S = n / (1:length(q))
	min(min(q*S), 1)
}

options('stringsAsFactors' = FALSE)
library(data.table)
library(BiocParallel)

# set variables and paths
gender = 'females'
data_path = '/../..'
out_path = '/../..'

# load genotype annotation
load(sprintf('%s/SNPLocations.Rdata', data_path))


# read output files
files = list.files(out_path, pattern = gender, full.names = TRUE)

# read all files
max.cores = 20
nworkers = min(length(files), max.cores)
BPPARAM = MulticoreParam(workers = nworkers, verbose = TRUE)
meQTLs = bplapply(files, FUN=function(f) {
	out = fread(f, header=TRUE, sep='\t')
	names(out) = c('SNP', 'gene', 'beta', 't.stat', 'p.value')
	return(out)
}, BPPARAM=BPPARAM)
meQTLs = do.call(rbind, meQTLs)
meQTLs = merge(meQTLs, SNPLocations[,c('SNPName', 'SNPChr')], by.x='SNP', by.y='SNPName')
meQTLs = split(meQTLs, meQTLs$SNPChr)

# calculate Simes p-value per SNP
n = 10286 # number of CpGs tested
nworkers = 5
BPPARAM = MulticoreParam(workers = nworkers, verbose = TRUE)
tmp = lapply(meQTLs, function(x) subset(x, p.value < .0001))
SNP.PValue = bplapply(tmp, FUN=function(x, n=10286) {
	SNP.PValue = with(x, tapply(p.value, SNP, simes, n=n))
	SNP.PValue = data.frame(
		SNPName = names(SNP.PValue),
		SimesPValue = SNP.PValue
	)
	SNP.PValue = merge(SNP.PValue, SNPLocations[,c('SNPName','SNPChr')], by='SNPName')
	return(SNP.PValue)
}, BPPARAM=BPPARAM)

SNP.PValue = do.call(rbind, SNP.PValue)
SNP.PValue = SNP.PValue[order(SNP.PValue$SimesPValue),]

write.table(SNP.PValue, file=sprintf('%s/%s.SimesPValue.PerSNP.txt', out_path, gender), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
