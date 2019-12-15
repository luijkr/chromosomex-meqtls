options('stringsAsFactors' = FALSE)
library(MatrixEQTL)
library(BiocParallel)
library(data.table)
library(matrixStats)
library(plyr)

# set variables and paths
gender = 'females'
data_path = '/../..'
out_path = '/../..'

# load methylation and covariate data
load(sprintf('%s/MethylationData/M.%s.Rdata', data_path, tolower(gender)))

# define covariate data
cvrt = t(model.matrix(~ . - ids - run_id - dna_id, data = tmp[, !sapply(tmp, function (x) any(is.na(x)))])[,-1])

# create sliced data
M = SlicedData$new(M)
cvrt = SlicedData$new(cvrt)

# function to load genotype data and run models
run <- function (chr) {
	# load genotypes
	load(sprintf('%s/GenotypeData.%s/chr.%s.Rdata', data_path, gender, chr))

	# split analysis up
	nsplits = round_any(nrow(genotypes)/MaxSNPs, 1, ceiling)
	splits = rep(1:nsplits, length = nrow(genotypes))

	for (i in 1:nsplits) {
		indx = splits == i
		geno = SlicedData$new(genotypes[indx,])

		# define output file
		outfile = sprintf('%s/%s.chr.%s.split.%s.of.%s.txt', out_path, gender, i, nsplits)

		# run models
		empty = Matrix_eQTL_engine(
			snps = geno,
			gene = M,
			cvrt = cvrt,
			output_file_name = outfile,
			pvOutputThreshold = 0.001,
			useModel = modelLINEAR,
			errorCovariance = numeric(),
			verbose = TRUE,
			pvalue.hist = FALSE,
			min.pv.by.genesnp = TRUE,
			noFDRsaveMemory = TRUE
		)

	}
}


chrs = 1:22
MaxSNPs = 45000 # max number of SNPs to test
max.cores = 8 # max cores to use
nworkers = min(max.cores, length(chrs))
BPPARAM = MulticoreParam(workers = nworkers, verbose = TRUE)
empty = bplapply(chrs, FUN=Run, BPPARAM=BPPARAM)
