# This script creates MEDIPS objects and carries out various QC and 
# exploratory plotting procedures

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(readxl)

# get windowsize
ws <- Sys.getenv("WINDOWSIZE")

outdir <- paste0("/scratch/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws)

bamdir <- "/arc/project/st-kdkortha-1/cfMeDIPseq/out/sortedbam_dup"
dir.create(outdir, showWarnings = FALSE)
ws <- as.numeric(ws)

bamdir_D <- "/arc/project/st-kdkortha-1/cfMeDIPseq/data/DeCarvalho"

BSgenome="BSgenome.Hsapiens.UCSC.hg19"

# To avoid artefacts caused by PCR over amplification MEDIPS determines a maximal allowed number of stacked reads per genomic position by a poisson distribution of stacked reads genome wide and by a given p-value. The smaller the p-value, the more reads at the same genomic position are potentially allowed. Alternatively, all reads mapping to exactly the same genomic position can be maintained (uniq = 0) or replaced by only one representative (uniq = 1).
uniq=1e-3

# All reads will be extended to a length of 300nt according to the given strand information:

extend=300

# As an alternative to the extend parameter, the shift parameter can be used. Here, the reads are not extended but shifted by the specified number of nucleotides with respect to the given strand infomation. One of the two parameters   or   has to be 0.
shift=0

# canonical chrs
chr.select <- paste0("chr", c(1:22))
#chr.select <- "chr21"

# JAN2020 samples

bam.jan2020 <- list.files(file.path(bamdir, "JAN2020"), "*.sorted.bam", 
	full.names = TRUE)
bam.jan2020 <- bam.jan2020[!grepl(".bai", bam.jan2020)]

#meta2 <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/20200108/20200108_Sample List.xlsx")


# create medip objects

if (!file.exists(file.path(outdir, "medip.jan2020.rds"))){
  medip.jan2020 = lapply(bam.jan2020, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
  	})


  # chop off lane number from sample names 
  for (j in seq_along(medip.jan2020)){
    medip.jan2020[[j]]@sample_name <- gsub("_L..sorted.bam", ".sorted.bam", medip.jan2020[[j]]@sample_name)
  }

  saveRDS(medip.jan2020, file = file.path(outdir, "medip.jan2020.rds"))
}else{
  medip.jan2020 <- readRDS(file.path(outdir, "medip.jan2020.rds"))
}

# exclude certain samples (reasons denoted below)
length_before <- length(medip.jan2020)
samp_name <- sapply(medip.jan2020, function(x) x@sample_name)
medip.jan2020 <- medip.jan2020[! grepl(paste0(c("R104_AN", # <1M reads
    "R41_LD009", # <1M reads
    "R91_233", # <1M reads / Oncocytoma
    "R93_250", # <1M reads / Oncocytoma
    "R54_2321", # Oncocytoma
    "R38_112", # Oncocytoma
    "R92_242", # Oncocytoma
    "R35_190", # Oncocytoma
    "R109_112", # Oncocytoma
    "R113_242", # Oncocytoma
    "R118_233", # Oncocytoma
    "R114_250", # Oncocytoma
    "R12_2350", # missing histology (per Jacob)
    "R67_168", # missing histology (per Jacob)
    "R82_205", # missing histology (per Jacob)
    "R11_2353_L7", # extra lane
    "R12_2350_L7", # extra lane   
    "R24_LD021_L7", # extra lane  
    "R42_LD013_L7", # extra lane   
    "R2_CG", # extra coverage (per Sandor/Pier)
    "R6_2366_L7"), collapse="|"), # extra lane
    samp_name)]

if (length_before > length(medip.jan2020))
  saveRDS(medip.jan2020, file = file.path(outdir, "medip.jan2020.rds"))


# For CpG density dependent normalization of MeDIP-seq data, we need to generate a coupling set. The coupling set must be created based on the same reference genome, the same set of chromosomes, and with the same window size used for the MEDIPS SETs. 
CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip.jan2020[[1]])


# saturation analysis
if (!file.exists(file.path(outdir, "saturation_jan2020.pdf"))){
	pdf(file.path(outdir, "saturation_jan2020.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.jan2020)){
	    sr = MEDIPS.saturation(file = bam.jan2020[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr, main = gsub(".sorted.bam", "", 
		            gsub("/arc/project/st-kdkortha-1/cfMeDIPseq/out/sortedbam_dup/JAN2020/", "", bam.jan2020))[i])
	  }
	dev.off()
}

# coverage
if (!file.exists(file.path(outdir, "CGcoverage_jan2020.pdf"))){
	pdf(file.path(outdir, "CGcoverage_jan2020.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.jan2020)){
	    cr = MEDIPS.seqCoverage(file = bam.jan2020[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("/arc/project/st-kdkortha-1/cfMeDIPseq/out/sortedbam_dup/JAN2020/", "", bam.jan2020))[i])
	  }
	dev.off()
}

# CpG enrichment
#if (!file.exists(file.path(outdir, "CGenrichment_scores_jan2020.pdf"))){
#	er.jan2020 <- vector("list", length(bam.jan2020))
#	for (i in seq_along(er.jan2020)){
#	  er.jan2020[[i]] <- MEDIPS.CpGenrich(file = bam.jan2020[[i]], 
#	  	BSgenome = BSgenome,
#		chr.select = chr.select, extend = extend, shift = shift, 
#		uniq = uniq)
#	}
	
#	pdf(file.path(outdir, "CGenrichment_scores_jan2020.pdf"))
#	hist(sapply(er.jan2020, function(x) x$enrichment.score.relH))
#	hist(sapply(er.jan2020, function(x) x$enrichment.score.GoGe))
#	dev.off()
#}


# CpG calibration plot
if (!file.exists(file.path(outdir, "CGcalibration_jan2020.pdf"))){
	pdf(file.path(outdir, "CGcalibration_jan2020.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.jan2020)){
	  	res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.jan2020[[i]]@sample_name)), 
	          	MSet = medip.jan2020[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.jan2020[[i]]@sample_name)), 
	          	ISet = medip.jan2020[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}
