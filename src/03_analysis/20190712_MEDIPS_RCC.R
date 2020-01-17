# This script creates MEDIPS objects and carries out various QC and 
# exploratory plotting procedures

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(readxl)

# get windowsize
ws <- Sys.getenv("WINDOWSIZE")

outdir <- paste0("/scratch/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws)

bamdir <- "../../out/sortedbam_dup"
dir.create(outdir, showWarnings = FALSE)
ws <- as.numeric(ws)

bamdir_D <- "../../data/DeCarvalho"

BSgenome="BSgenome.Hsapiens.UCSC.hg19"

# To avoid artefacts caused by PCR over amplification MEDIPS determines a maximal allowed number of stacked reads per genomic position by a poisson distribution of stacked reads genome wide and by a given p-value. The smaller the p-value, the more reads at the same genomic position are potentially allowed. Alternatively, all reads mapping to exactly the same genomic position can be maintained (uniq = 0) or replaced by only one representative (uniq = 1).
uniq=1e-3

# All reads will be extended to a length of 300nt according to the given strand information:

extend=300

# As an alternative to the extend parameter, the shift parameter can be used. Here, the reads are not extended but shifted by the specified number of nucleotides with respect to the given strand infomation. One of the two parameters   or   has to be 0.
shift=0

# canonical chrs
chr.select <- paste0("chr", c(1:22, "X", "Y", "M"))
#chr.select <- "chr21"

# bam files
bam.control <- list.files(file.path(bamdir, "CONTROL"), "*.sorted.bam", 
	full.names = TRUE)
bam.control <- bam.control[!grepl(".bai", bam.control)]

bam.rcc <- list.files(file.path(bamdir, "RCC"), "*.sorted.bam", 
	full.names = TRUE)
bam.rcc <- bam.rcc[!grepl(".bai", bam.rcc)]

bam.urineR <- list.files(file.path(bamdir, "URINE_RCC"), "*.sorted.bam", 
	full.names = TRUE)
bam.urineR <- bam.urineR[!grepl(".bai", bam.urineR)]

bam.urineC <- list.files(file.path(bamdir, "URINE_CONTROL"), "*.sorted.bam", 
	full.names = TRUE)
bam.urineC <- bam.urineC[!grepl(".bai", bam.urineC)]


# exclude RCC samples (per Sandor email)
meta.rcc <- read_excel("../../data/RCC/42 RCC fastq files.xlsx")
ids.rcc <- gsub(".sorted.bam", "", 
	            gsub("../../out/sortedbam_dup/RCC/", "", bam.rcc))
bam.rcc.extra <- bam.rcc[!(ids.rcc %in% meta.rcc$`Fastq name`)]
bam.rcc <- bam.rcc[ids.rcc %in% meta.rcc$`Fastq name`]


## Daniel's data
bam.control_D <- list.files(file.path(bamdir_D, "cfMeDIP-Control"), "*.bam", 
	full.names = TRUE)
bam.control_D <- bam.control_D[!grepl(".bai", bam.control_D)]

bam.rcc_D <- list.files(file.path(bamdir_D, "cfMeDIP-RCC"), "*.bam", 
	full.names = TRUE)
bam.rcc_D <- bam.rcc_D[!grepl(".bai", bam.rcc_D)]


## metastatic RCC data
bam.rcc_M <- list.files(file.path(bamdir, "RCCmet"), "*.bam", 
	full.names = TRUE)
bam.rcc_M <- bam.rcc_M[!grepl(".bai", bam.rcc_M)]

## bladder cancer data
bam.blca <- list.files(file.path(bamdir, "BLCA"), "*.bam", 
	full.names = TRUE)
bam.blca <- bam.blca[!grepl(".bai", bam.blca)]

## novartis data (rcc)
bam.nova <- list.files(file.path(bamdir, "Novartis"), "*.bam", 
	full.names = TRUE)
bam.nova <- bam.nova[!grepl(".bai", bam.nova)]

# create medip objects
if (!file.exists(file.path(outdir, "medip.rcc.rds"))){
  medip.rcc = lapply(bam.rcc, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
	  })
  saveRDS(medip.rcc, file = file.path(outdir, "medip.rcc.rds"))
}else{
  medip.rcc <- readRDS(file.path(outdir, "medip.rcc.rds"))
}

if (!file.exists(file.path(outdir, "medip.control.rds"))){
  medip.control = lapply(bam.control, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
  	})
  saveRDS(medip.control, file = file.path(outdir, "medip.control.rds"))
}else{
  medip.control <- readRDS(file.path(outdir, "medip.control.rds"))
}

if (!file.exists(file.path(outdir, "medip.urineR.rds"))){
  medip.urineR = lapply(bam.urineR, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
	})
  saveRDS(medip.urineR, file = file.path(outdir, "medip.urineR.rds"))
}else{
  medip.urineR <- readRDS(file.path(outdir, "medip.urineR.rds"))
}


if (!file.exists(file.path(outdir, "medip.urineC.rds"))){
  medip.urineC = lapply(bam.urineC, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
	})
  saveRDS(medip.urineC, file = file.path(outdir, "medip.urineC.rds"))
}else{
  medip.urineC <- readRDS(file.path(outdir, "medip.urineC.rds"))
}


if (!file.exists(file.path(outdir, "medip.control_D.rds"))){
  medip.control_D = lapply(bam.control_D, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select)
	})
  saveRDS(medip.control_D, file = file.path(outdir, "medip.control_D.rds"))
}else{
  medip.control_D <- readRDS(file.path(outdir, "medip.control_D.rds"))
}


if (!file.exists(file.path(outdir, "medip.rcc_D.rds"))){
  medip.rcc_D = lapply(bam.rcc_D, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select)
	})
  saveRDS(medip.rcc_D, file = file.path(outdir, "medip.rcc_D.rds"))
}else{
  medip.rcc_D <- readRDS(file.path(outdir, "medip.rcc_D.rds"))
}

if (!file.exists(file.path(outdir, "medip.rcc_M.rds"))){
  medip.rcc_M = lapply(bam.rcc_M, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select)
	})
  saveRDS(medip.rcc_M, file = file.path(outdir, "medip.rcc_M.rds"))
}else{
  medip.rcc_M <- readRDS(file.path(outdir, "medip.rcc_M.rds"))
}


if (!file.exists(file.path(outdir, "medip.blca.rds"))){
  medip.blca = lapply(bam.blca, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
  	})
  saveRDS(medip.blca, file = file.path(outdir, "medip.blca.rds"))
}else{
  medip.bcla <- readRDS(file.path(outdir, "medip.blca.rds"))
}



if (!file.exists(file.path(outdir, "medip.nova.rds"))){
  medip.nova = lapply(bam.nova, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
  	})
  saveRDS(medip.nova, file = file.path(outdir, "medip.nova.rds"))
}else{
  medip.nova <- readRDS(file.path(outdir, "medip.nova.rds"))
}



# For CpG density dependent normalization of MeDIP-seq data, we need to generate a coupling set. The coupling set must be created based on the same reference genome, the same set of chromosomes, and with the same window size used for the MEDIPS SETs. 
CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip.rcc[[1]])


# saturation analysis
if (!file.exists(file.path(outdir, "saturation_rcc.pdf"))){
	pdf(file.path(outdir, "saturation_rcc.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.rcc)){
	    sr = MEDIPS.saturation(file = bam.rcc[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr, main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/RCC/", "", bam.rcc))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "saturation_control.pdf"))){
	pdf(file.path(outdir, "saturation_control.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.control)){
	    sr = MEDIPS.saturation(file = bam.control[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr,  main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/CONTROL/", "", bam.control))[i])
	  }
	dev.off()
}


if (!file.exists(file.path(outdir, "saturation_urineC.pdf"))){
	pdf(file.path(outdir, "saturation_urineC.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.urineC)){
	    sr = MEDIPS.saturation(file = bam.urineC[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr,  main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/URINE_CONTROL/", "", bam.urineC))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "saturation_urineR.pdf"))){
	pdf(file.path(outdir, "saturation_urineR.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.urineR)){
	    sr = MEDIPS.saturation(file = bam.urineR[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr,  main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/URINE_RCC/", "", bam.urineR))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "saturation_rccMet.pdf"))){
	pdf(file.path(outdir, "saturation_rccMet.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.rcc_M)){
	    sr = MEDIPS.saturation(file = bam.rcc_M[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr, main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/RCC/", "", bam.rcc_M))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "saturation_blca.pdf"))){
	pdf(file.path(outdir, "saturation_blca.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.blca)){
	    sr = MEDIPS.saturation(file = bam.blca[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr, main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/BLCA/", "", bam.blca))[i])
	  }
	dev.off()
}


if (!file.exists(file.path(outdir, "saturation_nova.pdf"))){
	pdf(file.path(outdir, "saturation_nova.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.nova)){
	    sr = MEDIPS.saturation(file = bam.nova[i], BSgenome = BSgenome,
	     uniq = uniq, extend = extend, shift = shift, window_size = ws,
	     chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	     rank = FALSE)
	    MEDIPS.plotSaturation(sr, main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/Novartis/", "", bam.nova))[i])
	  }
	dev.off()
}


# coverage

if (!file.exists(file.path(outdir, "CGcoverage_rcc.pdf"))){
	pdf(file.path(outdir, "CGcoverage_rcc.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.rcc)){
	    cr = MEDIPS.seqCoverage(file = bam.rcc[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/RCC/", "", bam.rcc))[i])
	  }
	dev.off()
}


if (!file.exists(file.path(outdir, "CGcoverage_control.pdf"))){
	pdf(file.path(outdir, "CGcoverage_control.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.control)){
	    cr = MEDIPS.seqCoverage(file = bam.control[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/CONTROL/", "", bam.control))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcoverage_urineR.pdf"))){
	pdf(file.path(outdir, "CGcoverage_urineR.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.urineR)){
	    cr = MEDIPS.seqCoverage(file = bam.urineR[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/URINE_RCC/", "", bam.urineR))[i])
	  }
	dev.off()
}


if (!file.exists(file.path(outdir, "CGcoverage_urineC.pdf"))){
	pdf(file.path(outdir, "CGcoverage_urineC.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.urineC)){
	    cr = MEDIPS.seqCoverage(file = bam.urineC[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/URINE_CONTROL/", "", bam.urineC))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcoverage_rccMet.pdf"))){
	pdf(file.path(outdir, "CGcoverage_rccMet.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.rcc_M)){
	    cr = MEDIPS.seqCoverage(file = bam.rcc_M[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/RCCmet/", "", bam.rcc_M))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcoverage_blca.pdf"))){
	pdf(file.path(outdir, "CGcoverage_blca.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.blca)){
	    cr = MEDIPS.seqCoverage(file = bam.blca[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/BLCA/", "", bam.blca))[i])
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcoverage_nova.pdf"))){
	pdf(file.path(outdir, "CGcoverage_nova.pdf"), width = 10, height = 10) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.nova)){
	    cr = MEDIPS.seqCoverage(file = bam.nova[[i]], pattern = "CG",     
		  BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
		  shift = shift, uniq = uniq)
	    MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
		  cov.level = c(0,1, 2, 3, 4, 5), main = gsub(".sorted.bam", "", 
		            gsub("../../out/sortedbam_dup/Novartis/", "", bam.nova))[i])
	  }
	dev.off()
}

# CpG enrichment
if (!file.exists(file.path(outdir, "CGenrichment_scores_rcc.pdf"))){
	er.rcc <- vector("list", length(bam.rcc))
	for (i in seq_along(er.rcc)){
	  er.rcc[[i]] <- MEDIPS.CpGenrich(file = bam.rcc[[i]], BSgenome = BSgenome,
		chr.select = chr.select, extend = extend, shift = shift, 
		uniq = uniq)
	}
	
	pdf(file.path(outdir, "CGenrichment_scores_rcc.pdf"))
	hist(sapply(er.rcc, function(x) x$enrichment.score.relH))
	hist(sapply(er.rcc, function(x) x$enrichment.score.GoGe))
	dev.off()
}

if (!file.exists(file.path(outdir, "CGenrichment_scores_control.pdf"))){
	er.control <- vector("list", length(bam.control))
	for (i in seq_along(er.control)){
	  er.control[[i]] <- MEDIPS.CpGenrich(file = bam.control[[i]], BSgenome = BSgenome,
		chr.select = chr.select, extend = extend, shift = shift, 
		uniq = uniq)
	}

	pdf(file.path(outdir, "CGenrichment_scores_control.pdf"))
	hist(sapply(er.control, function(x) x$enrichment.score.relH))
	hist(sapply(er.control, function(x) x$enrichment.score.GoGe))
	dev.off()
}

if (!file.exists(file.path(outdir, "CGenrichment_scores_urineR.pdf"))){
	er.urineR <- vector("list", length(bam.urineR))
	for (i in seq_along(er.urineR)){
	  er.urineR[[i]] <- MEDIPS.CpGenrich(file = bam.urineR[[i]], 
	  	BSgenome = BSgenome,
		chr.select = chr.select, extend = extend, shift = shift, 
		uniq = uniq)
	}

	pdf(file.path(outdir, "CGenrichment_scores_urineR.pdf"))
	hist(sapply(er.urineR, function(x) x$enrichment.score.relH))
	hist(sapply(er.urineR, function(x) x$enrichment.score.GoGe))
	dev.off()
}

if (!file.exists(file.path(outdir, "CGenrichment_scores_urineC.pdf"))){
	er.urineC <- vector("list", length(bam.urineC))
	for (i in seq_along(er.urineC)){
	  er.urineC[[i]] <- MEDIPS.CpGenrich(file = bam.urineC[[i]], 
	  	BSgenome = BSgenome,
		chr.select = chr.select, extend = extend, shift = shift, 
		uniq = uniq)
	}

	pdf(file.path(outdir, "CGenrichment_scores_urineC.pdf"))
	hist(sapply(er.urineC, function(x) x$enrichment.score.relH))
	hist(sapply(er.urineC, function(x) x$enrichment.score.GoGe))
	dev.off()
}

if (!file.exists(file.path(outdir, "CGenrichment_scores_rccMet.pdf"))){
	er.rccM <- vector("list", length(bam.rcc_M))
	for (i in seq_along(er.rccM)){
	  er.rccM[[i]] <- MEDIPS.CpGenrich(file = bam.rcc_M[[i]], 
	  	BSgenome = BSgenome,
		chr.select = chr.select, extend = extend, shift = shift, 
		uniq = uniq)
	}
	
	pdf(file.path(outdir, "CGenrichment_scores_rccMet.pdf"))
	hist(sapply(er.rccM, function(x) x$enrichment.score.relH))
	hist(sapply(er.rccM, function(x) x$enrichment.score.GoGe))
	dev.off()
}

if (!file.exists(file.path(outdir, "CGenrichment_scores_blca.pdf"))){
	er.blca <- vector("list", length(bam.blca))
	for (i in seq_along(er.blca)){
	  er.blca[[i]] <- MEDIPS.CpGenrich(file = bam.blca[[i]], 
	  	BSgenome = BSgenome,
		chr.select = chr.select, extend = extend, shift = shift, 
		uniq = uniq)
	}
	
	pdf(file.path(outdir, "CGenrichment_scores_blca.pdf"))
	hist(sapply(er.blca, function(x) x$enrichment.score.relH))
	hist(sapply(er.blca, function(x) x$enrichment.score.GoGe))
	dev.off()
}

if (!file.exists(file.path(outdir, "CGenrichment_scores_nova.pdf"))){
	er.nova <- vector("list", length(bam.nova))
	for (i in seq_along(er.blca)){
	  er.nova[[i]] <- MEDIPS.CpGenrich(file = bam.nova[[i]], 
	  	BSgenome = BSgenome,
		chr.select = chr.select, extend = extend, shift = shift, 
		uniq = uniq)
	}
	
	pdf(file.path(outdir, "CGenrichment_scores_nova.pdf"))
	hist(sapply(er.nova, function(x) x$enrichment.score.relH))
	hist(sapply(er.nova, function(x) x$enrichment.score.GoGe))
	dev.off()
}


# CpG calibration plot

if (!file.exists(file.path(outdir, "CGcalibration_rcc.pdf"))){
	pdf(file.path(outdir, "CGcalibration_rcc.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.rcc)){
	  	res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.rcc[[i]]@sample_name)), 
	          	MSet = medip.rcc[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.rcc[[i]]@sample_name)), 
	          	ISet = medip.rcc[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcalibration_control.pdf"))){
	pdf(file.path(outdir, "CGcalibration_control.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.control)){
	    res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.control[[i]]@sample_name)),
	          	MSet = medip.control[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.control[[i]]@sample_name)),
	          	ISet = medip.control[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcalibration_urineR.pdf"))){
	pdf(file.path(outdir, "CGcalibration_urineR.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.urineR)){
	    res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.urineR[[i]]@sample_name)), 
	          	MSet = medip.urineR[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.urineR[[i]]@sample_name)),
	          	ISet = medip.urineR[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}


if (!file.exists(file.path(outdir, "CGcalibration_urineC.pdf"))){
	pdf(file.path(outdir, "CGcalibration_urineC.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.urineC)){
	    res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = ppaste0(gsub(".sorted.bam", "", medip.urineC[[i]]@sample_name)),
	          	MSet = medip.urineC[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.urineC[[i]]@sample_name)),
	          	ISet = medip.urineC[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}


if (!file.exists(file.path(outdir, "CGcalibration_rccMet.pdf"))){
	pdf(file.path(outdir, "CGcalibration_rccMet.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.rcc_M)){
	  	res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.rcc_M[[i]]@sample_name)), 
	          	MSet = medip.rcc_M[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.rcc_M[[i]]@sample_name)), 
	          	ISet = medip.rcc_M[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcalibration_blca.pdf"))){
	pdf(file.path(outdir, "CGcalibration_blca.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.blca)){
	  	res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.blca[[i]]@sample_name)), 
	          	MSet = medip.blca[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.bcla[[i]]@sample_name)), 
	          	ISet = medip.blca[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}

if (!file.exists(file.path(outdir, "CGcalibration_nova.pdf"))){
	pdf(file.path(outdir, "CGcalibration_nova.pdf")) 
	  par(mfrow=c(4,4))
	  for (i in seq_along(bam.nova)){
	  	res <- try(MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.nova[[i]]@sample_name)), 
	          	MSet = medip.nova[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE))
	    if(inherits(res, "try-error")){
	      MEDIPS.plotCalibrationPlot(CSet = CS, 
	          	main = paste0(gsub(".sorted.bam", "", medip.nova[[i]]@sample_name)), 
	          	ISet = medip.nova[[i]], 
	          	plot_chr = chr.select, rpkm = TRUE, xrange = TRUE)
	    }
	  }
	dev.off()
}
