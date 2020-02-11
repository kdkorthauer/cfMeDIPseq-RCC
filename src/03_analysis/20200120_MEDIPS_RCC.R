# This script creates MEDIPS objects and carries out various QC and 
# exploratory plotting procedures

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(readxl)
library(dplyr)

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
chr.select <- paste0("chr", c(1:22, "X", "Y", "M"))
#chr.select <- "chr21"

# JAN2020 samples

bam.jan2020 <- list.files(file.path(bamdir, "JAN2020"), "*.sorted.bam", 
	full.names = TRUE)
bam.jan2020 <- bam.jan2020[!grepl(".bai", bam.jan2020)]


fq.jan2020 <- list.files(file.path("/arc/project/st-kdkortha-1", 
	"cfMeDIPseq/data/20200108/128.120.88.251/H202SC19122450/Rawdata"), 
    "*_1.fq.gz", 
    full.names = TRUE,
    recursive = TRUE)

meta2 <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/20200108/20200108_Sample List.xlsx")


excl <- rbind(c("R104_AN",  "<1M reads"),
    c("R41_LD009", "<1M reads"),
    c("R91_233", "<1M reads / Oncocytoma"),
    c("R93_250", "<1M reads / Oncocytoma"),
    c("R54_2321", "Oncocytoma"),
    c("R38_112", "Oncocytoma"),
    c("R92_242", "Oncocytoma"),
    c("R35_190", "Oncocytoma"),
    c("R109_112", "Oncocytoma"),
    c("R113_242", "Oncocytoma"),
    c("R118_233", "Oncocytoma"),
    c("R114_250", "Oncocytoma"),
    c("R12_2350", "missing histology (per Jacob)"),
    c("R67_168", "missing histology (per Jacob)"),
    c("R82_205", "missing histology (per Jacob)"),
    c("R11_2353_L7", "extra lane"),
    c("R12_2350_L7", "missing histology (per Jacob); extra lane"),   
    c("R24_LD021_L7", "extra lane"),  
    c("R42_LD013_L7", "extra lane"),   
    c("R2_CG", "extra coverage (per Sandor/Pier)"),
    c("R6_2366_L7", "extra lane"))
colnames(excl) <- c("string", "reason")

manifest <- data.frame(filename=basename(fq.jan2020), stringsAsFactors=FALSE) %>%
  mutate(lane = substr(filename, nchar(filename)-8, nchar(filename)-8)) %>%
  mutate(ID = gsub("_CKDL.*", "", filename)) %>%
  mutate(ID2 = paste0(manifest$ID, "_L", manifest$lane)) %>%
  mutate(filesize_gb = file.info(fq.jan2020)$size/1e9) %>%
  mutate(barcode = sapply(fq.jan2020, function(x) system(paste0("zcat ", x, " | head -1"), 
  	intern=TRUE))) %>%
  mutate(barcode = substr(barcode, nchar(barcode)-5, nchar(barcode))) %>%
  mutate(inclusion = "include") %>%
  mutate(inclusion = ifelse(substr(filename, 1, 1)!="R", "exclude: Other project (per Jacob)", inclusion)) %>%
  mutate(inclusion = ifelse(grepl("Undetermined", filename), "exclude: Undetermined", inclusion)) %>%
  left_join(data.frame(excl) %>% mutate(ID=string), by = "ID") %>%
  mutate(inclusion = ifelse(!is.na(reason), 
  	paste0("exclude: ", as.character(reason)), inclusion)) %>%
  select(-string, -reason) %>%
  left_join(data.frame(excl) %>% mutate(ID2=string), by = "ID2") %>%
  mutate(inclusion = ifelse(!is.na(reason), 
  	paste0("exclude: ", as.character(reason)), inclusion)) %>%
  select(-string, -reason, -ID2)

write.table(manifest, 
	file=file.path("/arc/project/st-kdkortha-1", 
	"cfMeDIPseq/data/20200108/derived_sample_data.txt"), 
	quote=FALSE,
	row.names = FALSE,
	sep="\t")

# create medip objects

if (!file.exists(file.path(outdir, "medip.jan2020.rds"))){
  medip.jan2020 = lapply(bam.jan2020, function(x) {
	MEDIPS.createSet(file = x,
     BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
     window_size = ws, chr.select = chr.select, paired = TRUE)
  	})

  # exclude certain samples (reasons denoted below)
  samp_name <- sapply(medip.jan2020, function(x) x@sample_name)
  medip.jan2020 <- medip.jan2020[! grepl(paste0(excl[,1], collapse="|"), 
    samp_name)]

  # chop off lane number from sample names 
  for (j in seq_along(medip.jan2020)){
    medip.jan2020[[j]]@sample_name <- gsub("_L..sorted.bam", ".sorted.bam", medip.jan2020[[j]]@sample_name)
  }

  saveRDS(medip.jan2020, file = file.path(outdir, "medip.jan2020.rds"))
}else{
  medip.jan2020 <- readRDS(file.path(outdir, "medip.jan2020.rds"))
}

length(medip.jan2020)

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
