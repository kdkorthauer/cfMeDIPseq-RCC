# This script reads in MEDIPS objects (created in 20190712_MEDIPS_RCC.R), 
# identifies DMRs, and uses them to
# investigate training/test splits of the data.

library(MEDIPS)
library(ComplexHeatmap)
library(RColorBrewer)
library(annotatr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Morpho)
library(readxl)
library(ggdendro)
library(matrixStats)
library(stringr)
library(circlize)
library(edgeR)

theme_set(theme_bw())

# source modified MEDIPS functions
source("/arc/project/st-kdkortha-1/cfMeDIPseq/src/03_analysis/20190722_limmamedips.R")

# get windowsize
ws <- as.numeric(Sys.getenv("WINDOWSIZE"))
iter <- as.numeric(Sys.getenv("iter"))
ntop <- as.numeric(Sys.getenv("ntop"))

# dir where binned medips objects of all samples are saved
medipdir <- paste0("/arc/project/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws)
dir.create(medipdir, showWarnings = FALSE)

# outdir to leave-one-out output
outdir <- paste0("/scratch/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws, "/hold_", str_pad(iter, 3, pad = "0"))
dir.create(outdir, showWarnings = FALSE)

outdir_m <- paste0("/scratch/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws, "/hold_m_", str_pad(iter, 3, pad = "0"))
dir.create(outdir_m, showWarnings = FALSE)

bamdir <- "/arc/project/st-kdkortha-1/cfMeDIPseq/out/sortedbam_dup"

BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=1e-3
extend=300
shift=0

ws <- as.numeric(ws)

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

## bladder cancer data
bam.blca <- list.files(file.path(bamdir, "BLCA"), "*.bam", 
  full.names = TRUE)
bam.blca <- bam.blca[!grepl(".bai", bam.blca)]

bam.jan2020 <- list.files(file.path(bamdir, "JAN2020"), "*.sorted.bam", 
  full.names = TRUE)
bam.jan2020 <- bam.jan2020[!grepl(".bai", bam.jan2020)]

# exclude RCC samples (per Sandor email)
meta.rcc <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/RCC/42 RCC fastq files.xlsx")
ids.rcc <- gsub(".sorted.bam", "", 
              gsub("/arc/project/st-kdkortha-1/cfMeDIPseq/out/sortedbam_dup/RCC/", "", bam.rcc))
bam.rcc.extra <- bam.rcc[!(ids.rcc %in% meta.rcc$`Fastq name`)]
bam.rcc <- bam.rcc[ids.rcc %in% meta.rcc$`Fastq name`]


# exclude jan2020 samples failing qc (medip file already filtered)
bam.jan2020 <- bam.jan2020[!grepl(paste0(c("R104_AN", # <1M reads
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
    "R12_2350", # missing histology
    "R67_168", # missing histology
    "R82_205", # missing histology
    "R11_2353_L7", # extra lane
    "R12_2350_L7", # extra lane   
    "R24_LD021_L7", # extra lane  
    "R42_LD013_L7", # extra lane   
    "R6_2366_L7"), collapse="|"), # extra lane
    bam.jan2020)]


## metastatic RCC data
bam.rcc_M <- list.files(file.path(bamdir, "RCCmet"), "*.bam", 
  full.names = TRUE)
bam.rcc_M <- bam.rcc_M[!grepl(".bai", bam.rcc_M)]


set.seed(3874*as.numeric(iter))

# canonical chrs
# exclude X and Y chromosomes 
chr.select <- paste0("chr", c(1:22))
#chr.select <- "chr21"

# read in medip objects (created in 20190712_MEDIPS_RCC.R)
medip.rcc <- readRDS(file.path(medipdir, "medip.rcc.rds"))
medip.control <- readRDS(file.path(medipdir, "medip.control.rds"))
medip.urineR <- readRDS(file.path(medipdir, "medip.urineR.rds"))
medip.urineC <- readRDS(file.path(medipdir, "medip.urineC.rds"))
medip.jan2020 <- readRDS(file.path(medipdir, "medip.jan2020.rds")) # already filtered
medip.rcc_M <- readRDS(file.path(medipdir, "medip.rcc_M.rds")) # metastatic

# remove control plasma samples with failed qc (S040) 
medip.control <- medip.control[!sapply(medip.control, function(x) x@sample_name) %in% c("S040.sorted.bam")]
bam.control <- bam.control[!grepl("S040.sorted.bam", bam.control)]

# remove RCC plasma samples with failed qc (S011) and previous
# cancer (S007)
medip.rcc <- medip.rcc[!sapply(medip.rcc, function(x) x@sample_name) %in% c("S007.sorted.bam", "S011.sorted.bam", "S020.sorted.bam")]
bam.rcc <- bam.rcc[!grepl("S007.sorted.bam|S011.sorted.bam|S020.sorted.bam", bam.rcc)]

# remove two chromophobe RCC samples
medip.rcc <- medip.rcc[!sapply(medip.rcc, function(x) x@sample_name) %in% c("S015.sorted.bam", "S050.sorted.bam")]
bam.rcc <- bam.rcc[!grepl("S015.sorted.bam|S050.sorted.bam", bam.rcc)]

# remove samples with too few reads 
# (Urine control S56 too few reads)
# S4, S6, S37 failed qpcr efficiency assay
medip.urineC <- medip.urineC[!sapply(medip.urineC, function(x) x@sample_name) %in% c("S56.sorted.bam", "S4.sorted.bam", "S6.sorted.bam", "S37.sorted.bam")]
bam.urineC <- bam.urineC[!grepl("S56.sorted.bam|S4.sorted.bam|S6.sorted.bam|S37.sorted.bam", bam.urineC)]


# (Urine RCC S11, S12, S30, S33, S34, S48, S59, S60 too few reads) 
# (S62, S64 not saturated)
# (S42 S14 <25% mapping rate & < 40M reads)
# (S25 old urine)
medip.urineR <- medip.urineR[!sapply(medip.urineR, function(x) x@sample_name) %in%
    c("S11.sorted.bam", "S12.sorted.bam", "S30.sorted.bam", "S33.sorted.bam", 
	  "S34.sorted.bam", "S48.sorted.bam", "S59.sorted.bam", "S60.sorted.bam",
	  "S62.sorted.bam", "S64.sorted.bam", "S42.sorted.bam", "S25.sorted.bam",
    "S14.sorted.bam")]
bam.urineR <- bam.urineR[!grepl("S11.sorted.bam|S12.sorted.bam|S30.sorted.bam|S33.sorted.bam|S34.sorted.bam|S48.sorted.bam|S59.sorted.bam|S60.sorted.bam|S62.sorted.bam|S64.sorted.bam|S42.sorted.bam|S25.sorted.bam|S14.sorted.bam", bam.urineR)]

meta <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/RCC/Keegan - RCC plasma.xlsx", 
  sheet = 1)
meta2 <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/20200108/20200108_Sample List.xlsx")
CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip.jan2020[[1]])

m2 <- data.frame(ID=gsub(".sorted.bam", "", 
  sapply(medip.jan2020, function(x) x@sample_name))) %>%
  left_join(meta2, by = "ID")

# master metadata
master <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/20200108/20.01.31 - Final Sample List for NM Revisions.xlsx", 
  sheet = 3)
master <- master %>%
  na_if("N/A") %>%
  mutate(Histology = tolower(Histology)) %>%
  mutate(`Sample number` = ifelse(Batch == "Met", tolower(ID), ID)) %>%
  mutate(Histology = ifelse(Histology %in% c("collecting duct", "chrcc", "xptranslocation"), 
    "other", Histology))

# summary plots - make for all sets
# for comparison 
if (iter == 1){

# how much coverage in windows with no CpGs?
noCpGcov <- function(mdobjlist, CS, type){
  pctNoCpG <- sapply(mdobjlist, function(x){
  	  sum(x@genome_count[CS@genome_CF == 0], na.rm=TRUE) / 
    sum(x@genome_count)
  })
  pctNonzeroCov <- sapply(mdobjlist, function(x){
  	  sum(x@genome_count > 0) / length(x@genome_count)
  })

  data.frame(pctNoCpG = pctNoCpG,
  	pctNonzeroCov = pctNonzeroCov,
  	type = type)
}

df <- rbind(noCpGcov(medip.jan2020[m2$Source=="Plasma" & m2$Status == "RCC"],
                      CS=CS, "rcc_plasma"),
  noCpGcov(medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"], 
                      CS=CS, "ctrl_plasma"),
  noCpGcov(medip.jan2020[m2$Source=="Urine" & m2$Status == "RCC"],
                      CS=CS, "rcc_urine"),
  noCpGcov(medip.jan2020[m2$Source=="Urine" & m2$Status == "Control"],
                      CS=CS, "control_urine"))
ggplot(df, aes(x = pctNoCpG*100)) +
  geom_histogram() +
  facet_wrap(~type) +
  xlab("Percent of total reads mapping to bins with 0 CpGs")
ggsave(file.path(outdir, "../PctReadsinBinswithZeroCpGs_jan2020.pdf"))

ggplot(df, aes(x = pctNonzeroCov*100)) +
  geom_histogram() + 
  xlab("Percent of bins with nonzero coverage")+
  facet_wrap(~type) 
ggsave(file.path(outdir, "../PctBinswithNonzeroCov_jan2020.pdf"))


# pca plots
depths <- function(mdobjlist, CS, type){
  depth <- data.frame(sapply(mdobjlist, function(x){
  	  x@genome_count[CS@genome_CF > 0]
  }))
  colnames(depth) <- sapply(mdobjlist, function(x){ 
  	gsub(".bam|.sorted.bam", "", x@sample_name)})
  colnames(depth) <- gsub("_", "", colnames(depth))
  colnames(depth) <- paste0(type, "_", colnames(depth))
  return(depth)
}

if (!file.exists(file.path(outdir, "../PCA_2_3_combined_urine.pdf"))){
  
  # JAN2020 plasma

  df <- cbind(depths(medip.jan2020[m2$Source=="Plasma" & m2$Status == "RCC"],
                      CS=CS, "rcc_plasma"),
  depths(medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"], 
                      CS=CS, "control_plasma"))

  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  sf <- DESeq2::estimateSizeFactorsForMatrix(df)
  df <- sweep(df, MARGIN=2, sf, `/`)

  grp <- gsub("_R.*", "", colnames(df))
  ids <- gsub("rcc_urine_|control_urine_|rcc_plasma_|control_plasma_", "", colnames(df))

  pcs <- Morpho::prcompfast(t(log(df+0.01)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           id = ids) 
  
  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = type), size = 2)  
  ggsave(file.path(outdir, "../PCA_1_2_jan2020_plasma.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = type)) +
    geom_point(size = 2)
  ggsave(file.path(outdir, "../PCA_2_3_jan2020_plasma.pdf"), width = 4.5, height = 3.5)
  
  
  # JAN2020 urine

  df <- cbind(depths(medip.jan2020[m2$Source=="Urine" & m2$Status == "RCC"],
                      CS=CS, "rcc_urine"),
  depths(medip.jan2020[m2$Source=="Urine" & m2$Status == "Control"],
                      CS=CS, "control_urine"))

  # remove R99 (>16% outlier for pct reads in bins with no CpGs)
  df <- df[,-which(grepl("R99", colnames(df)))]  

  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  sf <- DESeq2::estimateSizeFactorsForMatrix(df)
  df <- sweep(df, MARGIN=2, sf, `/`)

  grp <- gsub("_R.*", "", colnames(df))
  ids <- gsub("rcc_urine_|control_urine_|rcc_plasma_|control_plasma_", "", colnames(df))

  pcs <- Morpho::prcompfast(t(log(df+0.01)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           id = ids) 
  
  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = type), size = 2)  
  ggsave(file.path(outdir, "../PCA_1_2_jan2020_urine.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = type)) +
    geom_point(size = 2)
  ggsave(file.path(outdir, "../PCA_2_3_jan2020_urine.pdf"), width = 4.5, height = 3.5)
  



  # JAN2020 Plus Original cohort (plasma)

   df <- cbind(depths(medip.jan2020[m2$Source=="Plasma" & m2$Status == "RCC"],
                      CS=CS, "rccNew"),
    depths(medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"], 
                      CS=CS, "controlNew"),
    depths(medip.rcc, CS, "rcc"),
    depths(medip.control, CS, "control"))
  
  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  sf <- DESeq2::estimateSizeFactorsForMatrix(df)
  df <- sweep(df, MARGIN=2, sf, `/`)

  grp <- gsub("_R.*|_S.*", "", colnames(df))
  ids <- gsub("rcc_urine_|control_urine_|rcc_plasma_|control_plasma_", "", colnames(df))

  pcs <- Morpho::prcompfast(t(log(df+0.01)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           id = ids) 
  
  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = type), size = 2)  
  ggsave(file.path(outdir, "../PCA_1_2_combined_plasma.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = type)) +
    geom_point(size = 2)
  ggsave(file.path(outdir, "../PCA_2_3_combined_plasma.pdf"), width = 4.5, height = 3.5)

  # combined urine
  df_urine <- cbind(depths(medip.jan2020[m2$Source=="Urine" & m2$Status == "RCC"],
                      CS=CS, "urineR_new"),
  depths(medip.jan2020[m2$Source=="Urine" & m2$Status == "Control"],
                      CS=CS, "urineC_new"),
    depths(medip.urineR, CS, "urineR"),
    depths(medip.urineC, CS, "urineC"))

  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  sf <- DESeq2::estimateSizeFactorsForMatrix(df)
  df <- sweep(df, MARGIN=2, sf, `/`)

  grp <- gsub("_R.*|_S.*", "", colnames(df))
  ids <- gsub("rcc_urine_|control_urine_|rcc_plasma_|control_plasma_", "", colnames(df))

  pcs <- Morpho::prcompfast(t(log(df+0.01)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           id = ids) 
  
  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = type), size = 2)  
  ggsave(file.path(outdir, "../PCA_1_2_combined_urine.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = type)) +
    geom_point(size = 2)
  ggsave(file.path(outdir, "../PCA_2_3_combined_urine.pdf"), width = 4.5, height = 3.5)

}
}


# obj1 = first medip object
# obj2 = second medip object
# lab1 = label of first medip obj
# lab2 = label of second medip obj
# validate1 = first validation medip obj (optional)
# validate2 = second validation medip obj (optional)
# sig.level = qvalue threshold for defining DMRs (optional, not used when top is specified)
# mdip.opt = logical whether to adjust for gc content bias in counts
# chrs = vector of chromosome names to include
# add_samp = sample to add to heatmap plotting
# lab_add = label for additional sample specified by add_samp
# training = proportion of obj1 and obj2 to use as training set (the rest will be used as test)
# regulatory = logical whether to restrict to regulatory regions
# top = threshold for number of DMRs (take top `top` DMRs ranked by qval) - will be appended to results files (heatmaps, results tables)
# colnames = logical whether to include sample names in heatmaps
# heatmap_name_addendum = character string to add to heatmap name if generating a nonstandard one (e.g. if `add_samp` is not null)
# merge = logical whether to allow for building larger regions by grouping together neighboring top regions
# out.dir = directory location to save results
# validate_lab = label for files saved in validation (heatmap, results table)
# saveprobs = logical whether to save sample-level probabilities in a table
# holdout = integer for which element (sample) to hold out; names() slot must contain either 'obj1' or 'obj2' designating which group the heldout sample is to be chosen.

# differential coverage
compute.diff <- function(obj1 = NULL, obj2 = NULL,
	                     lab1 = NULL, lab2 = NULL,
                       validate1 = NULL, validate2 = NULL,
	                     sig.level = NULL,
	                     mdip.opt = FALSE,
	                     chrs = chr.select,
	                     add_samp = NULL,
	                     lab_add = NULL,
	                     training = 1, 
	                     regulatory = FALSE,
	                     top = 300,
	                     colnames = TRUE,
	                     heatmap_name_addendum = NULL,
                       merge = FALSE,
                       out.dir = outdir, 
                       validate_lab = NULL,
                       saveprobs=TRUE,
                       holdout=NULL){

    set.seed(20190814/as.numeric(iter))

    if(grepl("urine", lab1) || grepl("urine", lab2)){
    	mdip.opt <- FALSE
    	if(mdip.opt)
    		message("Forcing mdip.opt to FALSE for urine samps ",
    			"since some are not well-behaved and throw an error ",
    			"when trying to fit the CG density curve")
    }

    if(!is.null(validate1) | !is.null(validate2)){
      heatmap.file.test <- file.path(out.dir, 
                            paste0(lab1, ".", lab2, ".diff.validate_",
                                   validate_lab, ".heatmap.top", top,".pdf"))
      if(training < 1)
        stop("Don't specify training percentage less than ",
                "1 if carrying out a validation")
    }

    if (!is.null(sig.level)){
    	message("Using ", sig.level, " FDR cutoff instead of taking top ", 
    		top, " regions.")
	    }else{
	    	message("Using top ", top , " regions instead of ",
	    		"FDR cutoff")
    }
    diff.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".diff.rds"))
    heatmap.file <- file.path(out.dir,paste0(lab1, ".", lab2, ".diff.heatmap.top", top,".pdf"))
    if (!is.null(heatmap_name_addendum))
      heatmap.file.add <- file.path(out.dir,paste0(lab1, ".", lab2, ".diff.heatmap.",
      	heatmap_name_addendum, ".pdf"))
    if (regulatory){
    	heatmap.file <- file.path(out.dir,
    		paste0(lab1, ".", lab2, ".diff.regulatory.heatmap.pdf"))
    	
    }
    n1 <- min(max(3, round(training*length(obj1))), length(obj1))
    n2 <- min(max(3, round(training*length(obj2))), length(obj2))
    ix1 <- 1:n1
    ix2 <- 1:n2

    prob.table.file <- file.path(out.dir,
      paste0("sampleprob_table_", lab1, "_", lab2, "_top", top, ".txt"))
    
    if (training < 1){
      diff.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".diff.training_",
      	                                    training, ".rds"))
      heatmap.file <- file.path(out.dir, 
      	                        paste0(lab1, ".", lab2, ".diff.training_",
      	                                    training, ".heatmap.top", top,".pdf"))
      heatmap.file.test <- file.path(out.dir, 
      	                        paste0(lab1, ".", lab2, ".diff.test_",
      	                                    1-training, ".heatmap.top", top,".pdf"))

      if(regulatory){
      	heatmap.file <- file.path(out.dir, 
      	                        paste0(lab1, ".", lab2, ".diff.training_",
      	                                    training, ".heatmap.regulatory.pdf"))
      	heatmap.file.test <- file.path(out.dir, 
      	                        paste0(lab1, ".", lab2, ".diff.test_",
      	                                    1-training, ".heatmap.regulatory.pdf"))
      }
      n1 <- min(n1, length(obj1)-1)
      n2 <- min(n2, length(obj2)-1)
      ix1 <- sample(1:length(obj1), n1)
      ix2 <- sample(1:length(obj2), n2)
    }else if (!is.null(holdout)){
      heatmap.file.test <- file.path(out.dir, 
                                paste0(lab1, ".", lab2, ".diff.test",
                                       ".heatmap.top", top,".pdf"))

      training <- 0
      if (names(holdout) == "obj1"){
        ix1 <- ix1[-as.numeric(holdout)]
        n1 <- n1-1
      }else if (names(holdout) == "obj2"){
        ix2 <- ix2[-as.numeric(holdout)]
        n2 <- n2-1
      }else{
        stop("names() of holdout must be either obj1 or obj2")
      }
    }

  if (merge)
        heatmap.file <- gsub(".pdf", ".merge.pdf", heatmap.file)
  if (merge & training < 1)
        heatmap.file.test <- gsub(".pdf", ".merge.pdf", heatmap.file.test)
 
	if (!file.exists(diff.file)){
    diff = MEDIPS.meth(MSet1 = obj1[ix1], MSet2 = obj2[ix2],
	          CSet = CS, diff.method = "limma", chr = chrs,
		        p.adj = "BH", MeDIP = mdip.opt, minRowSum = 0.2*(n1+n2))
	  saveRDS(diff, file = diff.file)
	}else{
	  diff <- readRDS(file=diff.file)
	}
  
  if (ncol(diff) != 3+2*(n1+n2)+11)
    message("WARNING; diff has ", ncol(diff), " columns. ",
      "Expecting ", 3+2*(n1+n2)+11, ".")

	colnames(diff)[3] <- "end"

	# restrict to regulatory regions
	if (regulatory){
      # get cpg and enhancer annots
      annots <- c("hg19_cpgs", "hg19_enhancers_fantom")
      annotations = build_annotations(genome = 'hg19', annotations = annots)
      annotations <- annotations[annotations$type != "hg19_cpg_inter",]
      ol <- findOverlaps(makeGRangesFromDataFrame(diff), annotations)
      diff <- diff[unique(ol@from),]
      diff$limma.adj.p.value <- p.adjust(diff$limma.p.value)
  }

	# make heatmap
  message("extracting top DMRs")
  message("there are ", 
    sum(diff$limma.adj.p.value < 0.05 & !is.na(diff$P.Value)),
    " dmrs at FDR 0.05.")
  message("there are ", 
    sum(diff$limma.adj.p.value < 0.01 & !is.na(diff$P.Value)),
    " dmrs at FDR 0.01.")

	if (!is.null(sig.level)){
	 	which.sig <- which(diff$limma.adj.p.value < sig.level & 
		               abs(diff$logFC) > 2 &
		               !is.na(diff$P.Value) )
	}else{# top ntop - half up and half down
    which.up <- which(diff$logFC > 0)
    which.down <- which(diff$logFC < 0)

		which.sig.up <- which(rank(diff$P.Value[which.up], 
      ties.method = "random") <= as.numeric(top)/2)

    which.sig.down <- which(rank(diff$P.Value[which.down], 
      ties.method = "random") <= as.numeric(top)/2)

    which.sig <- c(which.up[which.sig.up], which.down[which.sig.down])

    #which.sig <- which(rank(diff$limma.adj.p.value, 
    #  ties.method = "random") <= as.numeric(top))
	}
  
  message("qval threshold for top 300 is: ",
    max(diff$limma.adj.p.value[which.up[which.sig.up]], na.rm=TRUE), " (up)",
    max(diff$limma.adj.p.value[which.down[which.sig.down]], na.rm=TRUE), " (down)")
  print(which.sig)

  # exploratory - merge together those in top
  if (merge){
    which.sig <- which(rank(diff$P.Value, 
      ties.method = "random") <= 10000)
    regions <- makeGRangesFromDataFrame(diff[which.sig,])
    regions$FC <- diff$logFC[which.sig]
    regions$pval <- diff$P.Value[which.sig]
    
    dmrs <- diff[which.sig,
                 grepl("counts", colnames(diff))]
    dmrs <- dmrs[,!grepl("mean", colnames(dmrs))]


    # do separately for up and down
    up <- reduce(regions[regions$FC>0,], 
      min.gapwidth = round(0.1*as.numeric(ws)))
    dn <- reduce(regions[regions$FC<0,], 
      min.gapwidth = round(0.1*as.numeric(ws)))
    rr <- c(up,dn)

    # find overlaps
    ol <- findOverlaps(rr, regions)
    ct <- countOverlaps(rr, regions)  
    # summarize fc
    rr$FC <- by(regions$FC, ol@from, mean)
    rr$pval <- by(regions$pval, ol@from, 
      function(p) pchisq((sum(log(p))*-2), df=length(p)*2, lower.tail=F))
    rr$fdr <- p.adjust(rr$pval)
    rr$nw <- ct
    dmrs <- apply(dmrs, 2, function(x) by(x, ol@from, mean))

    rownames(dmrs) <- paste0(as.character(seqnames(rr)), ":", 
                           start(rr), "-",
                           end(rr))
    #which.top <- which(rank(rr$fdr, ties.method = "random") <= top)

    # enforce half up and half down
    which.up.r <- which(rr$FC > 0)
    which.down.r <- which(rr$FC < 0)
    which.sig.up.r <- which(rank(rr$pval[which.up.r], 
      ties.method = "random") <= as.numeric(top)/2)

    which.sig.down.r <- which(rank(rr$pval[which.down.r], 
      ties.method = "random") <= as.numeric(top)/2)

    which.top <- c(which.up.r[which.sig.up.r], which.down.r[which.sig.down.r])

    dmrs <- dmrs[which.top,]

    if (1 < 2){    
      data.frame(nw=rr$nw[which.top]) %>% ggplot(aes(x = nw)) +
        geom_histogram(bins=100) + 
        xlab("Number of windows")
      ggsave(file.path(out.dir, 
        paste0("nw_", lab1, ".", lab2, "_top", top, ".pdf")),
        width = 4, height =4)
    }
    
  }else{
      message("reformating dmr data.frame")
  	  dmrs <- diff[which.sig,
    	             grepl("counts", colnames(diff))]
    	dmrs <- dmrs[,!grepl("mean", colnames(dmrs))]

      rownames(dmrs) <- paste0(diff[which.sig,]$chr, ":", 
                               diff[which.sig,]$start, "-",
                               diff[which.sig,]$end)
  }
  
  message("simplifying sample names")
  lab1 <- gsub(iter, "", lab1)
  lab2 <- gsub(iter, "", lab2)
	colnames(dmrs)[1:n1] <- paste0(lab1, "_", colnames(dmrs)[1:n1])
	colnames(dmrs)[(n1 + 1):(n1 + n2)] <- paste0(lab2, "_", colnames(dmrs)[(n1 + 1):(n1 + n2)])
	colnames(dmrs) <- gsub(".counts|.rpkm", "", colnames(dmrs))
	colnames(dmrs) <- gsub(".sorted.bam", "", colnames(dmrs))
  colnames(dmrs) <- gsub("\\.1", "", colnames(dmrs))
  colnames(dmrs) <- gsub("_PDAC1", "", colnames(dmrs))

  # normalize 
  message("normalizing...")
  counts.train <- diff[, grepl("counts", colnames(diff))]
  counts.train <- counts.train[, !grepl("mean", colnames(counts.train))]
  which.norm <- which(rowSums(counts.train)>=0.25*ncol(counts.train))
  d_train <- DGEList(counts=counts.train)
  d_train <- calcNormFactors(d_train[which.norm,], refColumn = 1)
  dmrs <- sweep(dmrs, MARGIN=2, 
    (d_train$samples$norm.factors*d_train$samples$lib.size)/1e6, `/`)

	ecolors <- c("white",colorRampPalette( (brewer.pal(9, "Reds")) )(255)) # set colors
	type <- c("lightgrey", "black")
	names(type) <- c(lab1, lab2)
	ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2)),
	                              col = list(Type = type))

  # grab metadata from "master" spreadsheet
  x <- match(gsub(paste0(lab1, "_|", lab2, "_"), "", colnames(dmrs)), 
        master$`Sample number`)
  if (sum(is.na(x))>0)
    message("Warning: can't retrieve meta data for samples ", 
      gsub(paste0(lab1, "_|", lab2, "_"), "", colnames(dmrs))[which(is.na(x))])

  subtype = c(type[!grepl("urineR|rcc", names(type))], "#E69F00", "#56B4E9", "#009E73", "white")
  names(subtype) = c("control", "clear cell", "papillary", "chromophobe", "other")
  st <- as.character(master$Histology[x])
  st[is.na(st)] <- "control"

  bt_col <- c("#9CC6CF", "#8E2043", "#3E5496")
  names(bt_col) <- c("1", "2", "Met")
  bt <- as.character(master$Batch[x])

  inst_col <- c("#59C7EB", "#FEA090", "#0A9086", "#E0607E", "#9AA0A7")
  names(inst_col) <- c("BWH", "Fresh", "DFCI", "MGH", "Sue")
  inst <- as.character(master$Institution[x])
  
  if (length(unique(st))>2){
      ha_column = HeatmapAnnotation(df = 
            data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2), 
                       Subtype = st,
                       Batch = bt,
                       Institution = inst),
                       col = list(Type = type, 
                                  Subtype = subtype, 
                                  Batch = bt_col, 
                                  Institution = inst_col))
  }else{
      ha_column = HeatmapAnnotation(df = 
            data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2),
                       Batch = bt,
                       Institution = inst),
                       col = list(Type = type, Batch = bt_col, Institution = inst_col))
  }

  
  message("building heatmap")
	ht = Heatmap(log(dmrs+1), name = "log(CPM+1)", 
	             top_annotation = ha_column, col = ecolors,
	             show_row_names = FALSE, show_column_names = colnames,
	             column_names_gp = gpar(fontsize = 8),
               column_title = paste0("Top ", top, 
                ifelse(merge, " merged", "")))

  w = 8
  if((n1+n2)>75) w = 12
	pdf(heatmap.file, width=w)
	  draw(ht) 
  dev.off()

  if(training < 1){

   message("extracting test set")
   # overwrite ix1 and ix2 according to what is in diff (in case loading saved)
   grp1 <- gsub(paste0(lab1, "_"), "", 
    colnames(dmrs)[which(grepl(lab1, colnames(dmrs)))])
   grp2 <- gsub(paste0(lab2, "_"), "", 
    colnames(dmrs)[which(grepl(lab2, colnames(dmrs)))])
   ix1 <- match(grp1, gsub("PDAC1_", "", gsub(".sorted.bam", "", sapply(obj1, 
    function(x) x@sample_name))))
   ix1 <- na.omit(ix1)
   ix2 <- match(grp2, gsub("PDAC1_", "", gsub(".sorted.bam", "", sapply(obj2, 
    function(x) x@sample_name))))
   ix2 <- na.omit(ix2)

   m1 <- m2 <- NULL

   # make test heatmap
   if (regulatory){
    if (length(obj1[-ix1])>0)
     m1 <- MEDIPS.meth(MSet1 = obj1[-ix1], CSet = CS, chr = chrs)[unique(ol@from),]
    if (length(obj2[-ix2])>0)
     m2 <- MEDIPS.meth(MSet1 = obj2[-ix2], CSet = CS, chr = chrs)[unique(ol@from),]
   }else{
    if (length(obj1[-ix1])>0)
     m1 <- MEDIPS.meth(MSet1 = obj1[-ix1], CSet = CS, chr = chrs)
    if (length(obj2[-ix2])>0)
     m2 <- MEDIPS.meth(MSet1 = obj2[-ix2], CSet = CS, chr = chrs)
   }

   if (length(obj1[-ix1])>0){
    m1 <- m1[,!grepl("mean", colnames(m1)), drop = FALSE]
    m1 <- m1[,grepl("counts", colnames(m1)), drop = FALSE]
   }
   if (length(obj2[-ix2])>0){
    m2 <- m2[,!grepl("mean", colnames(m2)), drop = FALSE]
    m2 <- m2[,grepl("counts", colnames(m2)), drop = FALSE]
   }

   if (is.null(m1)){
     counts.test <- m2
   }else if (is.null(m2)){
     counts.test <- m1
   }else{
     counts.test <- cbind(m1,m2)
   }
   dmrs_new <- counts.test[which.sig,, drop=FALSE]

   # NORMALIZE - use same reference sample for train and test
   d_all <- DGEList(counts=cbind(counts.train, counts.test))
   d_all <- calcNormFactors(d_all[which.norm,], refColumn = 1)
   
   dmrs_new <- sweep(dmrs_new, MARGIN=2, 
    (d_all$samples$norm.factors[-(1:(n1+n2))]*d_all$samples$lib.size[-(1:(n1+n2))])/1e6, `/`)

  # mat <- cbind(m1[,grepl("rpkm", colnames(m1)), drop = FALSE],
  #              m2[,grepl("rpkm", colnames(m2)), drop = FALSE])
  
   #mat <- mat[rowSums(mat) > 0.2*ncol(mat),]
   #dmrs_new <- sweep(dmrs_new, MARGIN=2, DESeq2::estimateSizeFactorsForMatrix(mat), `/`)

   if (merge){
     dmrs_new <- apply(dmrs_new, 2, function(x) by(x, ol@from, mean))
     rownames(dmrs_new) <- paste0(as.character(seqnames(rr)), ":", 
       start(rr), "-",
       end(rr))

     dmrs_new <- dmrs_new[which.top,]
   }else{
     rownames(dmrs_new) <- paste0(diff[which.sig,]$chr, ":", 
       diff[which.sig,]$start, "-",
       diff[which.sig,]$end)
   }

   tn1 <- length(obj1[-ix1])
   tn2 <- length(obj2[-ix2])
    
   if (tn1 >= 1)
     colnames(dmrs_new)[1:tn1] <- paste0(lab1, "_", colnames(dmrs_new)[1:tn1])
   if (tn2 >= 1)
     colnames(dmrs_new)[(tn1 + 1):(tn1 + tn2)] <- paste0(lab2, "_", colnames(dmrs_new)[(tn1 + 1):(tn1 + tn2)])
   colnames(dmrs_new) <- gsub(".rpkm|.counts", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub(".sorted.bam", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub("\\.1", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub("_PDAC1", "", colnames(dmrs_new))
   rm(m1)
   rm(m2)

   # binary classification
   # glmnet 
   message("building predictive model for training set")
   library(glmnet)
   
   cvob1 = tryCatch({
     cv.glmnet(x=t(log(dmrs+1)),y=c(rep(1, n1), rep(0, n2)),
        family="binomial", alpha=0.5, 
        nfolds = 3, lambda = seq(0.01,0.05,by = 0.01), standardize=FALSE)
   }, error = function(e) {
       tryCatch({
         cv.glmnet(x=t(log(dmrs+1)),y=c(rep(1, n1), rep(0, n2)),
          family="binomial", alpha=0.5, 
          nfolds = 3, lambda = seq(0.01,0.05,by = 0.01), standardize=FALSE)
       }, error = function(e) {
         NULL
       })
   })

   message("evaluating test set")
   # best coefficient
   if (length(cvob1)>0){
     new <- predict(cvob1, newx=t(log(dmrs_new+1)), type = "response", 
      s = "lambda.min")
     
     if (ncol(dmrs_new)>1){
       library(ROCR)
       pred <- prediction(new, as.numeric(grepl(lab1, colnames(dmrs_new))))
       auc <- unlist(performance(pred,"auc")@y.values)
     }else{
      pred <- auc <- NA
     }
   }else{
     auc <- NA
     pred <- new <- rep(NA, ncol(dmrs_new))
   }

   ret_tab <- data.frame(sample_name = colnames(dmrs_new),
    true_label = ifelse(grepl(lab1, colnames(dmrs_new)), 
     lab1, lab2),
   class_prob = as.vector(new),
   auc = auc)


   # mean per group
   mns <- data.frame(grp1 = rowMeans(log(dmrs[,(1:n1)]+1)),
     grp2 = rowMeans(log(dmrs[,(n1+1):(n2+n1)]+1)))
   colnames(mns) <- c(lab1, lab2)
   newdist <- as.matrix(dist(t(cbind(log(dmrs_new+1), mns))))

   ret_tab$class_label <-  rep(NA, ncol(dmrs_new))

   for(col in seq_len(ncol(dmrs_new))){
     ret_tab$class_label[col] <- ifelse(newdist[col,ncol(newdist)] < newdist[col,ncol(newdist)-1],
       lab2, lab1)

     message("New sample ", colnames(dmrs_new)[col], 
       " classified to group ", ret_tab$class_label[col], "(prob= ", 
       ret_tab$class_prob[col], ")")
   }


   ecolors <- c("white",colorRampPalette( (brewer.pal(9, "Reds")) )(255)) # set colors

   type <- c("lightgrey", "black")
   classcol <- c("white", "darkblue")
   probcol <- colorRamp2(c(0, 1), c("darkblue", "white"))
   names(type) <- c(lab1, lab2)
   names(classcol) <- c(lab1, lab2)
   ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs_new)), lab1, lab2), 
     Prob = ret_tab$class_prob, Class = ret_tab$class_label),
     col = list(Type = type, Prob = probcol, Class = classcol))


  # grab metadata from "master" spreadsheet
  x <- match(gsub(paste0(lab1, "_|", lab2, "_"), "", colnames(dmrs)), 
        master$`Sample number`)
  if (sum(is.na(x))>0)
    message("Warning: can't retrieve meta data for samples ", 
      gsub(paste0(lab1, "_|", lab2, "_"), "", colnames(dmrs))[which(is.na(x))])

  subtype = c(type[!grepl("urineR|rcc", names(type))], "#E69F00", "#56B4E9", "#009E73", "white")
  names(subtype) = c("control", "clear cell", "papillary", "chromophobe", "other")
  st <- as.character(master$Histology[x])
  st[is.na(st)] <- "control"

  bt_col <- c("#9CC6CF", "#8E2043", "#3E5496")
  names(bt_col) <- c("1", "2", "Met")
  bt <- as.character(master$Batch[x])

  inst_col <- c("#59C7EB", "#FEA090", "#0A9086", "#E0607E", "#9AA0A7")
  names(inst_col) <- c("BWH", "Fresh", "DFCI", "MGH", "Sue")
  inst <- as.character(master$Institution[x])
  
  if (length(unique(st))>2){
      ha_column = HeatmapAnnotation(df = 
            data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2), 
                       Subtype = st,
                       Batch = bt,
                       Institution = inst,
                       Prob = ret_tab$class_prob,
                       Class = ret_tab$class_label),
                       col = list(Type = type, 
                                  Subtype = subtype, 
                                  Batch = bt_col, 
                                  Institution = inst_col,
                                  Prob = probcol,
                                  Class = classcol))
  }else{
      ha_column = HeatmapAnnotation(df = 
            data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2),
                       Batch = bt,
                       Institution = inst,
                       Prob = ret_tab$class_prob,
                       Class = ret_tab$class_label),
                       col = list(Type = type, Batch = bt_col, 
                        Institution = inst_col,Prob = probcol, Class = classcol))
  }

   ht = Heatmap(log(dmrs_new+1), name = "log(CPM+1)", 
     top_annotation = ha_column, col = ecolors,
     show_row_names = FALSE, show_column_names = colnames,
     column_names_gp = gpar(fontsize = 8),
     column_title = paste0("Top ", top, 
       ifelse(merge, " merged", "")))

   pdf(heatmap.file.test, width = 8)
   draw(ht)
   dev.off()

   if(saveprobs){
     message("saving sample level prob table")
     write.table(ret_tab, quote=FALSE, row.names=FALSE,
      file= prob.table.file, 
      sep = "\t")
    }else{
     message("not saving sample level prob table")
    }

   ret_tab_summary <- data.frame(auc = auc, 
    accuracy = sum(apply(ret_tab, 1, 
     function(x) x[1] == x[4])) / nrow(ret_tab),
    accuracy1 = sum(apply(ret_tab[ret_tab$true_label == lab1,], 1, 
     function(x) x[1] == x[4])) / sum(ret_tab$true_label == lab1),
    accuracy2 = sum(apply(ret_tab[ret_tab$true_label == lab2,], 1, 
     function(x) x[1] == x[4])) / sum(ret_tab$true_label == lab2),
    lab1 = lab1, lab2 = lab2,
    n1 = sum(ret_tab$true_label == lab1), 
    n2 = sum(ret_tab$true_label == lab2))
  
   return(ret_tab_summary)
  }

  if(!is.null(add_samp)){

    m1 <- MEDIPS.meth(MSet1 = add_samp, CSet = CS, chr = chrs)[which.sig,]

    m1 <- m1[,grepl("rpkm", colnames(m1)), drop = FALSE]
    m1 <- m1[,!grepl("mean", colnames(m1)), drop = FALSE]
    colnames(m1) <- paste0(lab_add, "_", colnames(m1)[1:ncol(m1)])

    dmrs <- cbind(dmrs,m1)
    colnames(dmrs) <- gsub("_HY2FCBBXX_L005","", colnames(dmrs))
    colnames(dmrs) <- gsub(".rpkm","", colnames(dmrs))
    colnames(dmrs) <- gsub(".sorted.bam", "", colnames(dmrs))
    colnames(dmrs) <- gsub("\\.1", "", colnames(dmrs))
    colnames(dmrs) <- gsub("_PDAC1", "", colnames(dmrs))
    ecolors <- c("white",colorRampPalette( (brewer.pal(9, "Reds")) )(255)) # set colors
    type <- c("lightgrey", "black", "blue")
    names(type) <- c(lab1, lab2, lab_add)
    ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, ifelse(grepl(lab2, colnames(dmrs)) ,lab2, lab_add))),
     col = list(Type = type))
    ht = Heatmap(log(dmrs+1), name = "log(RPKM+1)", 
      top_annotation = ha_column, col = ecolors,
      show_row_names = FALSE, show_column_names = colnames,
      column_names_gp = gpar(fontsize = 8))

    pdf(heatmap.file.add, width = 8)
    draw(ht)
    dev.off()

  }

  if (!is.null(validate1) || !is.null(validate2)){
  
   if (!is.null(validate1)){
     m1 <- MEDIPS.meth(MSet1 = validate1, CSet = CS, chr = chrs)
     m1 <- m1[,grepl("rpkm", colnames(m1)), drop = FALSE]
     m1 <- m1[,!grepl("mean", colnames(m1)), drop = FALSE]
   }else{
     m1 <- data.frame()
   }

   if (!is.null(validate2)){
     m2 <- MEDIPS.meth(MSet1 = validate2, CSet = CS, chr = chrs)
     m2 <- m2[,grepl("rpkm", colnames(m2)), drop = FALSE]
     m2 <- m2[,!grepl("mean", colnames(m2)), drop = FALSE]
   }else{
     m2 <- data.frame()
   }
  
   if (nrow(m1)==0){
    m <- m2
   }else if (nrow(m2)==0){
    m <- m1
   }else{
    m <- cbind(m1,m2)
   }

   # NORMALIZE - use same reference sample for train and test
   d_all <- DGEList(counts=cbind(counts.train, m))
   d_all <- calcNormFactors(d_all[which.norm,], refColumn = 1)
   
   dmrs_new <- sweep(m[which.sig,], MARGIN=2, 
    (d_all$samples$norm.factors[-(1:(n1+n2))]*d_all$samples$lib.size[-(1:(n1+n2))])/1e6, `/`)

   if (merge){
    dmrs_new <- apply(dmrs_new, 2, function(x) by(x, ol@from, mean))
    rownames(dmrs_new) <- paste0(as.character(seqnames(rr)), ":", 
     start(rr), "-",
     end(rr))
    dmrs_new <- dmrs_new[which.top,]
   }else{
    rownames(dmrs_new) <- paste0(diff[which.sig,]$chr, ":", 
     diff[which.sig,]$start, "-",
     diff[which.sig,]$end)
   }

   if (nrow(m1)>0)
     colnames(dmrs_new)[1:ncol(m1)] <- paste0(lab1, "_", colnames(dmrs_new)[1:ncol(m1)])
   if (nrow(m2)>0)
     colnames(dmrs_new)[(ncol(m1) + 1):(ncol(m1) + ncol(m2))] <- paste0(gsub("_partial", "",lab2), "_", colnames(dmrs_new)[(ncol(m1) + 1):(ncol(m1) + ncol(m2))])

   colnames(dmrs_new) <- gsub(".rpkm", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub(".sorted.bam", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub("\\.1", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub("_PDAC1", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub("_partial", "", colnames(dmrs_new))
   rm(m1)
   rm(m2)

   # binary classification
   # glmnet 

   library(glmnet)
   cvob1=cv.glmnet(x=t(log(dmrs+1)),y=as.numeric(grepl(lab1, colnames(dmrs))),
     family="binomial", alpha=0.5, 
     nfolds = 3, lambda = seq(0.01,0.05,by = 0.01), standardize=FALSE)

   # best coefficient
   new <- predict(cvob1, newx=t(dmrs_new), type = "response", s = "lambda.min")
  
   library(ROCR)
   fact <- as.factor(as.numeric(grepl(lab1, colnames(dmrs_new))))
   if (length(levels(fact)) == 1)
     levels(fact) <- c(1,0)
   pred <- prediction(new, fact)

   auc = tryCatch({
     unlist(performance(pred,"auc")@y.values)
   }, error = function(e) {
     NA
   })

   if(!is.na(auc)){
     pdf(file.path(out.dir, paste0("ROC_validation_", validate_lab, ".pdf")), 
       width =4, height=4)
     rocvals <- performance(pred, "tpr", "fpr")
     plot(rocvals); abline(0,1,lty=2,col="grey")
     dev.off()
   }

   ret_tab <- data.frame(sample_name = colnames(dmrs_new),
     true_label = ifelse(grepl(gsub("_partial","",lab1), colnames(dmrs_new)), 
     gsub("_partial","",lab1), gsub("_partial", "",lab2)),
   class_prob = as.vector(new),
   auc = auc)


   # mean per group
   mns <- data.frame(grp1 = rowMeans(log(dmrs[,grepl(lab1, colnames(dmrs))]+1)),
     grp2 = rowMeans(log(dmrs[,grepl(gsub("_partial", "",lab2), colnames(dmrs))]+1)))
   colnames(mns) <- c(lab1, gsub("_partial", "",lab2))
   newdist <- as.matrix(dist(t(cbind(log(dmrs_new+1), mns))))

   ret_tab$class_label <-  rep(NA, ncol(dmrs_new))

   for(col in seq_len(ncol(dmrs_new))){
     ret_tab$class_label[col] <- ifelse(newdist[col,ncol(newdist)] < newdist[col,ncol(newdist)-1],
       gsub("_partial", "",lab2), lab1)

     message("New sample ", colnames(dmrs_new)[col], 
       " classified to group ", ret_tab$class_label[col], "(prob= ", 
       ret_tab$class_prob[col], ")")
   }

   if(saveprobs){
     message("saving sample level prob table")
     write.table(ret_tab%>% select(-class_label), quote=FALSE, row.names=FALSE,
      file=file.path(out.dir, paste0("sampleprob_table_", validate_lab, 
      "_top", top, ".txt")), 
      sep = "\t")
    }else{
     message("not saving sample level prob table")
    }

  ret_tab <- mutate(ret_tab, true_label = ifelse(true_label == "rcc", 
    "RCC", "Control"))
  cols <- c("RCC" = "#56B4E9", "Control" = "#E69F00")
  if(grepl("met", validate_lab)){
    cols <- c("RCC Met" = "#56B4E9", "Control" = "#E69F00")
    ret_tab <- mutate(ret_tab, true_label = ifelse(true_label == "RCC", 
    "RCC Met", "Control"))
  }

  ret_tab %>% ggplot(aes(x=sample_name, y=class_prob, 
    fill = true_label)) + 
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size =7))+
    ylab("Probability RCC") +
    xlab("Sample") +
    labs(fill="True Class") +
    scale_fill_manual(values = cols) 
  ggsave(file.path(out.dir, paste0("sampleprob_barplot_", validate_lab, 
      "_top", top, ".pdf")), height=2.5, width=4.75)

   ecolors <- c("white",colorRampPalette( (brewer.pal(9, "Reds")) )(255)) # set colors

   type <- c("lightgrey", "black")
   classcol <- c("white", "darkblue")
   probcol <- colorRamp2(c(0, 1), c("darkblue", "white"))
   names(type) <- c(lab1, gsub("_partial", "",lab2))
   names(classcol) <- c(lab1, gsub("_partial", "",lab2))

   # grab metadata from "master" spreadsheet
  x <- match(gsub(paste0(lab1, "_|", lab2, "_"), "", colnames(dmrs)), 
        master$`Sample number`)
  if (sum(is.na(x))>0)
    message("Warning: can't retrieve meta data for samples ", 
      gsub(paste0(lab1, "_|", lab2, "_"), "", colnames(dmrs))[which(is.na(x))])

  subtype = c(type[!grepl("urineR|rcc", names(type))], "#E69F00", "#56B4E9", "#009E73", "white")
  names(subtype) = c("control", "clear cell", "papillary", "chromophobe", "other")
  st <- as.character(master$Histology[x])
  st[is.na(st)] <- "control"

  bt_col <- c("#9CC6CF", "#8E2043", "#3E5496")
  names(bt_col) <- c("1", "2", "Met")
  bt <- as.character(master$Batch[x])

  inst_col <- c("#59C7EB", "#FEA090", "#0A9086", "#E0607E", "#9AA0A7")
  names(inst_col) <- c("BWH", "Fresh", "DFCI", "MGH", "Sue")
  inst <- as.character(master$Institution[x])



  ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs_new)), lab1, gsub("_partial", "",lab2)), 
     Prob = ret_tab$class_prob, Class = ret_tab$class_label),
     col = list(Type = type, Prob = probcol, Class = classcol))
  
  if (length(unique(st))>2){
      ha_column = HeatmapAnnotation(df = 
            data.frame(Type = ifelse(grepl(lab1, colnames(dmrs_new)), lab1, lab2), 
                       Subtype = st,
                       Batch = bt,
                       Institution = inst,
                       Prob = ret_tab$class_prob,
                       Class = ret_tab$class_label),
                       col = list(Type = type, 
                                  Subtype = subtype, 
                                  Batch = bt_col, 
                                  Institution = inst_col,
                                  Prob = probcol,
                                  Class = classcol))
  }else{
      ha_column = HeatmapAnnotation(df = 
            data.frame(Type = ifelse(grepl(lab1, colnames(dmrs_new)), lab1, lab2),
                       Batch = bt,
                       Institution = inst,
                       Prob = ret_tab$class_prob,
                       Class = ret_tab$class_label),
                       col = list(Type = type, Batch = bt_col, 
                        Institution = inst_col,Prob = probcol, Class = classcol))
  }

   ht = Heatmap(log(dmrs_new+1), name = "log(RPKM+1)", 
     top_annotation = ha_column, col = ecolors,
     show_row_names = FALSE, show_column_names = colnames,
     column_names_gp = gpar(fontsize = 8),
     column_title = paste0("Top ", top, 
       ifelse(merge, " merged", "")))

  pdf(heatmap.file.test, width = 8)
   draw(ht)
  dev.off()

   ret_tab_summary <- data.frame(auc = auc, 
    accuracy = sum(apply(ret_tab, 1, 
     function(x) x[1] == x[4])) / nrow(ret_tab),
    accuracy1 = sum(apply(ret_tab[ret_tab$true_label == lab1,], 1, 
     function(x) x[1] == x[4])) / sum(ret_tab$true_label == lab1),
    accuracy2 = sum(apply(ret_tab[ret_tab$true_label == gsub("_partial", "",lab2),], 1, 
     function(x) x[1] == x[4])) / sum(ret_tab$true_label == gsub("_partial", "",lab2)),
    lab1 = lab1, lab2 = gsub("_partial", "",lab2),
    n1 = sum(ret_tab$true_label == lab1), 
    n2 = sum(ret_tab$true_label == gsub("_partial", "",lab2)))
  
   return(ret_tab_summary)
  }
}



if (iter == 105){
  if(FALSE){
  compute.diff(obj1 = medip.rcc, obj2 = medip.control,
	         lab1 = "rcc", lab2 = "control",
           out.dir = file.path(outdir, ".."), top = ntop)

  compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
	         lab1 = "urineR", lab2 = "urineC",
           out.dir = file.path(outdir, ".."), top = ntop)

  compute.diff(obj1 = medip.jan2020[m2$Source=="Plasma" & m2$Status == "RCC"],
           obj2 = medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"],
           lab1 = "rcc_New", lab2 = "control_New",
           out.dir = file.path(outdir, ".."), top = ntop)

  compute.diff(obj1 = medip.jan2020[m2$Source=="Urine" & m2$Status == "RCC"], 
           obj2 = medip.jan2020[m2$Source=="Urine" & m2$Status == "Control"],
           lab1 = "urineR_New", lab2 = "urineC_New",
           out.dir = file.path(outdir, ".."), top = ntop)

  ## train on orig, predict on new
  compute.diff(obj1 = medip.rcc, obj2 = medip.control,
           lab1 = "rcc", lab2 = "control",
           validate1 = medip.jan2020[m2$Source=="Plasma" & m2$Status == "RCC"],
           validate2 = medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"],
           out.dir = file.path(outdir, ".."), top = ntop,
           validate_lab = "validate_plasma")

 compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
           lab1 = "urineR", lab2 = "urineC",
           validate1 = medip.jan2020[m2$Source=="Urine" & m2$Status == "RCC"], 
           validate2 = medip.jan2020[m2$Source=="Urine" & m2$Status == "Control"],
           out.dir = file.path(outdir, ".."), top = ntop,
           validate_lab = "validate_urine")


  # rcc metastatic (train on batch 1 RCC/control, predict RCCmet/control batch2)
  met <- compute.diff(obj1 = medip.rcc, obj2 = medip.control,
           lab1 = "rcc", lab2 = "control",
           top = ntop, out.dir = file.path(outdir, ".."),
           validate1 = medip.rcc_M, 
           validate2 = medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"],
           validate_lab = "vRCCmet",
           saveprobs = TRUE)
  met$iteration <- NA
  met$window <- ws
  met$method <- c("original")
  
  write.table(met, quote=FALSE, row.names=FALSE,
  file=file.path(out.dir = file.path(outdir, ".."), 
    paste0("validation_accuracy_RCCmet_table_top", ntop, ".txt")), 
    sep = "\t")


  # repeat prev but with random (mixed) control group
  set.seed(125)
  y1 <- sample(1:sum(m2$Source=="Plasma" & m2$Status == "Control"), 
    ceiling(sum(m2$Source=="Plasma" & m2$Status == "Control")/2))
  y2 <- sample(1:length(medip.control), 
    ceiling(length(medip.control)/2))
  met <- compute.diff(obj1 = medip.rcc, 
           obj2 = c(medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"][y1],
            medip.control[y2]),
           lab1 = "rcc", lab2 = "controlMix",
           top = ntop, out.dir = file.path(outdir, ".."),
           validate1 = medip.rcc_M, 
           validate2 =  c(medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"][-y1],
            medip.control[-y2]),
           validate_lab = "vRCCmet_mix",
           saveprobs = TRUE)
  met$iteration <- NA
  met$window <- ws
  met$method <- c("original")
  
  write.table(met, quote=FALSE, row.names=FALSE,
  file=file.path(out.dir = file.path(outdir, ".."), 
    paste0("validation_accuracy_RCCmet_table_top", ntop, ".txt")), 
    sep = "\t")
  }


  # repeat prev but with FULL control group
  dir.create(file.path(outdir, "../pooled_c"))

  met <- compute.diff(obj1 = medip.rcc, 
           obj2 = c(medip.control,medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"]),
           lab1 = "rcc", lab2 = "control",
           top = ntop, out.dir = file.path(outdir, "../pooled_c"))

}


# iterative RCC met control splits (75% of both control batches combined for training)

if(iter <= 100){

  outdir_iterm <- paste0("/scratch/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws, "/iter_", 
    str_pad(iter, 3, pad = "0"))
  dir.create(outdir_iterm, showWarnings = FALSE)

  set.seed(iter*20200131)
  y1 <- sample(1:sum(m2$Source=="Plasma" & m2$Status == "Control"), 
    round(0.75*sum(m2$Source=="Plasma" & m2$Status == "Control")))
  y2 <- sample(1:length(medip.control), round(0.75*length(medip.control)))
  met <- compute.diff(obj1 = medip.rcc, 
           obj2 = c(medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"][y1],
            medip.control[y2]),
           lab1 = "rcc", lab2 = "control",
           top = ntop, 
           validate1 = medip.rcc_M, 
           validate2 =  c(medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"][-y1],
            medip.control[-y2]),
           validate_lab = "vRCCmet",
           out.dir = file.path(outdir_iterm),
           saveprobs = TRUE)
  met$iteration <- NA
  met$window <- ws
  met$method <- c("original")
  
  write.table(met, quote=FALSE, row.names=FALSE,
  file=file.path(out.dir = file.path(outdir, ".."), 
    paste0("validation_accuracy_RCCmet_table_top", ntop, "_iter", iter, ".txt")), 
    sep = "\t")

}

# create pooled set

# joint metadata for RCC samps
meta <- rbind(data.frame(`Sample number`=meta$`Sample number`,
                   Histology=meta$Histology),
        data.frame(`Sample number`=meta2$ID,
                   Histology=tolower(meta2$Histology)))
colnames(meta)[1] <- "Sample number"

# remove three RCCMet samples since don't have histology
metids <- gsub(".sorted.bam", "", sapply(medip.rcc_M, function(x) x@sample_name))
x <- match( metids, master$ID )
excl <- master$ID[x][grepl("Exclude", master$Inclusion[x])]
if(length(excl) > 0)
  medip.rcc_M <- medip.rcc_M[-which(metids %in% excl)]

# joint medips objs

medip.rcc <- c(medip.rcc, medip.jan2020[m2$Source=="Plasma" & m2$Status == "RCC"], medip.rcc_M)
medip.control <- c(medip.control, medip.jan2020[m2$Source=="Plasma" & m2$Status == "Control"])
medip.urineR <- c(medip.urineR, medip.jan2020[m2$Source=="Urine" & m2$Status == "RCC"])
medip.urineC <- c(medip.urineC, medip.jan2020[m2$Source=="Urine" & m2$Status == "Control"])


if (iter == 99){

dir.create(file.path(outdir, "../pooled_m"))
dir.create(file.path(outdir, "../pooled"))

compute.diff(obj1 = medip.rcc, obj2 = medip.control,
           lab1 = "rcc", lab2 = "control",
           out.dir = file.path(outdir, "../pooled_m"), top = ntop)

compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
           lab1 = "urineR", lab2 = "urineC",
           out.dir = file.path(outdir, "../pooled"), top = ntop)

}

# leave-one-out

# rcc plasma
iter
length(medip.rcc)
if(as.numeric(iter) <= length(medip.rcc)){
  names(iter) <- "obj1"
  compute.diff(obj1 = medip.rcc, obj2 = medip.control,
             holdout = iter,
             lab1 = paste0("rcc",iter), lab2 = "control",
             out.dir = file.path(outdir_m), top = ntop)
}

# control plasma
iter
length(medip.control)
if(as.numeric(iter) <= length(medip.control)){
  names(iter) <- "obj2"
  compute.diff(obj1 = medip.rcc, obj2 = medip.control,
             holdout = iter,
             lab1 = "rcc", lab2 = paste0("control",iter),
             out.dir = file.path(outdir_m), top = ntop)
}


# rcc urine
iter
length(medip.urineR)
if(as.numeric(iter) <= length(medip.urineR)){
  names(iter) <- "obj1"
  compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
             holdout = iter,
             lab1 = paste0("urineR",iter), lab2 = "urineC",
             out.dir = file.path(outdir), top = ntop)
}

# control urine
iter
length(medip.urineC)
if(as.numeric(iter) <= length(medip.urineC)){
  names(iter) <- "obj2"
  compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
             holdout = iter,
             lab1 = "urineR", lab2 = paste0("urineC",iter),
             out.dir = file.path(outdir), top = ntop)
}

