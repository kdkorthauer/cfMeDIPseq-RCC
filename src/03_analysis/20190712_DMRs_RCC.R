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
source("20190722_limmamedips.R")

# get windowsize
ws <- Sys.getenv("WINDOWSIZE")
iter <- Sys.getenv("iter")
ntop <- Sys.getenv("ntop")

# dir where binned medips objects of all samples are saved
medipdir <- paste0("../../out/MEDIPS_", ws)
dir.create(medipdir, showWarnings = FALSE)

# outdir to iteration specific for all output
outdir <- paste0("../../out/MEDIPS_", ws, "/iteration_", str_pad(iter, 3, pad = "0"))
dir.create(outdir, showWarnings = FALSE)

bamdir <- "../../out/sortedbam_dup"

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

# exclude RCC samples (per Sandor email)
meta.rcc <- read_excel("../../data/RCC/42 RCC fastq files.xlsx")
ids.rcc <- gsub(".sorted.bam", "", 
              gsub("../../out/sortedbam_dup/RCC/", "", bam.rcc))
bam.rcc.extra <- bam.rcc[!(ids.rcc %in% meta.rcc$`Fastq name`)]
bam.rcc <- bam.rcc[ids.rcc %in% meta.rcc$`Fastq name`]


## bladder cancer data
bam.blca <- list.files(file.path(bamdir, "BLCA"), "*.bam", 
  full.names = TRUE)
bam.blca <- bam.blca[!grepl(".bai", bam.blca)]


set.seed(3874*as.numeric(iter))

# canonical chrs
chr.select <- paste0("chr", c(1:22, "X", "Y", "M"))
#chr.select <- "chr21"

# read in medip objects (created in 20190712_MEDIPS_RCC.R)
medip.rcc <- readRDS(file.path(medipdir, "medip.rcc.rds"))
medip.control <- readRDS(file.path(medipdir, "medip.control.rds"))
medip.urineR <- readRDS(file.path(medipdir, "medip.urineR.rds"))
medip.urineC <- readRDS(file.path(medipdir, "medip.urineC.rds"))

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

meta <- read_excel("../../data/RCC/Keegan - RCC plasma.xlsx", 
  sheet = 1)
CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip.rcc[[1]])

# summary plots - make for all sets
# for comparison 
if (iter == 1){

# saturation analysis
if (!file.exists(file.path(medipdir, "saturation_rcc_filtered.pdf"))){
  pdf(file.path(medipdir, "saturation_rcc_filtered.pdf"), width = 10, height = 10) 
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

if (!file.exists(file.path(medipdir, "saturation_control_filtered.pdf"))){
  pdf(file.path(medipdir, "saturation_control_filtered.pdf"), width = 10, height = 10) 
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


if (!file.exists(file.path(medipdir, "saturation_urineC_filtered.pdf"))){
  pdf(file.path(medipdir, "saturation_urineC_filtered.pdf"), width = 10, height = 10) 
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

if (!file.exists(file.path(medipdir, "saturation_urineR_filtered.pdf"))){
  pdf(file.path(medipdir, "saturation_urineR_filtered.pdf"), width = 10, height = 10) 
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
df <- rbind(noCpGcov(medip.rcc, CS=CS, "rcc"),
	noCpGcov(medip.control, CS=CS, "ctrl"),
	noCpGcov(medip.urineR, CS=CS, "urineR"),
	noCpGcov(medip.urineC, CS=CS, "urineC"))
ggplot(df, aes(x = pctNoCpG*100)) +
  geom_histogram() +
  facet_wrap(~type) +
  xlab("Percent of total reads mapping to bins with 0 CpGs")
ggsave(file.path(medipdir, "PctReadsinBinswithZeroCpGs.pdf"))

ggplot(df, aes(x = pctNonzeroCov*100)) +
  geom_histogram() + 
  xlab("Percent of bins with nonzero coverage")+
  facet_wrap(~type) 
ggsave(file.path(medipdir, "PctBinswithNonzeroCov.pdf"))


df <- rbind(noCpGcov(medip.rcc_D, CS=CS, "rcc"),
  noCpGcov(medip.control_D, CS=CS, "ctrl"))
ggplot(df, aes(x = pctNoCpG*100)) +
  geom_histogram() +
  facet_wrap(~type) +
  xlab("Percent of total reads mapping to bins with 0 CpGs")
ggsave(file.path(medipdir, "PctReadsinBinswithZeroCpGs_D.pdf"))

ggplot(df, aes(x = pctNonzeroCov*100)) +
  geom_histogram() + 
  xlab("Percent of bins with nonzero coverage")+
  facet_wrap(~type) 
ggsave(file.path(medipdir, "PctBinswithNonzeroCov_D.pdf"))

# explore sample depths
depths <- function(mdobjlist, CS, type){
  depth <- data.frame(sapply(mdobjlist, function(x){
  	  x@genome_count[CS@genome_CF > 0]
  }))
  colnames(depth) <- sapply(mdobjlist, function(x){ 
  	gsub(".sorted.bam", "", x@sample_name)})
  depth <- depth %>% 
    filter(rowSums(depth == 0) < ncol(depth)) %>%
    tidyr::gather("sample", "count") 

  depth %>% filter(count > 0) %>%
  ggplot(aes(x = count)) +
    geom_histogram() +
    facet_wrap(~sample) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Coverage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(file.path(medipdir, paste0("NonzeroCovDistBySample_", type, ".pdf")), 
  	height = 9, width = 9)
}

depths(medip.rcc, CS, "rcc")
depths(medip.control, CS, "ctrl")
depths(medip.urineR, CS, "urineR")
depths(medip.urineC, CS, "urineC")

depths(medip.rcc_D, CS, "rcc_D")
depths(medip.control_D, CS, "control_D")

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

if (!file.exists(file.path(medipdir, "PCA_2_3_cf_rcc.pdf"))){
  ### all
  df <- cbind(depths(medip.rcc, CS, "rcc"),
  	depths(medip.control, CS, "ctrl"),
  	depths(medip.urineR, CS, "urineR"),
  	depths(medip.urineC, CS, "urineC"))
  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  sf <- DESeq2::estimateSizeFactorsForMatrix(df)
  df <- sweep(df, MARGIN=2, sf, `/`)

  grp <- gsub("_.*", "", colnames(df))
  ids = unlist(sapply(strsplit(colnames(df), "_"), function(x) x[[2]]))
  x <- match(ids, meta$`Sample number`)
  subtype <- ifelse(grp == "rcc" & ids %in% meta$`Sample number`,
    meta$Histology[x], NA)

 
  pcs <- Morpho::prcompfast(t(log(df+1)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) 
  
  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = type), size = 2)  
  ggsave(file.path(medipdir, "PCA_1_2_cf.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = type)) +
    geom_point(size = 2)
  ggsave(file.path(medipdir, "PCA_2_3_cf.pdf"), width = 4.5, height = 3.5)
  
  # dendrogram
  dm <- dist(t(log(as.matrix(df) + 1)))
  ggdendrogram(hclust(dm))
  hc <- hclust(dm)
  dendr <- dendro_data(hc, type="rectangle") 
  clust.gr<-data.frame(label = colnames(df), type=grp)
  text.df<-merge(label(dendr),clust.gr,by="label")

  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_text(data=text.df, aes(x=x, y=y, label=label, hjust=0,color=type), 
      size=2) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
    theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())
  ggsave(file.path(medipdir, "hclust_cf.pdf"), width = 7.5, height = 6.5)

 #### no urine
 df <- cbind(depths(medip.rcc, CS, "rcc"),
    depths(medip.control, CS, "ctrl"))
  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  df <- sweep(df, MARGIN=2, DESeq2::estimateSizeFactorsForMatrix(df), `/`)

  grp <- gsub("_.*", "", colnames(df))
  ids = unlist(sapply(strsplit(colnames(df), "_"), function(x) x[[2]]))
 
  pcs <- Morpho::prcompfast(t(log(df+1)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           id = ids) 
  
  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = type),size = 2)  
  ggsave(file.path(medipdir, "PCA_1_2_cf_plasma.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = type)) +
    geom_point(size = 2)
  ggsave(file.path(medipdir, "PCA_2_3_cf_plasma.pdf"), width = 4.5, height = 3.5)

  # dendrogram
  dm <- dist(t(log(as.matrix(df) + 1)))
  ggdendrogram(hclust(dm))
  hc <- hclust(dm)
  dendr <- dendro_data(hc, type="rectangle") 
  clust.gr<-data.frame(label = colnames(df), type=grp)
  text.df<-merge(label(dendr),clust.gr,by="label")

  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_text(data=text.df, aes(x=x, y=y, label=label, hjust=0,color=type), 
      size=2) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
    theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())
  ggsave(file.path(medipdir, "hclust_plasma.pdf"), width = 5.5, height = 4.5)


  # just urine
  df <- cbind(depths(medip.urineR, CS, "urineR"),
    depths(medip.urineC, CS, "urineC"))
  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  df <- sweep(df, MARGIN=2, DESeq2::estimateSizeFactorsForMatrix(df), `/`)

  grp <- gsub("_.*", "", colnames(df))
  ids = unlist(sapply(strsplit(colnames(df), "_"), function(x) x[[2]]))
 
  pcs <- Morpho::prcompfast(t(log(df+1)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           id = ids) 
  
  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = type),size = 2)  
  ggsave(file.path(medipdir, "PCA_1_2_cf_urine.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = type)) +
    geom_point(size = 2)
  ggsave(file.path(medipdir, "PCA_2_3_cf_urine.pdf"), width = 4.5, height = 3.5)
 
  # dendrogram
  dm <- dist(t(log(as.matrix(df) + 1)))
  ggdendrogram(hclust(dm))
  hc <- hclust(dm)
  dendr <- dendro_data(hc, type="rectangle") 
  clust.gr<-data.frame(label = colnames(df), type=grp)
  text.df<-merge(label(dendr),clust.gr,by="label")

  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_text(data=text.df, aes(x=x, y=y, label=label, hjust=0,color=type), 
      size=2) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
    theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())
  ggsave(file.path(medipdir, "hclust_urine.pdf"), width = 5.5, height = 4.5)


  # only rcc
  df <- cbind(depths(medip.rcc, CS, "rcc"),
    depths(medip.control, CS, "control"))
  df <- df %>%
    filter(rowMeans(df > 0) > 0.20)
  df <- df %>%
    filter(rowMeans(df) > 1) 

  # normalize
  df <- sweep(df, MARGIN=2, DESeq2::estimateSizeFactorsForMatrix(df), `/`)

  grp <- gsub("_.*", "", colnames(df))
  ids = unlist(sapply(strsplit(colnames(df), "_"), function(x) x[[2]]))
  x <- match(ids, meta$`Sample number`)
  subtype <- ifelse(grp == "rcc" & ids %in% meta$`Sample number`,
    meta$Histology[x], NA)

  pcs <- Morpho::prcompfast(t(log(df+1)), 
    center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) 
  tidydf <- mutate(tidydf, 
    subtype = ifelse(!is.na(subtype), subtype, "control"))
  
  ggplot() +
    geom_point(data=tidydf, 
      aes(x = PC1, y = PC2, color = subtype), size = 2)
  ggsave(file.path(medipdir, "PCA_1_2_cf_rcc.pdf"), width = 4.5, height = 3.5)

  ggplot() +
    geom_point(data=tidydf, 
      aes(x = PC2, y = PC3, color = subtype), size = 2)
  ggsave(file.path(medipdir, "PCA_2_3_cf_rcc.pdf"), width = 4.5, height = 3.5)
  
  ### hclust
  dm <- dist(t(log(as.matrix(df) + 1)))
  subtype <- ifelse(!is.na(subtype), subtype, "control")
  ggdendrogram(hclust(dm))
  hc <- hclust(dm)
  dendr <- dendro_data(hc, type="rectangle") 
  clust.gr<-data.frame(label = colnames(df), type=subtype)
  text.df<-merge(label(dendr),clust.gr,by="label")

  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_text(data=text.df, aes(x=x, y=y, label=label, hjust=0,color=type), 
      size=2) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
    theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())
  ggsave(file.path(medipdir, "hclust_rcc.pdf"), width = 5.5, height = 4.5)

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
                       saveprobs=FALSE){

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
    }
    ix1 <- sample(1:length(obj1), n1)
    ix2 <- sample(1:length(obj2), n2)
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
	if (!is.null(sig.level)){
	 	which.sig <- which(diff$limma.adj.p.value < sig.level & 
		               abs(diff$logFC) > 2 &
		               !is.na(diff$P.Value) )
	}else{# top ntop - half up and half down
    which.up <- which(diff$logFC > 0)
    which.down <- which(diff$logFC < 0)

		which.sig.up <- which(rank(diff$limma.adj.p.value[which.up], 
      ties.method = "random") <= as.numeric(top)/2)

    which.sig.down <- which(rank(diff$limma.adj.p.value[which.down], 
      ties.method = "random") <= as.numeric(top)/2)

    which.sig <- c(which.up[which.sig.up], which.down[which.sig.down])

    #which.sig <- which(rank(diff$limma.adj.p.value, 
    #  ties.method = "random") <= as.numeric(top))
	}

  # exploratory - merge together those in top
  if (merge){
    which.sig <- which(rank(diff$limma.adj.p.value, 
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
    rr$fdr <- p.adjust(by(regions$pval, ol@from, 
      function(p) pchisq((sum(log(p))*-2), df=length(p)*2, lower.tail=F)))
    rr$nw <- ct
    dmrs <- apply(dmrs, 2, function(x) by(x, ol@from, mean))

    rownames(dmrs) <- paste0(as.character(seqnames(rr)), ":", 
                           start(rr), "-",
                           end(rr))
    #which.top <- which(rank(rr$fdr, ties.method = "random") <= top)

    # enforce half up and half down
    which.up.r <- which(rr$FC > 0)
    which.down.r <- which(rr$FC < 0)
    which.sig.up.r <- which(rank(rr$fdr[which.up.r], 
      ties.method = "random") <= as.numeric(top)/2)

    which.sig.down.r <- which(rank(rr$fdr[which.down.r], 
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
  	  dmrs <- diff[which.sig,
    	             grepl("counts", colnames(diff))]
    	dmrs <- dmrs[,!grepl("mean", colnames(dmrs))]

      rownames(dmrs) <- paste0(diff[which.sig,]$chr, ":", 
                               diff[which.sig,]$start, "-",
                               diff[which.sig,]$end)
  }


	colnames(dmrs)[1:n1] <- paste0(lab1, "_", colnames(dmrs)[1:n1])
	colnames(dmrs)[(n1 + 1):(n1 + n2)] <- paste0(lab2, "_", colnames(dmrs)[(n1 + 1):(n1 + n2)])
	colnames(dmrs) <- gsub(".counts|.rpkm", "", colnames(dmrs))
	colnames(dmrs) <- gsub(".sorted.bam", "", colnames(dmrs))
  colnames(dmrs) <- gsub("\\.1", "", colnames(dmrs))
  colnames(dmrs) <- gsub("_PDAC1", "", colnames(dmrs))

  # normalize 
  counts.train <- diff[, grepl("counts", colnames(diff))]
  counts.train <- counts.train[, !grepl("mean", colnames(counts.train))]
  which.norm <- which(rowSums(counts.train)>=10)
  d_train <- DGEList(counts=counts.train)
  d_train <- calcNormFactors(d_train[which.norm,], refColumn = 1)
  dmrs <- sweep(dmrs, MARGIN=2, 
    (d_train$samples$norm.factors*d_train$samples$lib.size)/1e6, `/`)

	ecolors <- c("white",colorRampPalette( (brewer.pal(9, "Reds")) )(255)) # set colors
	type <- c("lightgrey", "black")
	names(type) <- c(lab1, lab2)
	ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2)),
	                              col = list(Type = type))

  if (lab1 == "rcc" || lab2 == "rcc"){
      rccsamps <- which(grepl("rcc", colnames(dmrs)))
      x <- match(gsub("rcc_", "", colnames(dmrs)[rccsamps]), 
        meta$`Sample number`)
      subtype = c(type[names(type) != "rcc"], "#E69F00", "#56B4E9", "#009E73")
      names(subtype) = c("control", "clear cell", "papillary", "chromophobe")
      st <- colnames(dmrs)
      st[rccsamps] <- meta$Histology[x]
      st[-rccsamps] <- "control"

      ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2), 
        Subtype = st),
                                col = list(Type = type, Subtype = subtype))
  }

	ht = Heatmap(log(dmrs+1), name = "log(CPM+1)", 
	             top_annotation = ha_column, col = ecolors,
	             show_row_names = FALSE, show_column_names = colnames,
	             column_names_gp = gpar(fontsize = 9),
               column_title = paste0("Top ", top, 
                ifelse(merge, " merged", "")))

	pdf(heatmap.file, width=8)
	  draw(ht)
  dev.off()

  if(training < 1){
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

   # make test heatmap
   if (regulatory){
     m1 <- MEDIPS.meth(MSet1 = obj1[-ix1])[unique(ol@from),]
     m2 <- MEDIPS.meth(MSet1 = obj2[-ix2])[unique(ol@from),]
   }else{
     m1 <- MEDIPS.meth(MSet1 = obj1[-ix1])
     m2 <- MEDIPS.meth(MSet1 = obj2[-ix2])
   }

   m1 <- m1[,!grepl("mean", colnames(m1)), drop = FALSE]
   m2 <- m2[,!grepl("mean", colnames(m2)), drop = FALSE]

   counts.test <- cbind(m1[,grepl("counts", colnames(m1)), drop = FALSE],
                        m2[,grepl("counts", colnames(m2)), drop = FALSE])
   dmrs_new <- counts.test[which.sig,]

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

   colnames(dmrs_new)[1:tn1] <- paste0(lab1, "_", colnames(dmrs_new)[1:tn1])
   colnames(dmrs_new)[(tn1 + 1):(tn1 + tn2)] <- paste0(lab2, "_", colnames(dmrs_new)[(tn1 + 1):(tn1 + tn2)])
   colnames(dmrs_new) <- gsub(".rpkm|.counts", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub(".sorted.bam", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub("\\.1", "", colnames(dmrs_new))
   colnames(dmrs_new) <- gsub("_PDAC1", "", colnames(dmrs_new))
   rm(m1)
   rm(m2)

   # binary classification
   # glmnet 

   library(glmnet)
   
   cvob1 = tryCatch({
     cv.glmnet(x=t(log(dmrs+1)),y=c(rep(1, n1), rep(0, n2)),
        family="binomial", alpha=0.5, type.measure="auc", 
        nfolds = 3, lambda = seq(0,0.05,by = 0.01), standardize=FALSE)
   }, error = function(e) {
       tryCatch({
         cv.glmnet(x=t(log(dmrs+1)),y=c(rep(1, n1), rep(0, n2)),
          family="binomial", alpha=0.5, type.measure="auc", 
          nfolds = 3, lambda = seq(0,0.05,by = 0.01), standardize=FALSE)
       }, error = function(e) {
         NULL
       })
   })

   # best coefficient
   if (length(cvob1)>0){
     new <- predict(cvob1, newx=t(log(dmrs_new+1)), type = "response", 
      s = "lambda.min")
  
     library(ROCR)
     pred <- prediction(new, as.numeric(grepl(lab1, colnames(dmrs_new))))
     auc <- unlist(performance(pred,"auc")@y.values)
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

   if (lab1 == "rcc" || lab2 == "rcc"){
    rccsamps <- which(grepl("rcc", colnames(dmrs_new)))
    x <- match(gsub("rcc_", "", colnames(dmrs_new)[rccsamps]), 
      meta$`Sample number`)
    subtype = c(type[names(type) != "rcc"], "#E69F00", "#56B4E9", "#009E73")
    names(subtype) = c("control", "clear cell", "papillary", "chromophobe")
    st <- colnames(dmrs_new)
    st[rccsamps] <- meta$Histology[x]
    st[-rccsamps] <- "control"

    ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs_new)), lab1, lab2), 
      Subtype = st, Prob = ret_tab$class_prob, Class = ret_tab$class_label),
    col = list(Type = type, 
      Subtype = subtype,
      Prob = probcol,
      Class = classcol))
   }

   ht = Heatmap(log(dmrs_new+1), name = "log(CPM+1)", 
     top_annotation = ha_column, col = ecolors,
     show_row_names = FALSE, show_column_names = colnames,
     column_names_gp = gpar(fontsize = 9),
     column_title = paste0("Top ", top, 
       ifelse(merge, " merged", "")))

   pdf(heatmap.file.test, width = 8)
   draw(ht)
   dev.off()

   if(saveprobs){
     message("saving sample level prob table")
     write.table(ret_tab, quote=FALSE, row.names=FALSE,
      file=file.path(out.dir, paste0("sampleprob_table_", lab1, "_", lab2,
      "_top", top, ".txt")), 
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

    m1 <- MEDIPS.meth(MSet1 = add_samp)[which.sig,]

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
      column_names_gp = gpar(fontsize = 9))

    pdf(heatmap.file.add, width = 8)
    draw(ht)
    dev.off()

  }

  if (!is.null(validate1) || !is.null(validate2)){
  
   if (!is.null(validate1)){
     m1 <- MEDIPS.meth(MSet1 = validate1)
     m1 <- m1[,grepl("rpkm", colnames(m1)), drop = FALSE]
     m1 <- m1[,!grepl("mean", colnames(m1)), drop = FALSE]
   }else{
     m1 <- data.frame()
   }

   if (!is.null(validate2)){
     m2 <- MEDIPS.meth(MSet1 = validate2)
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
   cvob1=cv.glmnet(x=t(dmrs),y=as.numeric(grepl(lab1, colnames(dmrs))),
     family="binomial", alpha=0.5, type.measure="auc", 
     nfolds = 3, lambda = seq(0,0.05,by = 0.01), standardize=FALSE)

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
   ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs_new)), lab1, gsub("_partial", "",lab2)), 
     Prob = ret_tab$class_prob, Class = ret_tab$class_label),
     col = list(Type = type, Prob = probcol, Class = classcol))

   ht = Heatmap(log(dmrs_new+1), name = "log(RPKM+1)", 
     top_annotation = ha_column, col = ecolors,
     show_row_names = FALSE, show_column_names = colnames,
     column_names_gp = gpar(fontsize = 9),
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

rcc.ctrl <- compute.diff(obj1 = medip.rcc, obj2 = medip.control,
           lab1 = "rcc", lab2 = "control",
           training = 0.8, top = ntop,
           saveprobs = TRUE)

urcc.uctrl <- compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
           lab1 = "urineR", lab2 = "urineC",
           training = 0.8, top = ntop,
           saveprobs = TRUE)

if (!file.exists(file.path(outdir, paste0("accuracy_table_top", ntop, ".txt")))){
## replication method
## full set
if (iter == 1){
  compute.diff(obj1 = medip.rcc, obj2 = medip.control,
	         lab1 = "rcc", lab2 = "control",
           out.dir = medipdir, top = ntop)

  compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
	         lab1 = "urineR", lab2 = "urineC",
           out.dir = medipdir, top = ntop)
}

#### training/test version of above (80%):
rcc.ctrl <- compute.diff(obj1 = medip.rcc, obj2 = medip.control,
	         lab1 = "rcc", lab2 = "control",
	         training = 0.8, top = ntop)

urcc.uctrl <- compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
	         lab1 = "urineR", lab2 = "urineC",
	         training = 0.8, top = ntop)

## merged
## full set
if (iter == 1){
  compute.diff(obj1 = medip.rcc, obj2 = medip.control,
           lab1 = "rcc", lab2 = "control", merge = TRUE,,
           out.dir = medipdir, top = ntop)

  compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
           lab1 = "urineR", lab2 = "urineC",
           merge = TRUE,
           out.dir = medipdir, top = ntop)
}

#### training/test version of above (80%):
rcc.ctrl.merge <- compute.diff(obj1 = medip.rcc, obj2 = medip.control,
           lab1 = "rcc", lab2 = "control", merge = TRUE, 
           training = 0.8, top = ntop)


urcc.uctrl.merge <- compute.diff(obj1 = medip.urineR, obj2 = medip.urineC,
           lab1 = "urineR", lab2 = "urineC",
           merge = TRUE, 
           training = 0.8, top = ntop)

# save results table
tab <- rbind(rcc.ctrl, urcc.uctrl, rcc.ctrl.merge, urcc.uctrl.merge)
tab$iteration <- iter
tab$window <- ws
tab$method <- c("original", "original", "merge", "merge") 

write.table(tab, quote=FALSE, row.names=FALSE,
  file=file.path(outdir, paste0("accuracy_table_top", ntop, ".txt")), 
    sep = "\t")
}else{
  message("Summary results file already exists")
}


# blca vs rcc - only for ws = 300, non-merge
if (ws == 300){
  medip.blca <- readRDS(file.path(medipdir, "medip.blca.rds"))
  rcc.blca <- compute.diff(obj1 = medip.rcc, obj2 = medip.blca,
           lab1 = "rcc", lab2 = "blca",
           training = 0.8, top = ntop)
  tab <- rcc.blca
  tab$iteration <- iter
  tab$window <- ws
  tab$method <- "original"

  write.table(tab, quote=FALSE, row.names=FALSE,
  file=file.path(outdir, paste0("accuracy_table_top", ntop, "_blca.txt")), 
    sep = "\t")


  # validation with blca
  if (!file.exists(file.path(medipdir, 
      paste0("validation_accuracy_blca_table_top", ntop, ".txt")))){

    met <- compute.diff(obj1 = medip.rcc, obj2 = medip.control,
             lab1 = "rcc", lab2 = "control",
             top = ntop, out.dir = medipdir,
             validate1 = medip.blca, 
             validate_lab = "vBLCA",
             colnames = TRUE)

    met$iteration <- iter
    met$window <- ws
    met$method <- "original"
    
    write.table(met, quote=FALSE, row.names=FALSE,
    file=file.path(medipdir, 
      paste0("validation_accuracy_blca_table_top", ntop, ".txt")), 
      sep = "\t")
  }
}

