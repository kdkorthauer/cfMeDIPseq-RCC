library(TCGAbiolinks)  
library(SummarizedExperiment)
library(minfi)
library(GEOquery)
library(filesstrings)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(ggplot2)

# doesn't find the supp files for some reason
#getGEOSuppFiles("GSE67393", 
#  baseDir="/scratch/st-kdkortha-1/tmp")

# get them manually
if (!file.exists("/scratch/st-kdkortha-1/tmp/GSE67393/idat/GSM1649832_9376537158_R06C02_Red.idat")){

  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67393/suppl/GSE67393_RAW.tar",
    destfile="/scratch/st-kdkortha-1/tmp/GSE67393_RAW.tar")

  dir.create("/scratch/st-kdkortha-1/tmp/GSE67393/idat")

  untar("/scratch/st-kdkortha-1/tmp/GSE67393_RAW.tar", 
    exdir = "/scratch/st-kdkortha-1/tmp/GSE67393/idat")

  idatFiles <- list.files("/scratch/st-kdkortha-1/tmp/GSE67393/idat", 
    pattern = "idat.gz", full.names = TRUE)

  sapply(idatFiles, gunzip, overwrite = TRUE)


  gset <- getGEO("GSE67393", GSEMatrix =TRUE, getGPL=FALSE)
  saveRDS(gset[[1]], file.path("/scratch/st-kdkortha-1/tmp/GSE67393/GSE67393_gset.rds"))
}else{
  gset <- readRDS(file.path("/scratch/st-kdkortha-1/tmp/GSE67393/GSE67393_gset.rds"))
}


# tcga processed (to get sex metadata)
tcga.dir <- "/scratch/st-kdkortha-1/tmp/TCGA"

if (!file.exists(file.path(tcga.dir, "TCGA.RData"))){
  ### Get TCGA KIRC data
  
  dir.create(tcga.dir, showWarnings = FALSE)

  query_meth.hg19 <- GDCquery(project= "TCGA-KIRC", 
   data.category = "DNA methylation", 
   platform = "Illumina Human Methylation 450",
   sample.type = c("Primary Tumor"),
   legacy = TRUE)

  GDCdownload(query_meth.hg19, 
  method = "api", 
  directory = tcga.dir,
    files.per.chunk = 20)
  dat <- GDCprepare(query_meth.hg19, 
    directory = tcga.dir, save = FALSE)
  saveRDS(dat, file = file.path(tcga.dir, "TCGA.rds"))
}else{
  dat <- readRDS(file=file.path(tcga.dir, "TCGA.rds"))
}



# Same for TCGA data
tcga.dir <- "/scratch/st-kdkortha-1/tmp/TCGA"

if (!file.exists("/scratch/st-kdkortha-1/tmp/TCGA/TCGA-KIRC/legacy/Raw_microarray_data/Raw_intensities/6042324043_R04C02_Grn.idat")){
  
  dir.create(tcga.dir, showWarnings = FALSE)

  query_meth.hg19 <- GDCquery(project= "TCGA-KIRC", 
    data.category = "Raw microarray data", 
    platform = "Illumina Human Methylation 450",
    sample.type = c("Primary Tumor"),
    file.type = ".idat",
    legacy = TRUE)

  GDCdownload(query_meth.hg19, 
    method = "api", 
    directory = tcga.dir,
    files.per.chunk = 20)

  # move everything up one dir
  idatFiles <- list.files("/scratch/st-kdkortha-1/tmp/TCGA/TCGA-KIRC/legacy/Raw_microarray_data/Raw_intensities/", 
    pattern = "idat", full.names = TRUE, recursive = TRUE)

  file.move(idatFiles, 
    file.path(dirname(idatFiles), ".."))

  # mapping from dir to R id
  tcga.dat <- query_meth.hg19$results[[1]]
  saveRDS(tcga.dat, file = file.path(tcga.dir, "querydat_TCGA.rds"))
}else{
  tcga.dat <- readRDS(file.path(tcga.dir, "querydat_TCGA.rds"))
}

# combine metadat for tcga
dup.samps <- names(which(table(tcga.dat$sample.submitter_id) > 2))
tcga.dat <- tcga.dat[!(tcga.dat$sample.submitter_id %in% dup.samps),]
tcga.dat$sample_id <- gsub("_Red.idat|_Grn.idat", "",tcga.dat$file_name)
tcga.dat <- tcga.dat %>% 
  select(sample.submitter_id, sample_id) %>% 
  distinct() %>%
  rename(bcr_patient_barcode=sample.submitter_id)

dat <- dat[,!(colData(dat)$bcr_patient_barcode %in% dup.samps)]
cdat <- colData(dat) %>% as.data.frame() %>%
 dplyr::select(bcr_patient_barcode, gender) %>%
 left_join(tcga.dat, by ="bcr_patient_barcode")

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)


library(RColorBrewer)
addalpha <- function(colors, alpha=1.0) {
    r <- col2rgb(colors, alpha=T)
    # Apply alpha
    r[4,] <- alpha*255
    r <- r/255.0
    return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

if (!file.exists(file.path(tcga.dir, "mset.rds"))){

  rgSet_pbmc <- read.metharray.exp("/scratch/st-kdkortha-1/tmp/GSE67393/idat")
  rgSet_rcc <- read.metharray.exp(base = "/scratch/st-kdkortha-1/tmp/TCGA/TCGA-KIRC/legacy/Raw_microarray_data/Raw_intensities/", recursive=TRUE)

  # check sample order (pbmc)
  sum(!sapply(1:ncol(gset), 
    function(x) grepl(colnames(gset)[x], colnames(rgSet_pbmc)[x])))

  # only keep samples 1-93
  rmv <- c("GSM1646161", "GSM1646162", "GSM1646163", "GSM1646164", "GSM1646165",
           "GSM1646166", "GSM1646167", "GSM1646168", "GSM1646169", "GSM1646170", 
           "GSM1646171", "GSM1646172", "GSM1646173", "GSM1646174", "GSM1646175", 
           "GSM1646176", "GSM1646177", "GSM1646178", "GSM1646179", "GSM1646180", 
           "GSM1646181", "GSM1646182", "GSM1646183", "GSM1646184")
  rgSet_pbmc <- rgSet_pbmc[, !grepl(paste0(rmv, collapse="|"), colnames(rgSet_pbmc))]
  gpd <- as(phenoData(gset), "data.frame") 
  gset <- gset[, !(phenoData(gset)$geo_accession %in% rmv)]

  cdat2 <- data.frame(sample_id = colnames(rgSet_rcc)) %>%
    left_join(cdat, by="sample_id")
  rgSet_rcc <- rgSet_rcc[, !is.na(cdat2$gender)]
  cdat2 <- na.omit(cdat2)

  pdat <- DataFrame(group = c( rep("pbmc", ncol(rgSet_pbmc)), 
                               rep("rcc", ncol(rgSet_rcc)) ),
                    sex = tolower(c(phenoData(gset)$"gender:ch1",
                            cdat2$gender)))

  rgSet <- cbind(rgSet_pbmc, rgSet_rcc)
  

  detP <- detectionP(rgSet)
  head(detP)
  keep <- colMeans(detP) < 0.05
  rgSet <- rgSet[,keep]
  detP <- detP[,keep]
  dim(detP)
  
  grSet <- preprocessNoob(rgSet) 

  pData(grSet)$group <- pdat$group
  pData(grSet)$knownSex <- pdat$sex

  # ensure probes are in the same order in the mSetSq and detP objects
  detP <- detP[match(featureNames(grSet),rownames(detP)),] 

  # remove any probes that have failed in one or more samples
  keep <- rowSums(detP < 0.01) == ncol(grSet) 
  table(keep)

  grSetFilt <- grSet[keep,]

  keep <- !(featureNames(grSetFilt) %in% ann450k$Name[ann450k$chr %in% 
                                                        c("chrX","chrY")])
  table(keep)
  grSetFilt <- grSetFilt[keep,]

  grSetFilt
  saveRDS(grSetFilt, file.path(tcga.dir, "mset.rds"))

  pdf(file.path(tcga.dir, "ssnoob_density.pdf"), width=8, height=4)
  par(mfrow=c(1,2))
  densityPlot(rgSet, sampGroups=pData(grSet)$group, main="Raw", legend=FALSE,
    pal = addalpha(RColorBrewer::brewer.pal(8, "Dark2"), alpha=0.1))
  densityPlot(getBeta(grSetFilt), sampGroups=pData(grSet)$group,
            main="Normalized", legend=FALSE, 
            pal = addalpha(RColorBrewer::brewer.pal(8, "Dark2"), alpha=0.1))
  dev.off()

}else{
  grSetFilt <- readRDS(file.path(tcga.dir, "mset.rds"))
}


# this is the factor of interest
pData(grSetFilt)$cellType <- factor(pData(grSetFilt)$group)
# this is the individual effect that we need to account for
pData(grSetFilt)$sex <- factor(pData(grSetFilt)$knownSex) 

# use the above to create a design matrix
design <- model.matrix(~cellType*sex, data=pData(grSetFilt))
colnames(design) <- c("intercept", "rcc","male", "rcc:male")
 
# fit the linear model 
library(limma)
mVals <- getM(grSetFilt)
bVals <- getBeta(grSetFilt)
fit <- lmFit(mVals, design)
fit <- eBayes(fit)

bVals_rcc <- bVals[,pData(grSetFilt)$group == "rcc"]
bVals_pbmc <- bVals[,pData(grSetFilt)$group == "pbmc"]
delta_beta <- rowMeans(bVals_rcc) - rowMeans(bVals_pbmc)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit))


ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
ann450kSub$deltaBeta <- delta_beta 
DMPs <- topTable(fit, num=Inf, coef="rcc", genelist=ann450kSub)
head(DMPs)

DMP_gr <- makeGRangesFromDataFrame(DMPs %>% mutate(start=pos, end=pos),
  keep.extra.columns=TRUE)



### Intersect with MEDIPS data

outdir <- "/scratch/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_300/revisions"

plasma <- readRDS(file.path(outdir, "pooled", "rcc.control.diff.rds"))
urine <- readRDS(file.path(outdir, "B1","urineR.urineC.diff.rds"))

keepCols <- c("chr", "start", "stop", 
            "MSets1.counts.mean", "MSets1.rpkm.mean",
            "MSets2.counts.mean", "MSets2.rpkm.mean",
            "t",
            "logFC", "P.Value", "adj.P.Val")
plasma <- plasma[,keepCols]
urine <- urine[,keepCols]

colnames(plasma)[3] <- "end"
colnames(urine)[3] <- "end"

plasma <- makeGRangesFromDataFrame(plasma, 
  keep.extra.columns=TRUE)
urine <- makeGRangesFromDataFrame(urine, 
  keep.extra.columns=TRUE)

# RCC vs array

plasma <- plasma[!is.na(plasma$logFC),]
urine <- urine[!is.na(urine$logFC),]
ol_plasma <- findOverlaps(plasma, DMP_gr)

plasma <- plasma[ol_plasma@from,]
dat_plasma <- DMP_gr[ol_plasma@to,]

cf <- ifelse(plasma$adj.P.Val < 0.01 & plasma$logFC > 1, "up",
             ifelse(plasma$adj.P.Val < 0.01 & plasma$logFC < -1, "down", 
              "no change"))

array <- ifelse(dat_plasma$adj.P.Val < 0.01 & dat_plasma$deltaBeta > 0.25, "up",
             ifelse(dat_plasma$adj.P.Val < 0.01 & dat_plasma$deltaBeta < -0.25, "down", 
              "no change"))

table(cf, array)

df <- data.frame(plasma = plasma$logFC,
  array = dat_plasma$deltaBeta)
ggplot(data=df, aes(x=array, y = plasma)) +
 geom_hex(bins=100) + 
 scale_fill_gradientn(colours=c("blue","orange"),
  name = "Frequency",na.value=NA,trans = "log10")+
 geom_hline(yintercept = 0) + 
 geom_vline(xintercept = 0) +
 geom_smooth(method="lm", se=FALSE, color="red") +
 theme_bw() +
 xlim(-1,1) + ylim(-3,4.5)+
 xlab("delta-Beta TCGA-KIRC vs PBMCs (450K array)") +
 ylab("logFC case vs control (cfMeDIPseq)") +
 annotate("text", -0.75, 4, 
  label = paste0("Correlation=",
    signif(cor.test(plasma$logFC, dat_plasma$deltaBeta, method="spearman")$estimate, 2)))

# remove sites that are methylated in healthy plasma
# file from Ankur (email 3/20/2020)
unmeth_plasma <- load(file.path(tcga.dir, "UnmethylatedInPlasma.RData")) 
keep <- which(dat_plasma$Name %in% Probes2.bloodLow$Probe)
df <- df[keep,]

ggplot(data=df, aes(x=array, y = plasma)) +
 geom_hex(bins=100) + 
 scale_fill_gradientn(colours=c("blue","orange"),
  name = "Frequency",na.value=NA,trans = "log10")+
 geom_hline(yintercept = 0) + 
 geom_vline(xintercept = 0) +
 geom_smooth(method="lm", se=FALSE, color="red") +
 theme_bw() +
 xlim(-1,1) + ylim(-3,4.5)+
 xlab("delta-Beta TCGA-KIRC vs PBMCs (450K array)") +
 ylab("logFC RCC vs control (cfMeDIPseq)") +
 annotate("text", -0.75, 4.5, 
  label = paste0("Correlation = ",
    signif(cor.test(plasma$logFC[keep], dat_plasma$deltaBeta[keep], 
      method="spearman")$estimate, 2)))

ggsave(file.path(tcga.dir,"tcga_pbmc_corr_plasma.pdf"), width=6, height=5)


# urine vs array
ol_urine <- findOverlaps(urine, DMP_gr)

urine <- urine[ol_urine@from,]
dat_urine <- DMP_gr[ol_urine@to,]

df <- data.frame(array = urine$logFC,
  urine = dat_urine$deltaBeta)
keep <- which(dat_urine$Name %in% Probes2.bloodLow$Probe)
df <- df[keep,]

ggplot(data=df, aes(x=array, y = urine)) +
 geom_hex(bins=100) + 
 scale_fill_gradientn(colours=c("blue","orange"),
  name = "Frequency",na.value=NA,trans = "log10")+
 geom_hline(yintercept = 0) + 
 geom_vline(xintercept = 0) +
 geom_smooth(method="lm", se=FALSE, color="red") +
 theme_bw() +
 xlab("delta-Beta TCGA-KIRC vs PBMCs (450K array)") +
 ylab("logFC RCC vs control (urine)") +
 annotate("text", -0.75, 1, 
  label = paste0("Correlation = ",
    signif(cor.test(urine$logFC[keep], dat_urine$deltaBeta[keep], 
      method="spearman")$estimate, 2)))

ggsave(file.path(tcga.dir,"tcga_pbmc_corr_urine.pdf"), width=6, height=5)




# plasma vs urine

plasma <- readRDS(file.path(outdir, "pooled", "rcc.control.diff.rds"))
urine <- readRDS(file.path(outdir, "B1","urineR.urineC.diff.rds"))
plasma <- plasma[,keepCols]
urine <- urine[,keepCols]
colnames(plasma)[3] <- "end"
colnames(urine)[3] <- "end"
plasma <- makeGRangesFromDataFrame(plasma, 
  keep.extra.columns=TRUE)
urine <- makeGRangesFromDataFrame(urine, 
  keep.extra.columns=TRUE)
plasma <- plasma[!is.na(plasma$logFC),]
urine <- urine[!is.na(urine$logFC),] 

ol_plasma <- findOverlaps(plasma, urine)

plasma <- plasma[ol_plasma@from,]
urine <- urine[ol_plasma@to,]

df <- data.frame(plasma = plasma$logFC,
  urine = urine$logFC)
ggplot(data=df, aes(x=urine, y = plasma)) +
 geom_hex(bins=100) + 
 scale_fill_gradientn(colours=c("blue","orange"),
  name = "Frequency",na.value=NA,trans = "log10")+
 geom_hline(yintercept = 0) + 
 geom_vline(xintercept = 0) +
 geom_smooth(method="lm", se=FALSE, color="red") +
 geom_smooth(se=FALSE, color="black", linetype="dashed") +
 theme_bw() +
 xlab("logFC case vs control (urine)") +
 ylab("logFC case vs control (plasma)") +
 annotate("text", -1.25, 5, 
  label = paste0("Correlation=",
    signif(cor.test(plasma$logFC, urine$logFC, method="spearman")$estimate, 2)))
ggsave(file.path(tcga.dir,"plasma_urine_corr.pdf"), width=6, height=5)
