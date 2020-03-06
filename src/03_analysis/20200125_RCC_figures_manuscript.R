# read in iter accuracy tables created by 20190712_DMRs_RCC.R
# plot results

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(cowplot)
library(readxl)
library(Rsamtools)
library(ComplexHeatmap)
library(DESeq2)
library(RColorBrewer)
library(edgeR)
library(gtable)
library(gridGraphics)
library(grid)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(scatterplot3d)
library(MEDIPS)
library(grDevices)


theme_set(theme_bw())

top <- 300
ws <- 300
merge <- FALSE
res <- NULL

# dir where binned medips objects of all samples are saved
outdir <- paste0("/arc/project/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws)
savedir <- paste0("/scratch/st-kdkortha-1/cfMeDIPseq/out/MEDIPS_", ws, "/revisions/pooled")

medip.rcc <- readRDS(file.path(outdir, "medip.rcc.rds"))
medip.rcc <- medip.rcc[!sapply(medip.rcc, function(x) x@sample_name) %in% c("S007.sorted.bam", "S011.sorted.bam", "S020.sorted.bam", "S015.sorted.bam", "S050.sorted.bam")]
medip.control <- readRDS(file.path(outdir, "medip.control.rds"))
medip.control <- medip.control[!sapply(medip.control, function(x) x@sample_name) %in% c("S040.sorted.bam")]


medip.urineR <- readRDS(file.path(outdir, "medip.urineR.rds"))
medip.urineC <- readRDS(file.path(outdir, "medip.urineC.rds"))
medip.urineC <- medip.urineC[!sapply(medip.urineC, function(x) x@sample_name) %in% c("S56.sorted.bam", "S4.sorted.bam", "S6.sorted.bam", "S37.sorted.bam")]
medip.urineR <- medip.urineR[!sapply(medip.urineR, function(x) x@sample_name) %in%
    c("S11.sorted.bam", "S12.sorted.bam", "S30.sorted.bam", "S33.sorted.bam", 
    "S34.sorted.bam", "S48.sorted.bam", "S59.sorted.bam", "S60.sorted.bam",
    "S62.sorted.bam", "S64.sorted.bam", "S42.sorted.bam", "S25.sorted.bam",
    "S14.sorted.bam")]

medip.jan2020 <- readRDS(file.path(outdir, "medip.jan2020.rds")) # already filtered


# create pooled set

# joint metadata for RCC samps
meta2 <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/20200108/20200108_Sample List.xlsx")
m2 <- data.frame(ID=gsub(".sorted.bam", "", sapply(medip.jan2020, function(x) x@sample_name))) %>%
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

table(master$Status, master$Institution, master$Source)

# remove three RCCMet samples since don't have histology
medip.rcc_M <- readRDS(file.path(outdir, "medip.rcc_M.rds")) # metastatic
metids <- gsub(".sorted.bam", "", sapply(medip.rcc_M, function(x) x@sample_name))
x <- match( metids, master$"Sample number" )
excl <- master$"Sample number"[x][grepl("Exclude", master$Inclusion[x])]
if(length(excl) > 0)
  medip.rcc_M <- medip.rcc_M[-which(metids %in% excl)]


# joint medips objs
medip.rcc <- c(medip.rcc, medip.rcc_M)
medip.urineR <- c(medip.urineR, medip.jan2020[m2$Source=="Urine" & m2$Status == "RCC"])
medip.urineC <- c(medip.urineC, medip.jan2020[m2$Source=="Urine" & m2$Status == "Control"])


################# volcanoes

# plasma
plotVolcano <- function(diff.file, sig=0.1){
  diff <- readRDS(file=diff.file)
  dmrs <- diff[, grepl("adj.p|FC|P.Value", colnames(diff))]
  dmrs <- dmrs[!is.na(dmrs$logFC),]

  #plot
  annotations <- data.frame(
        xpos = c(-Inf,Inf),
        ypos =  c(Inf,Inf),
        annotateText = c("Loss","Gain"),
        hjustvar = c(-0.8,1.5) ,
        vjustvar = c(2,2)) #<- adjust

  # reverse direction to reflect gain = gain in RCC 
  dmrs$logFC <- -1*dmrs$logFC

  message(sum(dmrs$limma.adj.p.value < sig & dmrs$logFC < 0), " Lost")
  message(sum(dmrs$limma.adj.p.value < sig & dmrs$logFC > 0), " Gained")


  which.up <- which(dmrs$logFC > 0 & !is.na(dmrs$logFC))
  which.down <- which(dmrs$logFC < 0 & !is.na(dmrs$logFC))

  which.sig.up <- which(rank(dmrs$P.Value[which.up], 
     ties.method = "random") <= as.numeric(top)/2)

  which.sig.down <- which(rank(dmrs$P.Value[which.down], 
      ties.method = "random") <= as.numeric(top)/2)

  message("there are ", 
    sum(diff$limma.adj.p.value < 0.05 & !is.na(diff$P.Value)),
    " dmrs at FDR 0.05.")
  message("there are ", 
    sum(diff$limma.adj.p.value < 0.01 & !is.na(diff$P.Value)),
    " dmrs at FDR 0.01.")

  message("qval threshold for top 300 is: ",
    max(dmrs$limma.adj.p.value[which.up][which.sig.up], na.rm=TRUE), " (up)",
    max(dmrs$limma.adj.p.value[which.down][which.sig.down], na.rm=TRUE), " (down)")

  allsig <- c(which.up[which.sig.up], which.down[which.sig.down])

  plt <- ggplot(dmrs, aes(x=logFC, y=-log10(limma.adj.p.value))) +
      stat_binhex(data=dmrs[-allsig,], bins=150) +
      scale_fill_viridis_c(direction=-1) +
      geom_point(data=dmrs[allsig,],
        color="red", size=0.6, alpha=0.4) +
      xlab("DNA Methylation (log2 Fold Change)") +
      ylab("Significance (-log10 q-value)") + 
      geom_hline(yintercept=-log10(sig), 
        colour="darkgrey", linetype="dashed") +
      geom_vline(xintercept=0) +
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,
        vjust=vjustvar,label=annotateText))+
      geom_text(data=annotations,aes(x=1.6,y=-log10(sig),
        label=paste0("FDR<", sig)))

  return(plt)
}


volcano_plasma <- plotVolcano(diff.file =file.path(savedir, "rcc.control.diff.rds"),
  sig=0.01)
volcano_plasma
ggsave(file.path(savedir, "volcano_plasma.pdf"), width=5, height=4)

volcano_urine <- plotVolcano(diff.file =file.path(savedir, "urineR.urineC.diff.rds"), 
  sig=0.05)
volcano_urine
ggsave(file.path(savedir, "volcano_urine.pdf"), width=5, height=4)



### PCA plots - normalize on all regions meeting low threshold for counts (mean 0.25)

depths <- function(mdobjlist, CS, type){
  depth <- data.frame(sapply(mdobjlist, function(x){
      x@genome_count[CS@genome_CF >= 0]
  }))
  colnames(depth) <- sapply(mdobjlist, function(x){ 
    gsub(".bam|.sorted.bam", "", x@sample_name)})
  colnames(depth) <- paste0(type, "_", colnames(depth))
  return(depth)
}

# avoid code repetition - package PCA plots into a function
# returns tidy data frame with pcs 1-4 and pct zero for each sample
makePCAplots <- function(diff.file, obj1, obj2, label, ntop=NULL){

  if(label=="plasma"){
    l1 <- "rcc"
    l2 <- "ctrl"
  }else if (label=="urine"){
    l1 <- "urineR"
    l2 <- "urineC"
  }

  df <- cbind(depths(obj1, CS, l1),
              depths(obj2, CS, l2))

  pctzero_overall <- colMeans2(df==0)
  total_overall <- colSums2(as.matrix(df))

  diff <- readRDS(file=diff.file)
  
  if(ncol(diff) != 3+2*(length(obj1)+length(obj2))+11)
   message("WARNING; diff has ", ncol(diff), " columns. ",
     "Expecting ", 3+2*(length(obj1)+length(obj2))+11, ".")

  if (!is.null(ntop)){
    which.up <- which(diff$logFC > 0)
    which.down <- which(diff$logFC < 0)
    which.sig.up <- which(rank(diff$P.Value[which.up], 
       ties.method = "random") <= as.numeric(top)/2)

    which.sig.down <- which(rank(diff$P.Value[which.down], 
       ties.method = "random") <= as.numeric(top)/2)

    which.sig <- c(which.up[which.sig.up], which.down[which.sig.down])

    rm(diff)
  }else{  
    keep <- which(!is.na(diff$P.Value))
    rm(diff)
  }

  # normalize on all genes
  dge <- DGEList(counts=df[which(rowSums(df,na.rm=TRUE)>=0.25*ncol(df)),])
  dge <- calcNormFactors(dge, refColumn = 1)

  # subset to dmrs
  if (!is.null(ntop)){
    df <- df[which.sig,] 
    pctzero_top <- colMeans2(df==0)
    total_top <- colSums2(as.matrix(df))
  }else{
    # still need to remove those where most counts are zero
    pctzero_top <- total_top <- NA
    df <- df[keep,] 
    df <- df[which(rowSums(df,na.rm=TRUE)>=0.25*ncol(df)),]
  }

  grp <- gsub("_.*", "", colnames(df))
  ids = colnames(df)

  if (label=="plasma"){
    x <- match(ids, paste0(ifelse(master$Status == "RCC", "rcc_", "ctrl_"), master$`Sample number`))
  }else if (label=="urine"){
    x <- match(ids, paste0(ifelse(master$Status == "RCC", "urineR_", "urineC_"), master$`Sample number`))
  }

  subtype <- as.character(master$Histology[x])
  subtype <- ifelse(is.na(subtype), master$Status[x], subtype)
  batch <- as.character(master$Batch[x])

  df <- sweep(df, MARGIN=2, 
    (dge$samples$norm.factors*dge$samples$lib.size)/1e6, `/`)

  pcs <- Morpho::prcompfast(t(log(df+1)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids,
           Batch = as.factor(batch)) %>%
    mutate(Type = ifelse(type == "rcc" | type == "urineR", "RCC", "Control"))%>%
    mutate(pctZero_top = pctzero_top,
           pctZero_overall = pctzero_overall,
           total_top = total_top,
           total_overall = total_overall)
    
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)


  if (!is.null(ntop)){
    label <- paste0(label, "_top", ntop)
  }
  tidydf %>% ggplot(aes(x=PC1, y=PC2, colour = Type, shape = Batch)) +
     geom_point(size = 2) +
     scale_color_manual(values = colors) 
  ggsave(file.path(savedir, paste0("PCA_1_2_", label, ".pdf")), width = 4.5, height = 3.5)

  tidydf %>% ggplot(aes(x=PC2, y=PC3, colour = Type, shape = Batch)) +
     geom_point(size = 2)+
     scale_color_manual(values = colors)
  ggsave(file.path(savedir, paste0("PCA_2_3_", label, ".pdf")), width = 4.5, height = 3.5)

  tidydf %>% ggplot(aes(x=PC1, y=PC3, colour = Type, shape = Batch)) +
     geom_point(size = 2)+
     scale_color_manual(values = colors)
  ggsave(file.path(savedir, paste0("PCA_1_3_", label, ".pdf")), width = 4.5, height = 3.5)

  pdf(file.path(savedir, paste0("PCA_1_2_3_", label, ".pdf")), width = 4.5, height = 4.5)
   colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
   colors <- colors[as.numeric(as.factor(tidydf$Type))]
   shape <- c(20,17)[as.numeric(batch)]
   scatterplot3d(tidydf[,1:3], pch = shape, 
     xlab="PC1", ylab="PC2", zlab="PC3", cex.symbols = 2,
     color=colors)
   legend(0,-5.3, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = c(20,17), 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
   dev.off()

  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)

  if (!is.null(ntop)){
    p1 <- ggplot() +
      geom_point(data = tidydf, aes(x=pctZero_top, y=PC1, colour = Type, shape=Batch), size = 2) +
      scale_color_manual(values = colors) + 
      theme(legend.position="none")

    p2 <- ggplot() +
      geom_point(data = tidydf, aes(x=pctZero_top, y=PC2, colour = Type, shape=Batch), size = 2) +
      scale_color_manual(values = colors)  + 
      theme(legend.position="none")
 
    p3 <- ggplot() +
      geom_point(data = tidydf, aes(x=pctZero_top, y=PC3, colour = Type, shape=Batch), size = 2) +
      scale_color_manual(values = colors) 

    leg <- get_legend(p3)
    plot_grid(p1,p2,p3 +  theme(legend.position="none"), leg, nrow=2)
    ggsave(file.path(savedir, paste0("PCA_pctZeroTop", ntop, "_", label, ".pdf")), width = 5, height = 5)
  

     p1 <- ggplot() +
      geom_point(data = tidydf, aes(x=total_top, y=PC1, colour = Type, shape=Batch), size = 2) +
      scale_color_manual(values = colors) + 
      theme(legend.position="none")

    p2 <- ggplot() +
      geom_point(data = tidydf, aes(x=total_top, y=PC2, colour = Type, shape=Batch), size = 2) +
      scale_color_manual(values = colors)  + 
      theme(legend.position="none")
 
    p3 <- ggplot() +
      geom_point(data = tidydf, aes(x=total_top, y=PC3, colour = Type, shape=Batch), size = 2) +
      scale_color_manual(values = colors) 

    leg <- get_legend(p3)
    plot_grid(p1,p2,p3 +  theme(legend.position="none"), leg, nrow=2)
    ggsave(file.path(savedir, paste0("PCA_totalTop", ntop, "_", label, ".pdf")), width = 5, height = 5)
   
    ggplot() +
      geom_point(data = tidydf, aes(x=total_top, y=pctZero_top, colour = Type, shape=Batch), size = 2) +
      scale_color_manual(values = colors) 
    ggsave(file.path(savedir, paste0("PCA_totalTop", ntop, "_vs_pctZero", ntop, "_", label, ".pdf")), width = 5, height = 5)

  }

  p1 <- ggplot() +
    geom_point(data = tidydf, aes(x=pctZero_overall, y=PC1, colour = Type, shape=Batch), size = 2) +
    scale_color_manual(values = colors) + 
    theme(legend.position="none")

  p2 <- ggplot() +
    geom_point(data = tidydf, aes(x=pctZero_overall, y=PC2, colour = Type, shape=Batch), size = 2) +
    scale_color_manual(values = colors)  + 
    theme(legend.position="none")
 
  p3 <- ggplot() +
    geom_point(data = tidydf, aes(x=pctZero_overall, y=PC3, colour = Type, shape=Batch), size = 2) +
    scale_color_manual(values = colors) 

  leg <- get_legend(p3)
  plot_grid(p1,p2,p3 +  theme(legend.position="none"), leg, nrow=2)
  ggsave(file.path(savedir, paste0("PCA_pctZeroOverall_", label, ".pdf")), width = 5, height = 5)


  p1 <- ggplot() +
    geom_point(data = tidydf, aes(x=total_overall, y=PC1, colour = Type, shape=Batch), size = 2) +
    scale_color_manual(values = colors) + 
    theme(legend.position="none")

  p2 <- ggplot() +
    geom_point(data = tidydf, aes(x=total_overall, y=PC2, colour = Type, shape=Batch), size = 2) +
    scale_color_manual(values = colors)  + 
    theme(legend.position="none")
 
  p3 <- ggplot() +
    geom_point(data = tidydf, aes(x=total_overall, y=PC3, colour = Type, shape=Batch), size = 2) +
    scale_color_manual(values = colors) 

  leg <- get_legend(p3)
  plot_grid(p1,p2,p3 +  theme(legend.position="none"), leg, nrow=2)
  ggsave(file.path(savedir, paste0("PCA_totalOverall_", label, ".pdf")), width = 5, height = 5)

  ggplot() +
      geom_point(data = tidydf, aes(x=total_overall, y=pctZero_overall, colour = Type, shape=Batch), 
        size = 2) +
      scale_color_manual(values = colors) 
    ggsave(file.path(savedir, paste0("PCA_totalTop_vs_pctZero_", label, ".pdf")), width = 5, height = 5)

  write.table(data.frame(PC=1:10, Proportion=(pcs$sdev/sum(pcs$sdev))[1:10]), 
      quote=FALSE, row.names=FALSE,
      file=file.path(savedir, paste0("PC_proportionVariation_", label, ".txt")), 
      sep = "\t")

  return(tidydf)
}

CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip.rcc[[1]])

# plasma

# top 300 DMRs
tidydf_plasma <- makePCAplots(diff.file=file.path(savedir, "rcc.control.diff.rds"), 
  obj1=medip.rcc, obj2=medip.control, label="plasma", ntop=top)

# urine

# top 300 DMRs
tidydf_urine <- makePCAplots(diff.file=file.path(savedir, "urineR.urineC.diff.rds"), 
  obj1=medip.urineR, obj2=medip.urineC, label="urine", ntop=top)

################# AUC summary 

# RCCmet
# file list

files <- list.files(file.path(savedir, ".."), pattern = "sampleprob_table", 
  recursive = TRUE, full.names = TRUE)

tmp <- files %>%
  purrr::map(read_tsv)
tmp <- lapply(seq_along(files), 
  function(x) {
    mutate(tmp[[x]], 
      filename=files[x],
      idx=paste0(true_label,1:nrow(tmp[[x]])))
  }) 
tmp <- tmp %>%
  do.call("rbind", .) 
idx <- stringr::str_locate(tmp$filename[1], "iter_")[2]
tmp$iteration = substr(tmp$filename, idx+1, idx+3)
tmp <- mutate(tmp, type = ifelse(true_label %in% c("rcc", "control", "blca"), "Plasma", "Urine"))
tmp <- mutate(tmp, grp = ifelse(grepl("blca", filename), "ubc", "rcc"))

# auc table
auc_summary <- tmp %>%
  group_by(iteration,type,grp) %>%
  summarize(auc=unique(auc)) %>%
  group_by(type, grp) %>%
  summarize(meanAUC = mean(auc, na.rm = TRUE),
            sdAUC = sd(auc, na.rm = TRUE),
            n = n(),
            lowerAUC = max(0, meanAUC - qnorm(0.975)*sdAUC/sqrt(n)),
            upperAUC = min(meanAUC + qnorm(0.975)*sdAUC/sqrt(n), 1))
auc_summary
write.table(auc_summary, quote=FALSE, row.names=FALSE,
  file=file.path(savedir, "..", "auc_summary_100iter.txt"))

# plot of AUC

tmp %>% 
  group_by(iteration,type,grp) %>%
  summarize(auc = mean(auc, na.rm = TRUE)) %>%
  ggplot(aes(x = interaction(type,grp), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=interaction(type,grp)), height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") +
  theme(legend.position = "none")
ggsave(file=file.path(savedir, "..", "boxplot_auc_summary_100iter.pdf"),
  width=4, height=3)


tmp %>% filter(type=="Plasma" & grp=="rcc") %>%
  group_by(iteration, type) %>%
  summarize(auc = mean(auc, na.rm = TRUE)) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") +
  theme(legend.position = "none")

ggsave(file=file.path(savedir, "..", "boxplot_auc_summary_100iter_plasma.pdf"),
  width=4, height=3)



tmp %>% filter(type=="Urine") %>%
  group_by(iteration, type) %>%
  summarize(auc = mean(auc, na.rm = TRUE)) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") +
  theme(legend.position = "none")

ggsave(file=file.path(savedir, "..", "boxplot_auc_summary_100iter_urine.pdf"),
  width=4, height=3)


tmp %>% filter(grp=="ubc") %>%
  group_by(iteration, type) %>%
  summarize(auc = mean(auc, na.rm = TRUE)) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") +
  theme(legend.position = "none")

ggsave(file=file.path(savedir, "..", "boxplot_auc_summary_100iter_urine.pdf"),
  width=4, height=3)




# sample-level score summary
# each sample


# Reorder
#######################################
# Functions for sorting factor levels #
# (janhove.github.io)                 #
#######################################
# Sort factor levels by the factor level mean of another covariate
sortLvlsByVar.fnc <- function(oldFactor, sortingVariable, ascending = TRUE) {
  
  require("dplyr")
  require("magrittr")
  
  # Combine into data frame
  df <- data.frame(oldFactor, sortingVariable)
  
  ###
  ### If you want to sort the levels by, say, the median, sd etc. instead of the mean,
  ### just change 'mean(sortingVariable)' below to, say, 'median(sortingVariable)'.
  ###
  
  # Compute average of sortingVariable and arrange (ascending)
  if (ascending == TRUE) {
    df_av <- df %>% group_by(oldFactor) %>% 
      summarise(meanSortingVariable = median(sortingVariable, na.rm=TRUE)) %>% 
      arrange(meanSortingVariable)
  }
  
  # Compute average of sortingVariable and arrange (descending)
  if (ascending == FALSE) {
    df_av <- df %>% group_by(oldFactor) %>%
      summarise(meanSortingVariable = median(sortingVariable, na.rm=TRUE)) %>% 
      arrange(desc(meanSortingVariable))
  }
  
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = df_av$oldFactor)
  return(newFactor)
}
####  end function to reorder


risk_summary <- tmp %>% group_by(sample_name) %>%
  summarize(mean_RCC_risk_score = mean(class_prob),
            sd_RCC_risk_score = sd(class_prob, na.rm = TRUE),
            n = n(),
            lower_RCC_risk_score = mean_RCC_risk_score - qnorm(0.975)*sd_RCC_risk_score/sqrt(n),
            upper_RCC_risk_score = mean_RCC_risk_score + qnorm(0.975)*sd_RCC_risk_score/sqrt(n),
            true_label = unique(true_label))
write.table(risk_summary, quote=FALSE, row.names=FALSE,
  file=file.path(savedir, "RCC_risk_score_summary_100iter_splitControls.txt"))


cols <- c("RCC" = "#56B4E9", "Control" = "#E69F00")
cols_ubc <- c("RCC" = "#56B4E9", "UBC" = "#E69F00")
sample_probs_all <- tmp %>% 
  mutate(true_label = ifelse(true_label %in% c("urineR", "rcc"), "RCC",
                      ifelse(true_label == "blca", "UBC", "Control"))) %>%
  mutate(sample_name = gsub("control_|rcc_|blca_", "", sample_name)) 


# rcc plasma

sample_probs <- filter(sample_probs_all, grp=="rcc" & type=="Plasma") %>%
  mutate(id = sample_name,
    sample_name = paste0(true_label, "_", sample_name))

sample_probs$sample_name <- as.factor(sample_probs$sample_name)
sample_probs$sample_name <- sortLvlsByVar.fnc(sample_probs$sample_name,
  sample_probs$class_prob)

p_rcc <- sample_probs %>%
  ggplot(aes(x = sample_name, y = class_prob, color = true_label)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  labs(color = "True class") + 
  scale_color_manual(values = cols) + 
  ylab("RCC Risk Score") + xlab("Sample") +
  ylim(0,1)
ggsave(p_rcc,
  file =file.path(savedir, "boxplot_risk_score_summary_100iter_rcc_plasma.pdf"),
  width = 6, height = 3)

rcc <- sample_probs

# rcc urine

sample_probs <- filter(sample_probs_all, grp=="rcc" & type=="Urine") %>%
  mutate(id = sample_name,
    sample_name = paste0(true_label, "_", sample_name))

sample_probs$sample_name <- as.factor(sample_probs$sample_name)
sample_probs$sample_name <- sortLvlsByVar.fnc(sample_probs$sample_name,
  sample_probs$class_prob)

p_urine <- sample_probs %>%
  ggplot(aes(x = sample_name, y = class_prob, color = true_label)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  labs(color = "True class") + 
  scale_color_manual(values = cols) + 
  ylab("RCC Risk Score") + xlab("Sample") +
  ylim(0,1)
ggsave(p_urine, 
  file = file.path(savedir, "boxplot_risk_score_summary_100iter_rcc_urine.pdf"),
  width = 6, height = 3)

urine <- sample_probs

# ubc plasma

sample_probs <- filter(sample_probs_all, grp=="ubc" & type=="Plasma") %>%
  mutate(id = sample_name,
    sample_name = paste0(true_label, "_", sample_name))

sample_probs$sample_name <- as.factor(sample_probs$sample_name)
sample_probs$sample_name <- sortLvlsByVar.fnc(sample_probs$sample_name,
  sample_probs$class_prob)

p_ubc <- sample_probs %>%
  ggplot(aes(x = sample_name, y = class_prob, color = true_label)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  labs(color = "True class") + 
  scale_color_manual(values = cols_ubc) + 
  ylab("RCC Risk Score") + xlab("Sample") +
  ylim(0,1)
ggsave(p_ubc, 
  file=file.path(savedir, "boxplot_risk_score_summary_100iter_ubc_plasma.pdf"),
  width = 6, height = 3)


ubc <- sample_probs

# Add stage/grade/histo annotation bar
library(viridis)
library(wesanderson)

add_annotations <- function(main, dat){
  dat <- dat %>% 
  mutate(Stage = as.factor(ifelse(Stage=="Control", NA, Stage)),
    Grade = as.factor(ifelse(Grade=="Control", NA, Stage)),
    Grade = as.factor(ifelse(is.na(Grade) & true_label == "RCC", "Unknown", Stage)))

  anno_stage <- ggplot(dat %>% 
    distinct(sample_name, Stage)) +
  geom_bar(mapping = aes(x = sample_name, y = 1,
   fill = Stage), 
  stat = "identity", 
  width = 1)+
  scale_fill_viridis_d(direction = -1, 
    guide = guide_legend(
    direction = "horizontal",
    title.position = "top",
    label.position = "bottom",
    label.theme = element_text(angle = 90, hjust = 1, vjust = 0.5, size=6),
    title.theme = element_text(size=8, hjust = 0.5))) +
  theme_void()

  anno_grade <- ggplot(dat %>% 
    distinct(sample_name, Grade)) +
  geom_bar(mapping = aes(x = sample_name, y = 1,
   fill = Grade), 
  stat = "identity", 
  width = 1)+
  scale_fill_manual(values = c(rev(viridis::viridis(5, option = "magma"))[-1], "grey"),
    guide = guide_legend(
    direction = "horizontal",
    nrow = 1,
    title.position = "top",
    label.position = "bottom",
    position= "bottom",
    label.theme = element_text(angle = 90, hjust = 1, vjust = 0.5, size=6),
    title.theme = element_text(size=8, hjust = 0.5))) +
  theme_void()

  rmv <- c(1,2,5)
  if (length(unique(dat$Histology)) == 2)
    rmv <- c(1,2,4,5)

  anno_hist <- ggplot(dat %>% 
    distinct(sample_name, Histology)) +
  geom_bar(mapping = aes(x = sample_name, y = 1,
   fill = Histology), 
  stat = "identity", 
  width = 1)+
  scale_fill_manual(values = c(wes_palette("Rushmore1", n = 5)[-c(1,2,5)], "#E69F00"),
    guide = guide_legend(
    direction = "horizontal",
    title.position = "top",
    label.position = "bottom",
    label.theme = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    title.theme = element_text(size=8, hjust = 0.5))) +
  theme_void()

  legend_anno <- plot_grid(
    get_legend(anno_stage), 
    get_legend(anno_grade), 
    get_legend(anno_hist), 
    get_legend(main), 
    ncol = 1)
  main <- main + theme(legend.position = "none")
  anno_stage <- anno_stage + theme(legend.position = "none")
  anno_grade <- anno_grade + theme(legend.position = "none")
  anno_hist <- anno_hist + theme(legend.position = "none")
  plot <- plot_grid(anno_stage, anno_grade, anno_hist,
    main, align = "v", ncol = 1, axis = "tb", 
    rel_heights = c(0.5, 0.5, 0.5, 15))
  plot_grid(plot, legend_anno, nrow = 1, rel_widths = c(10, 3))

}
 
# add stage/grade/hist to rcc, urine, ubc tibbles
master <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/20200108/20.01.31 - Final Sample List for NM Revisions.xlsx", 
  sheet = 3)
master <- master %>%
  na_if("N/A") %>%
  mutate(Histology = tolower(Histology)) %>%
  mutate(`Sample number` = ifelse(Batch == "Met", tolower(ID), ID)) %>%
  mutate(Histology = ifelse(Histology %in% c("collecting duct", "chrcc", "xptranslocation"), 
    "other", Histology))

# start here
x <- match(paste0(master$Status, "_", master$ID), rcc$sample_name)
master$Status[x]

#add_annotations(plot + 
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank()), plasma)
#ggsave(file.path(outdir, "Boxplot_testset_probs_plasma_annotated.pdf"),
#  width = 7, height = 3.5)



########################################################################
########################################################################
########################################################################
########################################################################



############ Supplementary


## unique vs duplicate reads for plasma & urine

outdir <- "/arc/project/st-kdkortha-1/cfMeDIPseq/out/multiqc"
dat.dir <- "/arc/project/st-kdkortha-1/cfMeDIPseq/out/sortedbam_dup"
setwd(outdir)

groups <- c("RCC", "CONTROL", "URINE_RCC", "URINE_CONTROL", "JAN2020")

tools <- c("fastqc_trimmed", "botwie")

dirs <- list.dirs()

tab_all <- NULL

meta <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/RCC/42 RCC fastq files.xlsx")
meta2 <- read_excel("/arc/project/st-kdkortha-1/cfMeDIPseq/data/20200108/20200108_Sample List.xlsx")


#####
for (group in groups){ 

# fastqc
tool <- tools[1]

thisdir <- dirs[grepl(tool, dirs) & grepl(group, dirs)]

if (length(thisdir)>1){
  thisdir <- thisdir[!grepl("URINE", thisdir)]
  thisdir <- thisdir[!grepl("RCCmet", thisdir)]
}

tab <- read_tsv(file.path(thisdir, "multiqc_general_stats.txt")) %>%
  mutate(id = gsub("_1_val_1|_2_val_2", "", Sample)) %>%
  mutate(mate = ifelse(grepl("_1_val_1", Sample), "1", "2")) %>%
  mutate(Duplicate = `FastQC_mqc-generalstats-fastqc-total_sequences` * 
                          `FastQC_mqc-generalstats-fastqc-percent_duplicates`/100,
         Unique = `FastQC_mqc-generalstats-fastqc-total_sequences` - 
                       Duplicate) %>%
  gather(type, count, 9:10) %>%
  mutate(group=group) %>%
  mutate(group = ifelse(group == "RCC", "PLASMA_RCC", 
    ifelse(group == "CONTROL", "PLASMA_CONTROL", group))) %>%
  mutate(id = gsub("HY2FCBBXX_L005|", "", id)) 

  if (group == "JAN2020"){
    cntk <- colnames(tab)
    tab <- tab %>%
       mutate(ID = gsub("_CKD.*", "", Sample)) %>%
       left_join(meta2, by="ID") %>%
       mutate(group = paste0(group, "_", Source, "_", Status)) %>%
       select(cntk)
  }

  if (group == "RCC"){
    tab <- tab[tab$id %in% meta$`Fastq name`,]
    # remove samples not meeting inclusion criterion
    tab <- filter(tab, ! tab$id %in% c("S007", "S015", "S050"))
  }

  if (group == "URINE_RCC"){
    # remove samples not meeting inclusion criterion
    tab <- filter(tab, ! tab$id %in% c("S25"))
  }

  tab$id <- paste0(group, "_",tab$id)
if (group == "RCC"){
  tab_all <- tab
} else{
  tab_all <- rbind(tab_all, tab)
}

}


tab_all <- tab_all %>% unique() %>%
  filter(!(group == "PLASMA_RCC" & grepl("S007|S011|S020|S015|S050", Sample))) %>%
  filter(!(group == "URINE_RCC" & grepl("S11|S12|S30|S33|S34|S48|S59|S60|
    S62|S64|S42|S25|S14", Sample))) %>% 
  filter(!(group == "PLASMA_CONTROL" & grepl("S040", Sample))) %>%
  filter(!(group == "URINE_CONTROL" & grepl("S56|S4.|S6.|S37", Sample))) %>%
  filter(!(grepl("JAN2020", group) & grepl(paste0(c("R104_AN", # <1M reads
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
    "R11_2353_.*L7", # extra lane
    "R12_2350_.*L7", # extra lane   
    "R24_LD021_.*L7", # extra lane  
    "R42_LD013_.*L7", # extra lane
    "R3_CS", # duplicate sample   
    "R2_CG", # extra coverage   
    "R6_2366_.*L7"), collapse="|"), # extra lane
    Sample))) 

# samples with fewer than 1M reads
tab_all %>% filter(type=="Unique") %>%
   bind_cols(tab_all %>% filter(type=="Duplicate") %>% 
   dplyr::rename(dup = count) %>% select(dup)) %>%
   mutate(total = dup + count) %>%
   mutate(ID=gsub("_CKD.*", "", Sample)) %>%
   filter(total < 1e6) %>%
   select(ID, group) %>% unique()

# samples using more than 1 lane
tab_all %>% filter(type=="Unique") %>%
   bind_cols(tab_all %>% filter(type=="Duplicate") %>% 
   dplyr::rename(dup = count) %>% select(dup)) %>%
   mutate(total = dup + count) %>%
  filter(grepl("val_1", Sample)) %>%
  mutate(ID=gsub("_CKD.*", "", Sample)) %>%
  group_by(ID, group) %>%
  summarize(count=n(), 
    tot_count=min(total)/1e6) %>%
  filter(count > 1) %>%
  select(ID, group, tot_count)


total <- tab_all %>% filter(type=="Unique") %>%
   bind_cols(tab_all %>% filter(type=="Duplicate") %>% 
   dplyr::rename(dup = count) %>% select(dup)) %>%
   mutate(total = dup + count) %>%
ggplot(aes(x = group, y = total, fill = group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Total number of reads") +
  labs(fill = "Sample Type") +
  xlab("Sample Type") 

tab_wide <- tab_all %>% filter(type=="Unique") %>%
   bind_cols(tab_all %>% filter(type=="Duplicate") %>% 
   dplyr::rename(dup = count) %>% select(dup)) %>%
   mutate(duprate = dup / (dup + count)) 
duprate <- ggplot(tab_wide,
  aes(x = group, y = duprate*100, fill = group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percentage of duplicate reads") +
  labs(fill = "Sample Type") +
  xlab("Sample Type") 


# bowtie
tool <- tools[2]

thisdir <- dirs[grepl(tool, dirs)]
thisdir <- thisdir[!grepl("BLCA|RCC|CONTROL|PRCA|PDAC|Novartis|JAN2020", thisdir)]

tab <- read_tsv(file.path(thisdir, "multiqc_general_stats.txt")) %>%
  mutate(group = ifelse(grepl("URINE_RCC", Sample), "URINE_RCC", 
                     ifelse(grepl("URINE_CONTROL", Sample), "URINE_CONTROL",
                      ifelse(grepl("RCCmet", Sample), "RCCmet",  
                        ifelse(grepl("BLCA", Sample), "BLCA", 
                         ifelse(grepl("PDAC", Sample), "PDAC", 
                          ifelse(grepl("PRCAcf", Sample), "PRCAcf", 
                            ifelse(grepl("Novartis", Sample), "Novartis",
                              ifelse(grepl("RCC", Sample), "RCC",
                                ifelse(grepl("PRCA", Sample), "PRCA", 
                                    "CONTROL")))))))))) %>%
  mutate(id = ifelse(group != "PDAC", 
    paste0(unlist(sapply(strsplit(tab$Sample, "_"), 
      function(x) x[[2]])), "_", group),
    paste0(unlist(sapply(strsplit(tab$Sample, "_"), 
      function(x) x[[3]])), "_", group))) %>%
  dplyr::rename(rate = `Bowtie 2_mqc-generalstats-bowtie_2-overall_alignment_rate`) %>%
  filter(group %in% c("RCC","CONTROL","URINE_RCC","URINE_CONTROL")) %>%
  filter(!(group == "RCC" & grepl("S007|S011|S020|S015|S050", Sample))) %>%
  filter(!(group == "URINE_RCC" & grepl("S11|S12|S30|S33|S34|S48|S59|S60|
    S62|S64|S42|S25|S14", Sample))) %>% 
  filter(!(group == "CONTROL" & grepl("S040", Sample))) %>%
  filter(!(group == "URINE_CONTROL" & grepl("S56|S4|S6|S37", Sample))) %>%
  filter(!(grepl("JAN2020", group) & grepl(paste0(c("R104_AN", # <1M reads
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
    "R3_CS", # duplicate sample
    "R2_CG", # extra coverage      
    "R6_2366_L7"), collapse="|"), # extra lane
    Sample))) %>%
  mutate(group = ifelse(group == "RCC", "PLASMA_RCC", 
    ifelse(group == "CONTROL", "PLASMA_CONTROL", group))) 

alignmentrate <- ggplot(tab, aes(x = group, y = rate, fill = group)) +
  geom_boxplot() +
  ylab("Overall Alignment percentage") +
  labs(fill = "Sample Type") +
  xlab("Sample Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


# jan2020 data 

thisdir <- dirs[grepl(tool, dirs)]
thisdir <- thisdir[grepl("JAN2020", thisdir)]

tab_jan2020 <- read_tsv(file.path(thisdir, "multiqc_general_stats.txt")) %>%
  mutate(ID = gsub("map_", "", Sample)) %>%
  mutate(ID = gsub("_L._JAN2020.err", "", ID)) %>%
  filter(!(grepl("JAN2020", group) & grepl(paste0(c("R104_AN", # <1M reads
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
    "R3_CS", # duplicate sample   
    "R2_CG", # extra coverage   
    "R42_LD013_L7", # extra lane   
    "R6_2366_L7"), collapse="|"), # extra lane
    Sample))) %>%
  left_join(meta2, by="ID") %>%
  dplyr::rename(rate = `Bowtie 2 / HiSAT2_mqc-generalstats-bowtie_2_hisat2-overall_alignment_rate`)

# no samples in jan2020 cohort with mapping rate below 25%
min(tab_jan2020$rate)

alignmentrate_jan2020 <- ggplot(tab_jan2020, 
  aes(x = Source, y = rate, fill = Status)) +
  geom_boxplot() +
  ylab("Overall Alignment percentage") +
  labs(fill = "Sample Type") +
  xlab("Sample Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

tab_jan2020 <- tab_jan2020 %>%
  mutate(group=paste0("JAN2020_", Source, "_", Status)) %>%
  rename(ID="id") %>%
  select(Sample, rate, group, id)
tab <- rbind(tab, tab_jan2020)


alignmentrate <- ggplot(tab, aes(x = group, y = rate, fill = group)) +
  geom_boxplot() +
  ylab("Overall Alignment percentage") +
  labs(fill = "Sample Type") +
  xlab("Sample Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


# get number of mapped reads

## num reads in bams (in house vs. toronto)
groups <- c("RCC", "CONTROL", "URINE_RCC", "URINE_CONTROL", "JAN2020")
names(groups) <- groups

tab <- NULL
for (g in seq_along(groups)){
  message("Counting group ", groups[g], "...")

  bams <- list.files(file.path(dat.dir, names(groups)[g]), 
      pattern = "sorted.bam$")
  
  if (groups[g] == "RCC"){
    meta.rcc <- read_excel("../../data/RCC/42 RCC fastq files.xlsx")
    ids.rcc <- gsub(".sorted.bam", "", bams)
    bams <- bams[ids.rcc %in% meta.rcc$`Fastq name`]
    bams <- bams[! bams %in% c("S007.sorted.bam", "S011.sorted.bam", "S020.sorted.bam", "S015.sorted.bam", "S050.sorted.bam")]
  }else if (groups[g] == "CONTROL"){
    bams <- bams[! bams %in% c("S040.sorted.bam")]
  }else if (groups[g] == "URINE_RCC"){
    bams <- bams[! bams %in% c("S56.sorted.bam", "S4.sorted.bam", "S6.sorted.bam", "S37.sorted.bam")]
  }else if (groups[g] == "URINE_CONTROL"){
    bams <- bams[! bams %in% c("S11.sorted.bam", "S12.sorted.bam", "S30.sorted.bam", "S33.sorted.bam", 
    "S34.sorted.bam", "S48.sorted.bam", "S59.sorted.bam", "S60.sorted.bam",
    "S62.sorted.bam", "S64.sorted.bam", "S42.sorted.bam", "S25.sorted.bam",
    "S14.sorted.bam")]
  }else if (groups[g] == "JAN2020"){
    bams <- bams[! grepl(paste0(c("R104_AN", # <1M reads
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
    "R3_CS", # duplicate sample    
    "R2_CG", # extra coverage   
    "R42_LD013_L7", # extra lane   
    "R6_2366_L7"), collapse="|"), # extra lane
    bams)]
  }

  for (b in bams){
    rw <- idxstatsBam(file.path(dat.dir, names(groups)[g], b[1]))
    rw <- rw %>% 
      summarize(mapped = sum(mapped)) %>%
      mutate(group = groups[g]) %>%
      mutate(group = ifelse(group == "RCC", "PLASMA_RCC", 
        ifelse(group == "CONTROL", "PLASMA_CONTROL", group))) %>%
      mutate(file=b)

    if (groups[g] == "JAN2020"){
    cntk <- colnames(rw)
    rw <- rw %>%
       mutate(ID = gsub("_L..sorted.bam", "", file)) %>%
       left_join(meta2, by="ID") %>%
       mutate(group = paste0(group, "_", Source, "_", Status)) %>%
       select(cntk)
    }

    tab <- rbind(tab, rw)
  }

}

mapped <- ggplot(tab, aes(x = group, y = mapped, fill = group)) +
  geom_boxplot() +
  ylab("Number of mapped reads") +
  labs(fill = "Sample Type") +
  xlab("Sample Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

## cowplot to combine
leg <- get_legend(total)
sf1 <- plot_grid(total + theme(legend.position = "none", 
                               axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank()),
                 duprate + theme(legend.position = "none",
                               axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank()),
                 alignmentrate + theme(legend.position = "none"),
                 mapped + theme(legend.position = "none"),
                 ncol = 2, align = "v", rel_heights=c(0.7,1),
                 labels=c("A", "B", "C", "D"))
plot_grid(sf1, leg, ncol=2, rel_widths = c(1,0.25))
ggsave(file.path(savedir, "SupplementaryFigure1_JAN2020.pdf"), width=12, height=8)

