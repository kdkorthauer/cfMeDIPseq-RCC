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
library(cowplot)
library(gtable)
library(gridGraphics)
library(grid)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(scatterplot3d)
library(MEDIPS)
library(grDevices)
library(Rtsne)

theme_set(theme_bw())

top <- 300
ws <- 300
merge <- FALSE
res <- NULL
savedir <- "../../out/MEDIPS_summary"
dir.create(savedir, showWarnings = FALSE)

# dir where binned medips objects of all samples are saved
outdir <- paste0("../../out/MEDIPS_", ws)


# read in spread sheet with stage/grade/histology
grade <- read_excel("../../data/RCC/Stage-Grade-Histology Analysis.xlsx")

meta <- read_excel("../../data/RCC/Keegan - RCC plasma.xlsx", 
  sheet = 1)

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

################# volcanoes

# plasma
plotVolcano <- function(diff.file, n1, n2, sig=0.1){
  diff <- readRDS(file=diff.file)
  dmrs <- diff[, grepl("adj.p|FC", colnames(diff))]
  dmrs <- dmrs[!is.na(dmrs$logFC),]

  #plot
  annotations <- data.frame(
        xpos = c(-Inf,Inf),
        ypos =  c(Inf,Inf),
        annotateText = c("Loss","Gain"),
        hjustvar = c(-0.8,1.5) ,
        vjustvar = c(2,2)) #<- adjust

  plt <- ggplot(dmrs, aes(x=logFC, y=-log10(limma.adj.p.value))) +
      stat_binhex(data=dmrs %>% filter(limma.adj.p.value > sig), bins=150) +
      #scale_fill_distiller(type="seq", direction = 1, palette="BuPu") +
      scale_fill_viridis_c(direction=-1) +
      geom_point(data=dmrs %>% filter(limma.adj.p.value < sig),
        color="red", size=0.6, alpha=0.4) +
      theme(legend.position = "none") +
      xlab("DNA Methylation (log2 Fold Change)") +
      ylab("Significance (-log10 q-value)") + 
      geom_hline(yintercept=-log10(sig)) +
      geom_vline(xintercept=0) +
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,
        vjust=vjustvar,label=annotateText))

  message(sum(dmrs$limma.adj.p.value < sig & dmrs$logFC < 0), " Lost")
  message(sum(dmrs$limma.adj.p.value < sig & dmrs$logFC > 0), " Gained")

  which.up <- which(dmrs$logFC > 0)
  which.down <- which(dmrs$logFC < 0)

  which.sig.up <- which(rank(dmrs$limma.adj.p.value[which.up], 
     ties.method = "random") <= as.numeric(top)/2)

  which.sig.down <- which(rank(dmrs$limma.adj.p.value[which.down], 
      ties.method = "random") <= as.numeric(top)/2)

  message(min(-log10(dmrs$limma.adj.p.value[which.down][which.sig.down])), " sig down")
  message(min(-log10(dmrs$limma.adj.p.value[which.up][which.sig.up])), " sig up")

  return(plt)
}


volcano_plasma <- plotVolcano(diff.file =file.path(outdir, "rcc.control.diff.rds"), 
  n1 = length(medip.rcc), 
  n2 = length(medip.control))

volcano_urine <- plotVolcano(diff.file =file.path(outdir, "urineR.urineC.diff.rds"), 
  n1 = length(medip.urineR), 
  n2 = length(medip.urineC), sig = 0.25)


################# heatmaps

# plasma
plotHeat <- function(diff.file, lab1, lab2, n1, n2, colnames = FALSE, 
	topSig = 1){
  diff <- readRDS(file=diff.file)
  which.up <- which(diff$logFC > 0)
  which.down <- which(diff$logFC < 0)

  if (topSig == 1){
    which.sig.up <- which(rank(diff$limma.adj.p.value[which.up], 
     ties.method = "random") <= as.numeric(top)/2)

    which.sig.down <- which(rank(diff$limma.adj.p.value[which.down], 
     ties.method = "random") <= as.numeric(top)/2)

    which.sig <- c(which.up[which.sig.up], which.down[which.sig.down])
  }else{
  	which.sig <- which(diff$limma.adj.p.value < topSig)
  }

  dmrs <- diff[which.sig, grepl("counts", colnames(diff))]
  dmrs <- dmrs[,!grepl("mean", colnames(dmrs))]

  rownames(dmrs) <- paste0(diff[which.sig,]$chr, ":", 
   diff[which.sig,]$start, "-",
   diff[which.sig,]$end)
  colnames(dmrs)[1:n1] <- paste0(lab1, "_", colnames(dmrs)[1:n1])
  colnames(dmrs)[(n1 + 1):(n1 + n2)] <- paste0(lab2, "_", colnames(dmrs)[(n1 + 1):(n1 + n2)])
  colnames(dmrs) <- gsub(".counts|.rpkm", "", colnames(dmrs))
  colnames(dmrs) <- gsub(".sorted.bam", "", colnames(dmrs))
  colnames(dmrs) <- gsub("\\.1", "", colnames(dmrs))

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


  ht = Heatmap(log(dmrs+1), name = "log(CPM+1)", 
   top_annotation = ha_column, col = ecolors,
   show_row_names = FALSE, show_column_names = colnames,
   column_names_gp = gpar(fontsize = 9),
   heatmap_legend_param = list(legend_direction="horizontal"))

  ht

}

# top 300

heat_plasma <- plotHeat(diff.file =file.path(outdir, "rcc.control.diff.rds"), 
  lab1="rcc", lab2="control",
  n1 = length(medip.rcc), 
  n2 = length(medip.control))

heat_urine <- plotHeat(diff.file =file.path(outdir, "urineR.urineC.diff.rds"), 
  lab1="urineR", lab2="urineC",
  n1 = length(medip.urineR), 
  n2 = length(medip.urineC))


## all significant

heat_urine_all <- plotHeat(diff.file =file.path(outdir, "urineR.urineC.diff.rds"), 
  lab1="urineR", lab2="urineC",
  n1 = length(medip.urineR), 
  n2 = length(medip.urineC), 
  topSig = 0.25)

heat_plasma_all <- plotHeat(diff.file =file.path(outdir, "rcc.control.diff.rds"), 
  lab1="rcc", lab2="control",
  n1 = length(medip.rcc), 
  n2 = length(medip.control),
  topSig = 0.05)



################# AUC

# file list
files <- list.files(outdir, pattern = "*.txt", recursive = TRUE, 
 full.names = TRUE)
files <- files[grepl(paste0("top", top), files)]
files <- files[!grepl(paste0("validation"), files)]
files_blca <- files[grepl("blca",files)]
files <- files[!grepl("blca",files)]
files <- files[!grepl("v",files)]
files <- files[!grepl("PC_",files)]
files <- files[!grepl("sampleprob",files)]

tmp <- files %>%
purrr::map(read_tsv) %>%
do.call("rbind", .)

res <- rbind(res, tmp)

res <- res %>% filter(lab1 != "rcc_D") %>%
  mutate(method = as.factor(method),
  	     window = as.factor(window),
  	     type = ifelse(lab1=="rcc", "plasma", "urine")) %>%
  filter(method == "original")


# summary table
res_summary <- res %>% 
  group_by(lab1, lab2, window, method) %>%
  summarize(meanAUC = mean(auc, na.rm = TRUE),
            sdAUC = sd(auc, na.rm = TRUE),
  	        n = n(),
            lowerAUC = meanAUC - qnorm(0.975)*sdAUC/sqrt(n),
            upperAUC = meanAUC + qnorm(0.975)*sdAUC/sqrt(n))

res_summary %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  select(lab1, lab2, window, method, 
         meanAUC, lowerAUC, upperAUC) %>%
  write.table(quote=FALSE, row.names=FALSE,
  file=file.path(savedir, 
    paste0("auc_table_top",top, "_ws", ws, "_original.txt")), sep = "\t")

# plot of AUC

auc_plasma <- res %>% filter(type=="plasma") %>%
  select(type, auc) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_boxplot(outlier.shape = NA,fill =  "#56B4E9") +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") 

auc_urine <- res %>% filter(type=="urine") %>%
  select(type, auc) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_boxplot(outlier.shape = NA, fill =  "#56B4E9") +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1)+
  xlab("")+ylab("AUC") 


auc_plasma_violin <- res %>% filter(type=="plasma") %>%
  select(type, auc) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_violin(fill =  "#56B4E9") +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") 
ggsave(auc_plasma_violin, file=file.path(savedir, "auc_plasma_violin.pdf"),
  width=3, height = 2.5)

auc_urine_violin <- res %>% filter(type=="urine") %>%
  select(type, auc) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_violin(fill =  "#56B4E9") +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1)+
  xlab("")+ylab("AUC") 
ggsave(auc_urine_violin, file=file.path(savedir, "auc_urine_violin.pdf"),
  width=3, height = 2.5)


############ Cowplot to combine

# figure 1
col1 <- plot_grid(volcano_plasma, grob(), 
     plot_grid(grob(), auc_plasma, grob(), rel_widths=c(0.3,1,0.3), nrow=1), 
  ncol=1, rel_heights = c(1, 0.05, 0.5), rel_widths = c(1, 0.5),
  labels = c("A", "C"))

plot_grid(col1, grid.grabExpr(draw(heat_plasma,
    annotation_legend_side = "right", 
  heatmap_legend_side="bottom"), wrap=TRUE), 
  rel_widths = c(1,1.7), labels=c("", "B"))
ggsave(file.path(savedir, "Figure1.pdf"), width=10, height=6)

# figure 2
col1 <- plot_grid(volcano_urine, grob(),
     plot_grid(grob(), auc_urine, grob(), rel_widths=c(0.3,1,0.3), nrow=1), 
  ncol=1, rel_heights = c(1, 0.05, 0.5), rel_widths = c(1, 0.5),
  labels = c("A", "C"))

plot_grid(col1, grid.grabExpr(draw(heat_urine,
  annotation_legend_side = "right", 
  heatmap_legend_side="bottom"), wrap=TRUE), 
  rel_widths = c(1,1.7), labels=c("", "B"))
ggsave(file.path(savedir, "Figure2.pdf"), width=10, height=6)


# figure 1
col1 <- plot_grid(volcano_plasma, grob(), 
     plot_grid(grob(), auc_plasma, grob(), rel_widths=c(0.3,1,0.3), nrow=1), 
  ncol=1, rel_heights = c(1, 0.05, 0.5), rel_widths = c(1, 0.5),
  labels = c("A", "C"))

plot_grid(col1, grid.grabExpr(draw(heat_plasma_all,
    annotation_legend_side = "right", 
  heatmap_legend_side="bottom"), wrap=TRUE), 
  rel_widths = c(1,1.7), labels=c("", "B"))
ggsave(file.path(savedir, "Figure1_allheat.pdf"), width=10, height=6)

# figure 2
col1 <- plot_grid(volcano_urine, grob(),
     plot_grid(grob(), auc_urine, grob(), rel_widths=c(0.3,1,0.3), nrow=1), 
  ncol=1, rel_heights = c(1, 0.05, 0.5), rel_widths = c(1, 0.5),
  labels = c("A", "C"))

plot_grid(col1, grid.grabExpr(draw(heat_urine_all,
  annotation_legend_side = "right", 
  heatmap_legend_side="bottom"), wrap=TRUE), 
  rel_widths = c(1,1.7), labels=c("", "B"))
ggsave(file.path(savedir, "Figure2_allheat.pdf"), width=10, height=6)


############ Supplementary

## unique vs duplicate reads for plasma & urine

outdir <- "../../out/multiqc"
dat.dir <- "../../out/sortedbam_dup"
setwd(outdir)

groups <- c("RCC", "CONTROL", "URINE_RCC", "URINE_CONTROL")

tools <- c("fastqc_trimmed", "botwie")

dirs <- list.dirs()

tab_all <- NULL

meta <- read_excel("../../data/RCC/42 RCC fastq files.xlsx")


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
  filter(!(group == "URINE_CONTROL" & grepl("S56|S4.|S6.|S37", Sample))) 


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
thisdir <- thisdir[!grepl("BLCA|RCC|CONTROL|PRCA|PDAC|Novartis", thisdir)]

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
  mutate(group = ifelse(group == "RCC", "PLASMA_RCC", 
    ifelse(group == "CONTROL", "PLASMA_CONTROL", group))) 

alignmentrate <- ggplot(tab, aes(x = group, y = rate, fill = group)) +
  geom_boxplot() +
  ylab("Overall Alignment percentage") +
  labs(fill = "Sample Type") +
  xlab("Sample Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


# get number of mapped reads

## num reads in bams (in house vs. toronto)
groups <- c("RCC", "CONTROL", "URINE_RCC", "URINE_CONTROL")
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
  }else{
    bams <- bams[! bams %in% c("S11.sorted.bam", "S12.sorted.bam", "S30.sorted.bam", "S33.sorted.bam", 
    "S34.sorted.bam", "S48.sorted.bam", "S59.sorted.bam", "S60.sorted.bam",
    "S62.sorted.bam", "S64.sorted.bam", "S42.sorted.bam", "S25.sorted.bam",
    "S14.sorted.bam")]
  }

  for (b in bams){
    rw <- idxstatsBam(file.path(dat.dir, names(groups)[g], b[1]))
    rw <- rw %>% 
      summarize(mapped = sum(mapped)) %>%
      mutate(group = groups[g]) %>%
      mutate(group = ifelse(group == "RCC", "PLASMA_RCC", 
        ifelse(group == "CONTROL", "PLASMA_CONTROL", group)))
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
ggsave(file.path(savedir, "SupplementaryFigure1.pdf"), width=9, height=7)


### AUC RCC vx BLCA
res <- NULL
tmp <- files_blca %>%
purrr::map(read_tsv) %>%
do.call("rbind", .)

res <- rbind(res, tmp)

res <- res %>%
  mutate(method = as.factor(method),
  	     window = as.factor(window),
  	     type = ifelse(lab1=="rcc", "plasma", "urine")) %>%
  filter(method == "original")


# summary table
res_summary <- res %>% 
  group_by(lab1, lab2, window, method) %>%
  summarize(meanAUC = mean(auc, na.rm = TRUE),
            sdAUC = sd(auc, na.rm = TRUE),
  	        n = n(),
            lowerAUC = meanAUC - qnorm(0.975)*sdAUC/sqrt(n),
            upperAUC = meanAUC + qnorm(0.975)*sdAUC/sqrt(n))

res_summary %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  select(lab1, lab2, window, method, 
         meanAUC, lowerAUC, upperAUC) %>%
  write.table(quote=FALSE, row.names=FALSE,
  file=file.path(savedir, 
    paste0("auc_table_top",top, "_ws", ws, "_original_blca.txt")), sep = "\t")

# plot of AUC

auc_plasma <- res %>% filter(type=="plasma") %>%
  select(type, auc) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_boxplot(outlier.shape = NA,fill =  "#56B4E9") +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") 

auc_plasma
ggsave(file.path(savedir, "AUC_rcc_bcla.pdf"), width=3, height=2.75)


auc_plasma <- res %>% filter(type=="plasma") %>%
  select(type, auc) %>%
  ggplot(aes(x = type, y = auc)) +
  geom_violin(fill =  "#56B4E9") +
  geom_jitter(height = 0, alpha = 0.3)+
  ylim(0,1) +
  xlab("")+ylab("AUC") 
auc_plasma
ggsave(file.path(savedir, "auc_rcc_bcla_violin.pdf"), width=3, height=2.5)



### PCA plots - all sig

depths <- function(mdobjlist, CS, type){
  depth <- data.frame(sapply(mdobjlist, function(x){
  	  x@genome_count[CS@genome_CF >= 0]
  }))
  colnames(depth) <- sapply(mdobjlist, function(x){ 
  	gsub(".bam|.sorted.bam", "", x@sample_name)})
  colnames(depth) <- gsub("_", "", colnames(depth))
  colnames(depth) <- paste0(type, "_", colnames(depth))
  return(depth)
}

CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip.rcc[[1]])

if (!file.exists(file.path(outdir, "PCA_1_2_3_plasma_all.pdf"))){
  ### plasma
  df <- cbind(depths(medip.rcc, CS, "rcc"),
  	depths(medip.control, CS, "ctrl"))
  # include all sig rows
  diff.file =file.path(outdir, "rcc.control.diff.rds")
  diff <- readRDS(file=diff.file)
  which.sig <- which(diff$limma.adj.p.value < 0.05)
  rm(diff)
  df <- df[which.sig,] 

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
           id = ids) %>%
    mutate(Type = ifelse(type == "rcc", "RCC", "Control"))
  
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)

  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = Type), size = 2)  +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_1_2_plasma_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_2_3_plasma_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC1, y=PC3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_1_3_plasma_all.pdf"), width = 4.5, height = 3.5)

  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
  colors <- colors[as.numeric(as.factor(tidydf$Type))]
  pdf(file.path(outdir, "PCA_1_2_3_plasma_all.pdf"), width = 4.5, height = 4.5)
   scatterplot3d(tidydf[,1:3], pch =20, 
  	 xlab="PC1", ylab="PC2", zlab="PC3", cex.symbols = 2,
  	 color=colors)
   legend(0,-3.75, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
  dev.off()

  write.table(data.frame(PC=1:10, Proportion=(pcs$sdev/sum(pcs$sdev))[1:10]), 
      quote=FALSE, row.names=FALSE,
      file=file.path(outdir, paste0("PC_proportionVariation_plasma_all.txt")), 
      sep = "\t")

  # tsne
  tsne <- Rtsne(t(log(df+1)), dims = 3, perplexity = 10)
  tidydf <- select(data.frame(tsne$Y), "X1", "X2", "X3") %>%
    mutate(tSNE1 = X1,
           tSNE2 = X2,
           tSNE3 = X3) %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) %>%
    mutate(Type = ifelse(type == "rcc", "RCC", "Control"))
  
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)
 
  ggplot() +
    geom_point(data = tidydf, aes(x=tSNE1, y=tSNE2, colour = Type), size = 2)  +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_2_plasma_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE2, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_2_3_plasma_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE1, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_3_plasma_all.pdf"), width = 4.5, height = 3.5)

  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
  colors <- colors[as.numeric(as.factor(tidydf$Type))]
  pdf(file.path(outdir, "tSNE_1_2_3_plasma_all.pdf"), width = 4.5, height = 4.5)
   scatterplot3d(tidydf[,1:3], pch =20, 
     xlab="tSNE1", ylab="tSNE2", zlab="tSNE3", cex.symbols = 2,
     color=colors)
   legend(0,-4.95, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
  dev.off()


 }


if (!file.exists(file.path(outdir, "PCA_1_2_3_urine_all.pdf"))){
  ### urine
  df <- cbind(depths(medip.urineR, CS, "urineR"),
  	depths(medip.urineC, CS, "urineC"))
  # include all sig rows
  diff.file =file.path(outdir, "urineR.urineC.diff.rds")
  diff <- readRDS(file=diff.file)
  which.sig <- which(diff$limma.adj.p.value < 0.25)
  rm(diff)
  df <- df[which.sig,] 

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
           id = ids) %>%
    mutate(Type = ifelse(type == "urineR", "RCC", "Control"))
  
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)

  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = Type), size = 2)  +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_1_2_urine_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = Type)) +
    geom_point(size = 2)+
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_2_3_urine_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC1, y=PC3, colour = Type)) +
    geom_point(size = 2)+
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_1_3_urine_all.pdf"), width = 4.5, height = 3.5)


  pdf(file.path(outdir, "PCA_1_2_3_urine_all.pdf"), width = 4.5, height = 4.5)
   colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
   colors <- colors[as.numeric(as.factor(tidydf$Type))]
   scatterplot3d(tidydf[,1:3], pch =20, 
  	 xlab="PC1", ylab="PC2", zlab="PC3", cex.symbols = 2,
  	 color=colors)
   legend(0,-6.75, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
  dev.off()


  write.table(data.frame(PC=1:10, Proportion=(pcs$sdev/sum(pcs$sdev))[1:10]), 
      quote=FALSE, row.names=FALSE,
      file=file.path(outdir, paste0("PC_proportionVariation_urine_all.txt")), 
      sep = "\t")

  tsne <- Rtsne(t(log(df+1)), dims = 3, perplexity = 10)
  tidydf <- select(data.frame(tsne$Y), "X1", "X2", "X3") %>%
    mutate(tSNE1 = X1,
           tSNE2 = X2,
           tSNE3 = X3) %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) %>%
    mutate(Type = ifelse(type == "urineR", "RCC", "Control"))
  
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)
 
  ggplot() +
    geom_point(data = tidydf, aes(x=tSNE1, y=tSNE2, colour = Type), size = 2)  +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_2_urine_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE2, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_2_3_urine_all.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE1, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_3_urine_all.pdf"), width = 4.5, height = 3.5)

  pdf(file.path(outdir, "tSNE_1_2_3_urine_all.pdf"), width = 4.5, height = 4.5)
   colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
   colors <- colors[as.numeric(as.factor(tidydf$Type))]
   scatterplot3d(tidydf[,1:3], pch =20, 
     xlab="tSNE1", ylab="tSNE2", zlab="tSNE3", cex.symbols = 2,
     color=colors)
   legend(0,-4, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
  dev.off()

}



### PCA plots - top300 - don't normalize on just 300 genes

if (!file.exists(file.path(outdir, "PCA_1_2_3_plasma_top300.pdf"))){
  ### plasma
  df <- cbind(depths(medip.rcc, CS, "rcc"),
  	depths(medip.control, CS, "ctrl"))
  # include all sig rows
  diff.file =file.path(outdir, "rcc.control.diff.rds")
  diff <- readRDS(file=diff.file)
  which.up <- which(diff$logFC > 0)
  which.down <- which(diff$logFC < 0)
  which.sig.up <- which(rank(diff$limma.adj.p.value[which.up], 
     ties.method = "random") <= as.numeric(top)/2)

  which.sig.down <- which(rank(diff$limma.adj.p.value[which.down], 
     ties.method = "random") <= as.numeric(top)/2)

  which.sig <- c(which.up[which.sig.up], which.down[which.sig.down])
 
  rm(diff)
  df <- df[which.sig,] 

  grp <- gsub("_.*", "", colnames(df))
  ids = unlist(sapply(strsplit(colnames(df), "_"), function(x) x[[2]]))
  x <- match(ids, meta$`Sample number`)
  subtype <- ifelse(grp == "rcc" & ids %in% meta$`Sample number`,
    meta$Histology[x], NA)

 
  pcs <- Morpho::prcompfast(t(log(df+1)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) %>%
    mutate(Type = ifelse(type == "rcc", "RCC", "Control"))
    
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)

  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = Type), size = 2) +
    scale_color_manual(values = colors) 
  ggsave(file.path(outdir, "PCA_1_2_plasma_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = Type)) +
    geom_point(size = 2)+
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_2_3_plasma_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC1, y=PC3, colour = Type)) +
    geom_point(size = 2)+
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_1_3_plasma_top300.pdf"), width = 4.5, height = 3.5)

  pdf(file.path(outdir, "PCA_1_2_3_plasma_top300.pdf"), width = 4.5, height = 4.5)
   colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
   colors <- colors[as.numeric(as.factor(tidydf$Type))]
   scatterplot3d(tidydf[,1:3], pch =20, 
  	 xlab="PC1", ylab="PC2", zlab="PC3", cex.symbols = 2,
  	 color=colors)
   legend(0,-5.3, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
   dev.off()


  write.table(data.frame(PC=1:10, Proportion=(pcs$sdev/sum(pcs$sdev))[1:10]), 
      quote=FALSE, row.names=FALSE,
      file=file.path(outdir, paste0("PC_proportionVariation_plasma_top300.txt")), 
      sep = "\t")

  tsne <- Rtsne(t(log(df+1)), dims = 3, perplexity = 10)
  tidydf <- select(data.frame(tsne$Y), "X1", "X2", "X3") %>%
    mutate(tSNE1 = X1,
           tSNE2 = X2,
           tSNE3 = X3) %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) %>%
    mutate(Type = ifelse(type == "rcc", "RCC", "Control"))
  
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)
 
  ggplot() +
    geom_point(data = tidydf, aes(x=tSNE1, y=tSNE2, colour = Type), size = 2)  +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_2_plasma_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE2, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_2_3_plasma_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE1, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_3_plasma_top300.pdf"), width = 4.5, height = 3.5)

  pdf(file.path(outdir, "tSNE_1_2_3_plasma_top300.pdf"), width = 4.5, height = 4.5)
   colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
   colors <- colors[as.numeric(as.factor(tidydf$Type))]
   scatterplot3d(tidydf[,1:3], pch =20, 
     xlab="tSNE1", ylab="tSNE2", zlab="tSNE3", cex.symbols = 2,
     color=colors)
   legend(0,-5.3, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
   dev.off()

}


if (!file.exists(file.path(outdir, "PCA_1_2_3_urine_top300.pdf"))){
  ### urine
  df <- cbind(depths(medip.urineR, CS, "urineR"),
  	depths(medip.urineC, CS, "urineC"))
  # include all sig rows
  diff.file =file.path(outdir, "urineR.urineC.diff.rds")
  diff <- readRDS(file=diff.file)
  which.up <- which(diff$logFC > 0)
  which.down <- which(diff$logFC < 0)
  which.sig.up <- which(rank(diff$limma.adj.p.value[which.up], 
     ties.method = "random") <= as.numeric(top)/2)

  which.sig.down <- which(rank(diff$limma.adj.p.value[which.down], 
     ties.method = "random") <= as.numeric(top)/2)

  which.sig <- c(which.up[which.sig.up], which.down[which.sig.down])
  rm(diff)
  df <- df[which.sig,] 

  grp <- gsub("_.*", "", colnames(df))
  ids = unlist(sapply(strsplit(colnames(df), "_"), function(x) x[[2]]))
  x <- match(ids, meta$`Sample number`)
  subtype <- ifelse(grp == "rcc" & ids %in% meta$`Sample number`,
    meta$Histology[x], NA)

 
  pcs <- Morpho::prcompfast(t(log(df+1)), center = TRUE, scale. = TRUE)
  tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) %>%
    mutate(Type = ifelse(type == "urineR", "RCC", "Control"))
    
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)

  ggplot() +
    geom_point(data = tidydf, aes(x=PC1, y=PC2, colour = Type), size = 2) +
    scale_color_manual(values = colors) 
  ggsave(file.path(outdir, "PCA_1_2_urine_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC2, y=PC3, colour = Type)) +
    geom_point(size = 2)+
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_2_3_urine_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=PC1, y=PC3, colour = Type)) +
    geom_point(size = 2)+
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "PCA_1_3_urine_top300.pdf"), width = 4.5, height = 3.5)

  pdf(file.path(outdir, "PCA_1_2_3_urine_top300.pdf"), width = 4.5, height = 4.5)
   colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
   colors <- colors[as.numeric(as.factor(tidydf$Type))]
   scatterplot3d(tidydf[,1:3], pch =20, 
  	 xlab="PC1", ylab="PC2", zlab="PC3", cex.symbols = 2,
  	 color=colors)
   legend(0,-3.92, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
  dev.off()


  write.table(data.frame(PC=1:10, Proportion=(pcs$sdev/sum(pcs$sdev))[1:10]), 
      quote=FALSE, row.names=FALSE,
      file=file.path(outdir, paste0("PC_proportionVariation_urine_top300.txt")), 
      sep = "\t")

  tsne <- Rtsne(t(log(df+1)), dims = 3, perplexity = 10)
  tidydf <- select(data.frame(tsne$Y), "X1", "X2", "X3") %>%
    mutate(tSNE1 = X1,
           tSNE2 = X2,
           tSNE3 = X3) %>%
    mutate(type = grp,
           subtype = subtype,
           id = ids) %>%
    mutate(Type = ifelse(type == "urineR", "RCC", "Control"))
  
  colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.8)
 
  ggplot() +
    geom_point(data = tidydf, aes(x=tSNE1, y=tSNE2, colour = Type), size = 2)  +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_2_urine_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE2, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_2_3_urine_top300.pdf"), width = 4.5, height = 3.5)

  tidydf %>%ggplot(aes(x=tSNE1, y=tSNE3, colour = Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors)
  ggsave(file.path(outdir, "tSNE_1_3_urine_top300.pdf"), width = 4.5, height = 3.5)

  pdf(file.path(outdir, "tSNE_1_2_3_urine_top300.pdf"), width = 4.5, height = 4.5)
   colors <-  adjustcolor(c("#E69F00", "#56B4E9"), alpha=0.5)
   colors <- colors[as.numeric(as.factor(tidydf$Type))]
   scatterplot3d(tidydf[,1:3], pch =20, 
     xlab="tSNE1", ylab="tSNE2", zlab="tSNE3", cex.symbols = 2,
     color=colors)
   legend(0,-3.8, legend = levels(as.factor(tidydf$Type)),
      col =  c("#E69F00", "#56B4E9"), pch = 20, 
      inset = -0.25, xpd = TRUE, horiz = TRUE, bty="n")
  dev.off()
}

### boxplot of sample probabilities for those in test set

# list files
files <- list.files(pattern = "sampleprob*", recursive = TRUE,
	path = outdir, full.names = TRUE)
files <- files[!grepl("pdf", files)]
files <- files[!grepl("v", files)]

tmp <- files %>%
       purrr::map(read_tsv) %>%
       do.call("rbind", .)

grade$true_label <- "RCC"
plasma <- tmp %>% filter(!grepl("urine", true_label)) %>%
  mutate(true_label = ifelse(true_label=="rcc", "RCC", "Control")) %>%
  mutate(Sample = gsub("rcc_|control_", "", sample_name)) %>%
  left_join(grade, by = c("Sample", "true_label"))
urine <- tmp %>% filter(grepl("urine", true_label))%>%
  mutate(true_label = ifelse(true_label=="urineR", "RCC", "Control")) %>%
  mutate(Sample = gsub("urineR_|urineC_", "", sample_name)) %>%
  left_join(grade, by = c("Sample", "true_label"))
cols <- c("RCC" = "#56B4E9", "Control" = "#E69F00")


plasma <- plasma %>%
  mutate(Stage = ifelse(true_label == "RCC", Stage, "Control"),
    Grade = ifelse(true_label == "RCC", Grade, "Control"),
    Histology = ifelse(true_label == "RCC", Histology, "Control"))

urine <- urine %>%
  mutate(Stage = ifelse(true_label == "RCC", Stage, "Control"),
    Grade = ifelse(true_label == "RCC", Grade, "Control"),
    Histology = ifelse(true_label == "RCC", Histology, "Control"))

plasma$Stage <- as.factor(plasma$Stage)
urine$Stage <- as.factor(urine$Stage)

plasma$Grade <- as.factor(plasma$Grade)
urine$Grade <- as.factor(urine$Grade)

plasma <- mutate(plasma, 
  Histology = ifelse(grepl("Clear", Histology), "Clear Cell", Histology))
plasma$Histology <- as.factor(plasma$Histology)
plasma$Histology <- factor(plasma$Histology, levels = levels(plasma$Histology)[c(1,3,2)])


# each sample

plasma %>% ggplot(aes(x = sample_name, y = class_prob, fill = true_label)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
  labs(fill = "True class") + 
  scale_fill_manual(values = cols) + 
  ylab("Test set probability of RCC") + xlab("Sample")
ggsave(file.path(outdir, "Boxplot_testset_probs_plasma.pdf"),
  width = 5, height = 3)

plasma %>% ggplot(aes(x = sample_name, y = class_prob, fill = true_label)) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = "True class") + 
  scale_fill_manual(values = cols) + 
  ylab("Test set probability of RCC") + xlab("Sample")
ggsave(file.path(outdir, "Boxplot_testset_probs_plasma_nolabels.pdf"),
  width = 5, height = 2.7)

# group by stage/grade/histo
plasma %>% filter(!is.na(Stage)) %>%
  ggplot(aes(x = Stage, y = class_prob)) +
  geom_boxplot() + 
  ylab("Test set probability of RCC") + xlab("Stage")
ggsave(file.path(outdir, "Boxplot_testset_probs_plasma_stage.pdf"),
  width = 4, height = 3)

plasma %>% filter(!is.na(Grade)) %>%
  ggplot(aes(x = Grade, y = class_prob)) +
  geom_boxplot() +
  ylab("Test set probability of RCC") + xlab("Grade")
ggsave(file.path(outdir, "Boxplot_testset_probs_plasma_grade.pdf"),
  width = 4, height = 3)

plasma %>% filter(!is.na(Histology)) %>%
  ggplot(aes(x = Histology, y = class_prob)) +
  geom_boxplot() +
  ylab("Test set probability of RCC") + xlab("Histology")
ggsave(file.path(outdir, "Boxplot_testset_probs_plasma_hist.pdf"),
  width = 4, height = 3)


plasma %>% filter(!is.na(Stage)) %>%
  ggplot(aes(x = Stage, y = class_prob)) +
  geom_violin() + 
  ylab("Test set probability of RCC") + xlab("Stage")
ggsave(file.path(outdir, "Violin_testset_probs_plasma_stage.pdf"),
  width = 4, height = 3)

plasma %>% filter(!is.na(Grade)) %>%
  ggplot(aes(x = Grade, y = class_prob)) +
  geom_violin() +
  ylab("Test set probability of RCC") + xlab("Grade")
ggsave(file.path(outdir, "Violin_testset_probs_plasma_grade.pdf"),
  width = 4, height = 3)

plasma %>% filter(!is.na(Histology)) %>%
  ggplot(aes(x = Histology, y = class_prob)) +
  geom_violin() +
  ylab("Test set probability of RCC") + xlab("Histology")
ggsave(file.path(outdir, "Violin_testset_probs_plasma_hist.pdf"),
  width = 4, height = 3)

# each sample

urine %>% ggplot(aes(x = sample_name, y = class_prob, fill = true_label)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
  labs(fill = "True class") + 
  scale_fill_manual(values = cols) + 
  ylab("Test set probability of RCC") + xlab("Sample")
ggsave(file.path(outdir, "Boxplot_testset_probs_urine.pdf"),
  width = 5, height = 3)

urine %>% ggplot(aes(x = sample_name, y = class_prob, fill = true_label)) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = "True class") + 
  scale_fill_manual(values = cols) + 
  ylab("Test set probability of RCC") + xlab("Sample")
ggsave(file.path(outdir, "Boxplot_testset_probs_urine_nolabels.pdf"),
  width = 5, height = 2.7)


# group by stage/grade/histo
urine %>% filter(!is.na(Stage)) %>%
  ggplot(aes(x = Stage, y = class_prob)) +
  geom_boxplot() + 
  ylab("Test set probability of RCC") + xlab("Stage")
ggsave(file.path(outdir, "Boxplot_testset_probs_urine_stage.pdf"),
  width = 4, height = 3)

urine %>% filter(!is.na(Grade)) %>%
  ggplot(aes(x = Grade, y = class_prob)) +
  geom_boxplot() +
  ylab("Test set probability of RCC") + xlab("Grade")
ggsave(file.path(outdir, "Boxplot_testset_probs_urine_grade.pdf"),
  width = 4, height = 3)


urine %>% filter(!is.na(Stage)) %>%
  ggplot(aes(x = Stage, y = class_prob)) +
  geom_violin() + 
  ylab("Test set probability of RCC") + xlab("Stage")
ggsave(file.path(outdir, "Violin_testset_probs_urine_stage.pdf"),
  width = 4, height = 3)

urine %>% filter(!is.na(Grade)) %>%
  ggplot(aes(x = Grade, y = class_prob)) +
  geom_violin() +
  ylab("Test set probability of RCC") + xlab("Grade")
ggsave(file.path(outdir, "Violin_testset_probs_urine_grade.pdf"),
  width = 4, height = 3)

## repeat above but sorted

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
    df_av <- df %>% group_by(oldFactor) %>% summarise(meanSortingVariable = median(sortingVariable, na.rm=TRUE)) %>% 
      arrange(meanSortingVariable)
  }
  
  # Compute average of sortingVariable and arrange (descending)
  if (ascending == FALSE) {
    df_av <- df %>% group_by(oldFactor) %>% summarise(meanSortingVariable = median(sortingVariable, na.rm=TRUE)) %>% 
      arrange(desc(meanSortingVariable))
  }
  
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = df_av$oldFactor)
  return(newFactor)
}

# Sort factor levels by their frequency of occurrence
sortLvlsByN.fnc <- function(oldFactor, ascending = TRUE) {

  require("magrittr")
  
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = table(oldFactor)  %>% sort(., decreasing = !ascending)  %>% names())
  return(newFactor)
}

# Sort factor levels arbitrarily
sortLvls.fnc <- function(oldFactor, levelOrder) {
  if(!is.factor(oldFactor)) stop("The variable you want to reorder isn't a factor.")
  
  if(!is.numeric(levelOrder)) stop("'order' should be a numeric vector.")
  
  if(max(levelOrder) > length(levels(oldFactor))) stop("The largest number in 'order' can't be larger than the number of levels in the factor.")
  
  if(length(levelOrder) > length(levels(oldFactor))) stop("You can't have more elements in 'order' than there are levels in the factor.")
  
  if(length(levelOrder) == length(levels(oldFactor))) {
    reorderedFactor <- factor(oldFactor, levels = levels(oldFactor)[levelOrder])
  }
  
  if(length(levelOrder) < length(levels(oldFactor))) {
    levelOrderAll <- c(levelOrder, (1:length(levels(oldFactor)))[-levelOrder])
    reorderedFactor <- factor(oldFactor, levels = levels(oldFactor)[levelOrderAll])
  }
  
  return(reorderedFactor)
}
####  end functions to reorder

plasma$sample_name <- as.factor(plasma$sample_name)
plasma$sample_name <- sortLvlsByVar.fnc(plasma$sample_name, plasma$class_prob)

urine$sample_name <- as.factor(urine$sample_name)
urine$sample_name <- sortLvlsByVar.fnc(urine$sample_name, urine$class_prob)


plasma %>% ggplot(aes(x = sample_name, y = class_prob, fill = true_label)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
  labs(fill = "True class") + 
  scale_fill_manual(values = cols) + 
  ylab("Test set probability of RCC") + xlab("Sample")
ggsave(file.path(outdir, "Boxplot_testset_probs_plasma_sorted.pdf"),
  width = 5, height = 3)

urine %>% ggplot(aes(x = sample_name, y = class_prob, fill = true_label)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
  labs(fill = "True class") + 
  scale_fill_manual(values = cols) + 
  ylab("Test set probability of RCC") + xlab("Sample")
ggsave(file.path(outdir, "Boxplot_testset_probs_urine_sorted.pdf"),
  width = 5, height = 3)