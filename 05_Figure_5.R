#Figure 5 plots
#save all necessary files if needed with ggsave
#bulkseq from GSE135251
set.seed(1234567)
setwd("C:/Users/u0125188/Desktop/Liver") ##change to output directory
library(Seurat)
library(SeuratObject)
library(dplyr)
library (ggplot2)
library(cowplot)
library(clustree)
library(patchwork)
library(pheatmap)
library(Matrix)
library(readxl)
library(stringr)
library(SingleCellExperiment)
library(RColorBrewer)
library(biomaRt)
library(openxlsx)
library(ggpubr)
library("xlsx")
library(openxlsx)
library(gprofiler2)
library(ggpubr)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(GSEABase)
library(limma)
library(pheatmap)
library(dendsort)

#setup
fontsize<-14.5
colours_disease<-c("palegreen4", "steelblue3", "sienna3", "red4" )
colours_Liver<-c("lightpink3" , "slategray3","mediumpurple2","navy","#00B6EB","turquoise4", "burlywood","#FB61D7")
colours_Hepato_Chol<-c("deeppink4", "orange3", "darkslategrey", "turquoise3")
colours_Mesenchymal<-c("deeppink4", "orange3", "darkslategrey")
colours_Lymphocyte<-c("slategray3","lightpink3","navy","darkgreen","darkmagenta","mediumorchid2","slateblue2","dodgerblue3","#00B6EB" ,"aquamarine4","steelblue4")
colours_Endothelial<-c("lightpink3" , "slategray3","mediumpurple2","navy","#00B6EB","turquoise4")
colours_Myeloid_new<-c("slategray3","#00B6EB", "navy","darkgreen","darkmagenta","mediumorchid2","slateblue2","dodgerblue3")

#load in data
Myeloid_postclean <- readRDS("data/Myeloid_postclean.rds")
Myeloid_postclean$disease<-factor(x=Myeloid_postclean$disease, levels = c("lean", "obese", "MASL","MASH"))
Myeloid_postclean$new_subcluster<-factor(x=Myeloid_postclean$new_subcluster, levels = c("KC", "GPNMB Mac", "TransMac","preMac",  "Monocyte",  "cDC1", "cDC2", "migDC"))
Idents(Myeloid_postclean)<-Myeloid_postclean$new_subcluster

Myeloid_postclean$disease_MASH<-Myeloid_postclean$disease
Idents(Myeloid_postclean)<-Myeloid_postclean$disease_MASH
Myeloid_postclean<-RenameIdents(Myeloid_postclean,"MASL"="No MASH","lean"="No MASH","obese"="No MASH")
Myeloid_postclean$disease_MASH<-Idents(Myeloid_postclean)

##Cibersort
#Disable quantile normalization true
#Min expression 0.25
#Filter non-hematopoetic true

#Batch correction enabled, s-mode
#Disable quantily normalization true (microarray off)
#Run mode relative
#100 permutations
sub<-subset(Myeloid_postclean, subset= new_subcluster %in% c("KC", "GPNMB Mac",  "Monocyte"))
sub$new_subcluster<-factor(x=sub$new_subcluster, levels = c("KC", "GPNMB Mac",  "Monocyte"))
Idents(sub)<-sub$new_subcluster
cell_counts_df =data.frame(sub[["RNA"]]@counts)
celltype_labels = t(sub@meta.data$new_subcluster)
colnames(cell_counts_df) = celltype_labels
combined<-rbind(names(cell_counts_df),cell_counts_df)
rownames(combined)[rownames(combined) == "1"] <- "GeneSymbol"
df <- data.frame(GeneSymbol = row.names(combined), combined)

data.table::fwrite(df, file = "data/Myeloid_postclean_Mac_Mon_cibersort.txt", sep = "/t", col.names=F, row.names=F)


#bulk-seq
bulk <- read_excel("RNAseq_n216_normalised_batch_sex_IDcleaned.xlsx")
bulk2<-bulk[-1,]
library('biomaRt')
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
genes <- bulk2$`Sample ID`
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

merged<-merge(bulk2,G_list,by.x="Sample ID",by.y="ensembl_gene_id")
merged<-merged[,-2]
merged2<-merged[,c(1,218,2:217)]
names(merged2)[names(merged2) == 'hgnc_symbol'] <- 'Gene_Name'
merged2<-merged2[,-1]
data.table::fwrite(merged2, file = "data/O_bulk_cibersort.txt", sep = "/t", col.names=T, row.names=F)

#cibersort graphs
library("readxl")
ciber_results<-data.table::fread("data/CIBERSORTx_Job15_Results.csv",header = T, sep=",")

rownames(ciber_results)<-ciber_results$Mixture
ciber_results$Mixture<-NULL
ciber_results$`P-value`<-NULL
ciber_results$Correlation<-NULL
ciber_results$RMSE<-NULL
table<-as.data.frame.matrix(ciber_results)
table.perc<-apply(table[],2,function (x){x*100})
table.perc.round<-as.data.frame(round(table.perc, digits = 2)) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-as.character(bulk[1,][3:218])
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("normal", "NAFL",  "MASH_F0","MASH_F1","MASH_F2","MASH_F3", "cirrhosis"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")
#longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("LAM"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_disease<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

Myeloid_box_disease_stats<-Myeloid_box_disease+geom_pwc(
  aes(group = disease), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_disease, filename="Cibersort_deconvolute_Myeloid_MASH.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_disease_stats, filename="Cibersort_deconvolute_Myeloid_MASH_stats.pdf" ,height=5, width=5, units="in", dpi=320)


table.perc.round$disease_MASH<-table.perc.round$disease
table.perc.round$disease_MASH<-gsub("MASH_F0", "MASH_F0_1", table.perc.round$disease_MASH)
table.perc.round$disease_MASH<-gsub("MASH_F1", "MASH_F0_1", table.perc.round$disease_MASH)

table.perc.round$disease_MASH<-factor(table.perc.round$disease_MASH, levels = c("normal", "NAFL",  "MASH_F0_1","MASH_F2","MASH_F3", "cirrhosis"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")
#longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("LAM"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_disease_MASH<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease_MASH", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_disease_MASH_stats<-Myeloid_box_disease_MASH+geom_pwc(
  aes(group = disease_MASH), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_disease_MASH, filename="Cibersort_deconvolute_Myeloid_disease_MASH.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_disease_MASH_stats, filename="Cibersort_deconvolute_Myeloid_disease_MASH_stats.pdf" ,height=5, width=5, units="in", dpi=320)



table.perc.round$Fibrosis<-table.perc.round$disease
table.perc.round$Fibrosis<-gsub("normal", "F0", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("NAFL", "F0", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("MASH_F0", "F0", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("MASH_F1", "F1", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("MASH_F2", "F2-4", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("MASH_F3", "F2-4", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("cirrhosis", "F2-4", table.perc.round$Fibrosis)

table.perc.round$Fibrosis<-factor(table.perc.round$Fibrosis, levels = c("F0", "F1",  "F2-4"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")
#longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("LAM"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_fibrosis<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="Fibrosis", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_fibrosis_stats<-Myeloid_box_fibrosis+geom_pwc(
  aes(group = Fibrosis), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_fibrosis, filename="Cibersort_deconvolute_Myeloid_fibrosis.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_fibrosis_stats, filename="Cibersort_deconvolute_Myeloid_fibrosis_stats.pdf" ,height=5, width=5, units="in", dpi=320)

#cirrhosis Histology      NAFL   MASH_F0   MASH_F1   MASH_F2   MASH_F3    normal 
# 14         1            51         5        29        53        54        10


table.perc.round$MASH<-table.perc.round$disease
table.perc.round$MASH<-gsub("MASH_F0", "MASH", table.perc.round$MASH)
table.perc.round$MASH<-gsub("MASH_F1", "MASH", table.perc.round$MASH)
table.perc.round$MASH<-gsub("MASH_F2", "MASH", table.perc.round$MASH)
table.perc.round$MASH<-gsub("MASH_F3", "MASH", table.perc.round$MASH)
table.perc.round$MASH<-gsub("cirrhosis", "MASH", table.perc.round$MASH)

table.perc.round$MASH<-factor(table.perc.round$MASH, levels = c("normal", "NAFL",  "MASH"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")
#longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("LAM"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_MASH<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="MASH", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
#Myeloid_box_MASH_stats<-
Myeloid_box_MASH+geom_pwc(
  aes(group = MASH), tip.length = 0,
  method = "t_test", label = "p.adj.format", p.adjust.method = "holm",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_MASH, filename="Cibersort_deconvolute_Myeloid_MASH.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_MASH_stats, filename="Cibersort_deconvolute_Myeloid_MASH_stats.pdf" ,height=5, width=5, units="in", dpi=320)


bulk_saf<-read_excel("RNAseq_n216_normalised_batch_sex_IDcleaned.xlsx", sheet = "Sheet1")
table.perc.round<-merge(table.perc.round, bulk_saf, by.x="patient", by.y="Filename")
table.perc.round$NAS_score<-table.perc.round$`NAFLD Activity Score (NAS) Kleiner`
table.perc.round$NAS_score<-gsub("1|2|3", "1-3", table.perc.round$NAS_score)
table.perc.round$NAS_score<-gsub("4|5|6|7|8", "4-8", table.perc.round$NAS_score)

table.perc.round$NAS_score<-factor(table.perc.round$NAS_score, levels = c("0", "1-3",  "4-8"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")
#longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("LAM"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_NAS_score<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="NAS_score", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_NAS_score_stats<-Myeloid_box_NAS_score+geom_pwc(
  aes(group = NAS_score), tip.length = 0,
  method = "t_test", label = "p.adj.format", p.adjust.method = "holm",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_NAS_score, filename="Cibersort_deconvolute_Myeloid_NAS_score.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_NAS_score_stats, filename="Cibersort_deconvolute_Myeloid_NAS_score_stats.pdf" ,height=5, width=5, units="in", dpi=320)

table.perc.round$FLIP_score<-table.perc.round$`FLIP Activity score`

table.perc.round$FLIP_score<-factor(table.perc.round$FLIP_score, levels = c("0", "1", "2", "3","4"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")
#longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("LAM"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_FLIP_score<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="FLIP_score", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_FLIP_score_stats<-Myeloid_box_FLIP_score+geom_pwc(
  aes(group = FLIP_score), tip.length = 0,
  method = "t_test", label = "p.adj.format", p.adjust.method = "holm",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_FLIP_score, filename="Cibersort_deconvolute_Myeloid_FLIP_score.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_FLIP_score_stats, filename="Cibersort_deconvolute_Myeloid_FLIP_score_stats.pdf" ,height=5, width=5, units="in", dpi=320)

