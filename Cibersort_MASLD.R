##Cibersortx
#prepare snRNA-seq data
#GPNMB Mac was later changed to MetMac

sub<-subset(Myeloid_downstream, subset= new_subcluster %in% c("KC", "GPNMB Mac",  "Monocyte"))
sub$new_subcluster<-factor(x=sub$new_subcluster, levels = c("KC", "GPNMB Mac",  "Monocyte"))
Idents(sub)<-sub$new_subcluster
cell_counts_df =data.frame(sub[["RNA"]]@counts)
celltype_labels = t(sub@meta.data$new_subcluster)
colnames(cell_counts_df) = celltype_labels
combined<-rbind(names(cell_counts_df),cell_counts_df)
rownames(combined)[rownames(combined) == "1"] <- "GeneSymbol"
df <- data.frame(GeneSymbol = row.names(combined), combined)

data.table::fwrite(df, file = "data/Myeloid_downstream_Mac_Mon_cibersort.txt", sep = "/t", col.names=F, row.names=F)


#prepare bulk-seq data
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

#upload files and rub cibersortX
#Disable quantile normalization true
#Min expression 0.25
#Filter non-hematopoetic true

#Batch correction enabled, s-mode
#Disable quantily normalization true (microarray off)
#Run mode relative
#100 permutations


#Load-in CibersortX output
#Group together samples by different Conditions
#Generate boxplot, run t-test and adjust with BH
#GPNMB Mac was later changed to MetMac

library("readxl")
library(ggpubr)
library(ggplot2)
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
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("normal", "NAFL",  "NASH_F0","NASH_F1","NASH_F2","NASH_F3", "cirrhosis"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_disease<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

#run stats
Myeloid_box_disease_stats<-Myeloid_box_disease+geom_pwc(
  aes(group = disease), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_disease, filename="Cibersort_deconvolute_Myeloid_NASH.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_disease_stats, filename="Cibersort_deconvolute_Myeloid_NASH_stats.pdf" ,height=5, width=5, units="in", dpi=320)

#Combine F0 and F1 NASH samples
table.perc.round$disease_NASH<-table.perc.round$disease
table.perc.round$disease_NASH<-gsub("NASH_F0", "NASH_F0_1", table.perc.round$disease_NASH)
table.perc.round$disease_NASH<-gsub("NASH_F1", "NASH_F0_1", table.perc.round$disease_NASH)

table.perc.round$disease_NASH<-factor(table.perc.round$disease_NASH, levels = c("normal", "NAFL",  "NASH_F0_1","NASH_F2","NASH_F3", "cirrhosis"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_disease_NASH<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease_NASH", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_disease_NASH_stats<-Myeloid_box_disease_NASH+geom_pwc(
  aes(group = disease_NASH), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_disease_NASH, filename="Cibersort_deconvolute_Myeloid_disease_NASH.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_disease_NASH_stats, filename="Cibersort_deconvolute_Myeloid_disease_NASH_stats.pdf" ,height=5, width=5, units="in", dpi=320)


#Seperate based on F-score
table.perc.round$Fibrosis<-table.perc.round$disease
table.perc.round$Fibrosis<-gsub("normal", "F0", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("NAFL", "F0", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("NASH_F0", "F0", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("NASH_F1", "F1", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("NASH_F2", "F2-4", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("NASH_F3", "F2-4", table.perc.round$Fibrosis)
table.perc.round$Fibrosis<-gsub("cirrhosis", "F2-4", table.perc.round$Fibrosis)

table.perc.round$Fibrosis<-factor(table.perc.round$Fibrosis, levels = c("F0", "F1",  "F2-4"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_fibrosis<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="Fibrosis", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_fibrosis_stats<-Myeloid_box_fibrosis+geom_pwc(
  aes(group = Fibrosis), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_fibrosis, filename="Cibersort_deconvolute_Myeloid_fibrosis.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_fibrosis_stats, filename="Cibersort_deconvolute_Myeloid_fibrosis_stats.pdf" ,height=5, width=5, units="in", dpi=320)


#Combine based on disease, despite F-score
table.perc.round$NASH<-table.perc.round$disease
table.perc.round$NASH<-gsub("NASH_F0", "NASH", table.perc.round$NASH)
table.perc.round$NASH<-gsub("NASH_F1", "NASH", table.perc.round$NASH)
table.perc.round$NASH<-gsub("NASH_F2", "NASH", table.perc.round$NASH)
table.perc.round$NASH<-gsub("NASH_F3", "NASH", table.perc.round$NASH)
table.perc.round$NASH<-gsub("cirrhosis", "NASH", table.perc.round$NASH)

table.perc.round$NASH<-factor(table.perc.round$NASH, levels = c("normal", "NAFL",  "NASH"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_NASH<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="NASH", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

Myeloid_box_NASH_stats<-Myeloid_box_NASH+geom_pwc(
  aes(group = NASH), tip.length = 0,
  method = "t_test", label = "p.adj.format", p.adjust.method = "holm",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_NASH, filename="Cibersort_deconvolute_Myeloid_NASH.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_NASH_stats, filename="Cibersort_deconvolute_Myeloid_NASH_stats.pdf" ,height=5, width=5, units="in", dpi=320)

#Combine based on Activity/NAS score
bulk_saf<-read_excel("RNAseq_n216_normalised_batch_sex_IDcleaned.xlsx", sheet = "Sheet1")
table.perc.round<-merge(table.perc.round, bulk_saf, by.x="patient", by.y="Filename")
table.perc.round$NAS_score<-table.perc.round$`NAFLD Activity Score (NAS) Kleiner`
table.perc.round$NAS_score<-gsub("1|2|3", "1-3", table.perc.round$NAS_score)
table.perc.round$NAS_score<-gsub("4|5|6|7|8", "4-8", table.perc.round$NAS_score)

table.perc.round$NAS_score<-factor(table.perc.round$NAS_score, levels = c("0", "1-3",  "4-8"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_NAS_score<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="NAS_score", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_NAS_score_stats<-Myeloid_box_NAS_score+geom_pwc(
  aes(group = NAS_score), tip.length = 0,
  method = "t_test", label = "p.adj.format", p.adjust.method = "holm",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_NAS_score, filename="Cibersort_deconvolute_Myeloid_NAS_score.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_NAS_score_stats, filename="Cibersort_deconvolute_Myeloid_NAS_score_stats.pdf" ,height=5, width=5, units="in", dpi=320)

#Combine based on FLIP Activity score
table.perc.round$FLIP_score<-table.perc.round$`FLIP Activity score`

table.perc.round$FLIP_score<-factor(table.perc.round$FLIP_score, levels = c("0", "1", "2", "3","4"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac",  "Monocyte"),  names_to = "Cluster", values_to = "Percentage")

Myeloid_box_FLIP_score<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="FLIP_score", ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
Myeloid_box_FLIP_score_stats<-Myeloid_box_FLIP_score+geom_pwc(
  aes(group = FLIP_score), tip.length = 0,
  method = "t_test", label = "p.adj.format", p.adjust.method = "holm",
  bracket.nudge.y = -0.08,  hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=Myeloid_box_FLIP_score, filename="Cibersort_deconvolute_Myeloid_FLIP_score.pdf" ,height=5, width=5, units="in", dpi=320)
ggsave(plot=Myeloid_box_FLIP_score_stats, filename="Cibersort_deconvolute_Myeloid_FLIP_score_stats.pdf" ,height=5, width=5, units="in", dpi=320)
