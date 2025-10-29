#Figure 2 plots
#save all necessary files if needed with ggsave
set.seed(1234567)
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
library(SpatialDecon)

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
Liver_clean <- readRDS("data/Liver_clean.rds")
Liver_clean$disease<-factor(x=Liver_clean$disease, levels = c("lean", "obese", "MASL","MASH"))
Liver_clean$Cluster<-factor(x=Liver_clean$Cluster, levels = c("Hepatocytes",  "Cholangiocytes",    "Endothelial", "Mesenchymal", "Lymphocytes","Myeloid-like cells", "B cells", "pDC"))
DefaultAssay(Liver_clean)<-"RNA"

#GeoMX deconvolution whole Liver
geomx_oliver_annot<-read_xlsx("030_sDAS_Govaere_Count_Matrices.adj.xlsx", "Annotations", col_names = F)
geomx_oliver_annot <- geomx_oliver_annot[!is.na(geomx_oliver_annot$...1),]
geomx_oliver_annot<-as.data.frame(geomx_oliver_annot)
rownames(geomx_oliver_annot)<-geomx_oliver_annot$...1
geomx_oliver_annot<-geomx_oliver_annot[,-c(1,3,8,10,11,12,13,14,15,16)]
colnames(geomx_oliver_annot)<-c("Cell", "Marker", "Hep_Area","Area_short","ROI", "ROI_number")
geomx_oliver_annot<-na.omit(geomx_oliver_annot)

geomx_oliver_norm<-read_xlsx("030_sDAS_Govaere_Count_Matrices.adj.xlsx", "Q3_Normalized_Counts")
geomx_oliver_norm<-geomx_oliver_norm[-c(1,2),]
rownames(geomx_oliver_norm)<-geomx_oliver_norm$...1
geomx_oliver_norm<-as.matrix(geomx_oliver_norm)
geomx_oliver_norm<-geomx_oliver_norm[,-c(1)]
geomx_oliver_norm<-geomx_oliver_norm[,colnames(geomx_oliver_norm) %in% rownames(geomx_oliver_annot)]
class(geomx_oliver_norm) <- "numeric"

geomx_oliver_raw<-read_xlsx("030_sDAS_Govaere_Count_Matrices.adj.xlsx", "Raw_probe_data")
duplicates <- which(duplicated(geomx_oliver_raw$TargetName) | duplicated(geomx_oliver_raw$TargetName, fromLast = TRUE))
duplicated_negs<-geomx_oliver_raw$TargetName[duplicates]
for (i in duplicates) {
  positions <- which(geomx_oliver_raw$TargetName == geomx_oliver_raw$TargetName[i])
  if (length(positions) > 1) {
    for (j in seq_along(positions)) {
      geomx_oliver_raw$TargetName[positions[j]] <- paste(geomx_oliver_raw$TargetName[positions[j]], j, sep = "_")
    }
  }
}
rownames(geomx_oliver_raw)<-geomx_oliver_raw$TargetName
geomx_oliver_raw<-as.matrix(geomx_oliver_raw)
geomx_oliver_raw<-geomx_oliver_raw[,-c(1,2,3)]
geomx_oliver_raw<-geomx_oliver_raw[,colnames(geomx_oliver_raw) %in% colnames(geomx_oliver_norm)]
class(geomx_oliver_raw) <- "numeric"

colnames(geomx_oliver_norm) = colnames(geomx_oliver_raw) = rownames(geomx_oliver_annot) = 
  paste0(geomx_oliver_annot$ROI,geomx_oliver_annot$Cell)

bg2 = derive_GeoMx_background(norm = geomx_oliver_norm,
                              probepool = rep(1, nrow(geomx_oliver_norm)),
                              negnames = duplicated_negs)

Idents(Liver_clean)<-Liver_clean$Cluster
sub<-Liver_clean
annots <- data.frame(cbind(cellType=as.character(Idents(sub)), 
                           cellID=names(Idents(sub))))

custom_mtx_seurat <- create_profile_matrix(mtx = Seurat::GetAssayData(object = sub, 
                                                                      assay = "RNA", 
                                                                      slot = "counts"), 
                                           cellAnnots = annots, 
                                           cellTypeCol = "cellType", 
                                           cellNameCol = "cellID", 
                                           matrixName = "Liver_MASH",
                                           outDir = NULL, 
                                           normalize = T, 
                                           minCellNum = 5, 
                                           minGenes = 10)
head(custom_mtx_seurat)

res = spatialdecon(norm = geomx_oliver_norm,
                   bg = bg2,
                   X = custom_mtx_seurat,
                   align_genes = TRUE)

ordering<-c("Hepatocytes",  "Cholangiocytes",    "Endothelial", "Mesenchymal", "Lymphocytes","Myeloid-like cells", "B cells", "pDC")
res$beta<-res$beta[rev(ordering),]
colnames(res$beta)<- gsub("cd68", "CD68", colnames(res$beta))
colnames(res$beta)<- gsub("Cd68", "CD68", colnames(res$beta))
colnames(res$beta)<- gsub("panCK", "PanCK", colnames(res$beta))
colnames(res$beta)<- gsub("PanCk", "PanCK", colnames(res$beta))
colnames(res$beta)<- gsub("panCk", "PanCK", colnames(res$beta))
colnames(res$beta)<- gsub("cd45", "CD45", colnames(res$beta))
colnames(res$beta)<- gsub("Cd45", "CD45", colnames(res$beta))
col_names<-colnames(res$beta)

df<-data.frame("position"=1:77, "Names"=substr(col_names, nchar(col_names) - 1, nchar(col_names)),
               "Marker"=unlist(regmatches(colnames(res$beta), gregexpr(("CD45|CD68|PanCK"), colnames(res$beta)))))
# Select sort order for X1: y, x
Marker.values = c("PanCK","CD45", "CD68")
# Select sort order for X2: c, a, b
Names.values = c("PT","HN","HL")

# Place call to factor for both columns
df = df[order(factor(df$Marker, levels = Marker.values), factor(df$Names, levels = Names.values)),]

res$beta <- res$beta[, df$position]


cols_Marker<-c(rep("lightpink3",27), rep("00B6EB",17), rep("turquoise4",33))
cols_Names<-c(rep("sienna3",9),rep("red4",18),rep("steelblue3",7),rep("sienna3",1),rep("red4",9),rep("steelblue3",7),rep("sienna3",8),rep("red4",18))
myCols = cbind(cols_Marker, cols_Names)
rownames(myCols)<-colnames(res$beta)

pdf("GeoMX_deconv_all.pdf", height = 10, width = 10)
heatmap(res$beta, cexCol = 0.5, cexRow = 1, margins = c(10,10),Colv = NA, Rowv = NA)
dev.off()

res$prop_of_all<-res$prop_of_all[rev( ordering),]
colnames(res$prop_of_all)<- gsub("cd68", "CD68", colnames(res$prop_of_all))
colnames(res$prop_of_all)<- gsub("Cd68", "CD68", colnames(res$prop_of_all))
colnames(res$prop_of_all)<- gsub("panCK", "PanCK", colnames(res$prop_of_all))
colnames(res$prop_of_all)<- gsub("PanCk", "PanCK", colnames(res$prop_of_all))
colnames(res$prop_of_all)<- gsub("panCk", "PanCK", colnames(res$prop_of_all))
colnames(res$prop_of_all)<- gsub("cd45", "CD45", colnames(res$prop_of_all))
colnames(res$prop_of_all)<- gsub("Cd45", "CD45", colnames(res$prop_of_all))
res$prop_of_all <- res$prop_of_all[, df$position]

pdf("GeoMX_deconv_all_barplot_gapped.pdf", height = 10, width = 10)
TIL_barplot(res$prop_of_all,draw_legend = F, cex.names = 0.3, col = rev(c(colours_Liver)),space=c(rep(0.05,8),0.2,rep(0.05,17), 0.7, rep(0.05,6),0.2,0.2,rep(0.05,8), 0.7,rep(0.05,6),0.2,rep(0.05,7),0.2,rep(0.05,18)))
dev.off()


