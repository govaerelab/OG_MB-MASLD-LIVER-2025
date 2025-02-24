#Figure 3 plots
#save all necessary files if needed with ggsave
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
Myeloid_postclean <- readRDS("data/Myeloid_postclean.rds")
Myeloid_postclean$disease<-factor(x=Myeloid_postclean$disease, levels = c("lean", "obese", "MASL","MASH"))
Myeloid_postclean$new_subcluster<-factor(x=Myeloid_postclean$new_subcluster, levels = c("KC", "GPNMB Mac", "TransMac","preMac",  "Monocyte",  "cDC1", "cDC2", "migDC"))
Idents(Myeloid_postclean)<-Myeloid_postclean$new_subcluster

Myeloid_postclean$disease_MASH<-Myeloid_postclean$disease
Idents(Myeloid_postclean)<-Myeloid_postclean$disease_MASH
Myeloid_postclean<-RenameIdents(Myeloid_postclean,"MASL"="No MASH","lean"="No MASH","obese"="No MASH")
Myeloid_postclean$disease_MASH<-Idents(Myeloid_postclean)

#GeoMX deconvolution Myeloids
geomx_oliver_annot<-read_xlsx("030_sDAS_Govaere_Count_Matrices.adj.xlsx", "Annotations", col_names = F)
geomx_oliver_annot <- geomx_oliver_annot[!is.na(geomx_oliver_annot$...1),]
geomx_oliver_annot<-as.data.frame(geomx_oliver_annot)
rownames(geomx_oliver_annot)<-geomx_oliver_annot$...1
geomx_oliver_annot<-geomx_oliver_annot[,-c(1,3,8,10,11,12,13,14,15,16)]
colnames(geomx_oliver_annot)<-c("Cell", "Marker", "Hep_Area","Area_short","ROI", "ROI_number")
geomx_oliver_annot<-geomx_oliver_annot[geomx_oliver_annot$Marker%in% c("Cd68","cd68", "CD68"),]
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

Idents(Myeloid_downstream)<-Myeloid_downstream$new_subcluster
sub<-subset(Myeloid_downstream, subset=new_subcluster %in% c("cDC2","migDC","cDC1", "preMac","TransMac"), invert=T)
annots <- data.frame(cbind(cellType=as.character(Idents(sub)), 
                           cellID=names(Idents(sub))))

custom_mtx_seurat <- create_profile_matrix(mtx = Seurat::GetAssayData(object = sub, 
                                                                      assay = "RNA", 
                                                                      slot = "counts"), 
                                           cellAnnots = annots, 
                                           cellTypeCol = "cellType", 
                                           cellNameCol = "cellID", 
                                           matrixName = "Myeloids_MASH",
                                           outDir = NULL, 
                                           normalize = T, 
                                           minCellNum = 5, 
                                           minGenes = 10)
head(custom_mtx_seurat)

res = spatialdecon(norm = geomx_oliver_norm,
                   bg = bg2,
                   X = custom_mtx_seurat,
                   align_genes = TRUE)


res$beta<-res$beta[c(  "GPNMB Mac","KC","Monocyte"),]

col_names<-colnames(res$beta)
df<-data.frame("position"=1:33, "Names"=substr(col_names, nchar(col_names) - 1, nchar(col_names)))
PTs<-df[df$Names %in% "PT",]$position
HLs<-df[df$Names %in% "HL",]$position
HNs<-df[df$Names %in% "HN",]$position


# Reorder columns based on sorted indices
res$beta <- res$beta[, c(PTs,HNs,HLs)]

colnames(res$beta)<- gsub("cd68", "CD68", colnames(res$beta))
colnames(res$beta)<- gsub("Cd68", "CD68", colnames(res$beta))
colnames(res$prop_of_all)<- gsub("cd68", "CD68", colnames(res$prop_of_all))
colnames(res$prop_of_all)<- gsub("Cd68", "CD68", colnames(res$prop_of_all))

pdf("GeoMX_deconv_myeloid_coloured.pdf", height = 10, width = 10)
heatmap(res$beta, cexCol = 0.5, cexRow = 1, margins = c(10,10),Colv = NA, Rowv = NA,ColSideColors = c(rep("steelblue3",7),rep("sienna3",8),rep("red4",18)))
dev.off()

pdf("GeoMX_deconv_myeloid_2.pdf", height = 13, width = 10)
heatmap(res$beta, cexCol = 1, cexRow = 1, margins = c(10,7),Colv = NA, Rowv = NA,col=colorRampPalette(c("white", "#0092b5", "#0092b5", "#a6ce39", "#a6ce39"))(101))
dev.off()



res$prop_of_all<-res$prop_of_all[c(  "GPNMB Mac","KC","Monocyte"),]

col_names<-colnames(res$prop_of_all)
df<-data.frame("position"=1:33, "Names"=substr(col_names, nchar(col_names) - 1, nchar(col_names)))
PTs<-df[df$Names %in% "PT",]$position
HLs<-df[df$Names %in% "HL",]$position
HNs<-df[df$Names %in% "HN",]$position
# Reorder columns based on sorted indices
res$prop_of_all <- res$prop_of_all[, c(PTs,HNs,HLs)]

pdf("GeoMX_deconv_myeloid_barplot_gaps.pdf", height = 10, width = 10)
TIL_barplot(res$prop_of_all,draw_legend = F, cex.names = 0.3,col = rev(c("darkmagenta", "slategray3","#00B6EB")),space=c(rep(0.05,7),0.5,rep(0.05,7), 0.5, rep(0.05,17)))
dev.off()

#DGE Volcano GPNMB Macs
Idents(Myeloid_postclean)<-Myeloid_postclean$disease_MASH
sub<-subset(Myeloid_postclean, subset=new_subcluster %in% "GPNMB Mac")
DefaultAssay(sub)<-"RNA"
sub<-NormalizeData(object = sub, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
subMM_volc <- FindMarkers(object = sub,ident.1 = "No MASH", ident.2 = "MASH",  min.pct = 0.10, logfc.threshold = 0,  slot="data", test.use = "LR", latent.vars = c("nCount_RNA","percent.mito"), verbose=F)

subMM_volc_sig = subMM_volc[which(subMM_volc$p_val_adj<0.05),]
subMM_volc_sig<-subMM_volc_sig %>% dplyr::rename(No_MASH=pct.1, MASH=pct.2)
keyvals <- rep("grey", nrow(subMM_volc))
# set the base name/label as "NS"
names(keyvals) <- rep("NS", nrow(subMM_volc))
# modify keyvals for variables with fold change >= 1
keyvals[which(subMM_volc$avg_log2FC >= 0.5 & subMM_volc$p_val_adj <= 0.05)] <- "sienna3"
names(keyvals)[which(subMM_volc$avg_log2FC >= 0.5 & subMM_volc$p_val_adj <= 0.05)] <- "no MASH"
# modify keyvals for variables with fold change <= -1
keyvals[which(subMM_volc$avg_log2FC <= -0.5 & subMM_volc$p_val_adj <= 0.05)] <- "red4"
names(keyvals)[which(subMM_volc$avg_log2FC <= -0.5 & subMM_volc$p_val_adj <= 0.05)] <- "MASH"
library(EnhancedVolcano)
genes<-c("LPL", "GPNMB", "FABP5", "CD9", "LYZ", "HLA-DRA", "ACP5", "SPP1")
EnhancedVolcano(subMM_volc,
                lab = rownames(subMM_volc),
                x = 'avg_log2FC',
                y = 'p_val_adj', parseLabels = T,drawConnectors = T, widthConnectors = 0.5,arrowheads=F,labFace = 'bold', boxedLabels = T,
                title = "", subtitle=NULL, gridlines.major = F, gridlines.minor = F, selectLab = genes,
                colCustom = keyvals, FCcutoff = 0.5, pCutoff = 0.05,ylim=c(0,120), xlim =c(-4,3), caption = "")+ 
  FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))

#CosMX, switch to Seurat_5 for spatial analysis
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize= 8912896000)

Leuven1<-LoadNanostring("CosMX/Leuven_1/", fov = "Leuven1")
Leuven1$orig.ident<-"Leuven1"
Leuven2<-LoadNanostring("CosMX/Leuven_2/", fov = "Leuven2")
Leuven2$orig.ident<-"Leuven2"
Leuven3<-LoadNanostring("CosMX/Leuven_3/", fov = "Leuven3")
Leuven3$orig.ident<-"Leuven3"
Leuven4<-LoadNanostring("CosMX/Leuven_4/", fov = "Leuven4")
Leuven4$orig.ident<-"Leuven4"

names<-c("Leuven1","Leuven2","Leuven3","Leuven4")
merged_CosMX <- merge(Leuven1, y = list(Leuven2,Leuven3,Leuven4), add.cell.ids =names, project = "CosMX_Liver")
DefaultAssay(merged_CosMX)<-"Nanostring"
merged_CosMX <-SCTransform(merged_CosMX, assay = "Nanostring",clip.range = c(-10, 10))
#slot(object = Xenium.big@assays$SCT@SCTModel.list[[2]], name="umi.assay")
merged_CosMX <- RunPCA(merged_CosMX, npcs = 30)
merged_CosMX <- RunUMAP(merged_CosMX, dims = 1:30, min.dist = 0.1)
merged_CosMX <- FindNeighbors(merged_CosMX, reduction = "pca", dims=1:30, k.param = 30)
merged_CosMX <- FindClusters(merged_CosMX, resolution = 0.3)

Idents(merged_CosMX)<-merged_CosMX$SCT_snn_res.0.3
new.cluster.ids <- c("Hepatocytes", "Mesenchymal", "Myeloid-like cells", "Hepatocytes",
                     "Cholangiocytes", "Endothelial", "Hepatocytes", "Lympho-/Granulocytes", "B cells", 
                     "Mesenchymal", "Mesenchymal", "B cells", "Mesenchymal", "Hepatocytes",
                     "B cells", "Hepatocytes", "B cells", "Lympho-/Granulocytes", "Hepatocytes", 
                     "Hepatocytes", "B cells", "Lympho-/Granulocytes", "Mesenchymal"
)
names(new.cluster.ids) <- seq(0,22)
merged_CosMX <- RenameIdents(merged_CosMX, new.cluster.ids)
merged_CosMX$Cluster<-Idents(merged_CosMX)

colours_CosMX<-c("lightpink3" , "slategray3","mediumpurple2","navy","#00B6EB","turquoise4", "burlywood")
merged_CosMX$Cluster<-factor(x=merged_CosMX$Cluster, levels = c("Hepatocytes",  "Cholangiocytes",    "Endothelial", "Mesenchymal", "Lympho-/Granulocytes","Myeloid-like cells", "B cells"))
Idents(merged_CosMX)<-merged_CosMX$Cluster
DefaultAssay(merged_CosMX)<-"Nanostring"

names<-c("Leuven1", "Leuven2","Leuven3","Leuven4")
for (i in names){
  DefaultBoundary(merged_CosMX[[i]]) <- "centroids"
}

cropped.coords <- Crop(merged_CosMX[["Leuven1"]], y = c(40000, 125000), x = c(10000, 72000), coords = "plot")
merged_CosMX[["Leuven1lower"]] <- cropped.coords
cropped.coords <- Crop(merged_CosMX[["Leuven1"]], y = c(40000, 75000), x = c(50000, 72000), coords = "plot")
merged_CosMX[["Leuven1zoom1"]] <- cropped.coords
cropped.coords <- Crop(merged_CosMX[["Leuven1"]], y = c(60000, 75000), x = c(57000, 72000), coords = "plot")
merged_CosMX[["Leuven1zoom2"]] <- cropped.coords
cropped.coords <- Crop(merged_CosMX[["Leuven1"]], y = c(61000, 70000), x = c(65000, 72000), coords = "plot")
merged_CosMX[["Leuven1zoom3"]] <- cropped.coords

#plot1<-ImageDimPlot(merged_CosMX,fov="Leuven1lower",group.by = "Cluster",  alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F, axes = F)+theme(panel.grid = element_blank())
plot2<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom1",group.by = "Cluster", size=1, alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F, axes = F)+theme(panel.grid = element_blank())
plot3<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom2",group.by = "Cluster",size=1.5,  alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F, axes = F)+theme(panel.grid = element_blank())
plot4<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom3",group.by = "Cluster",size=2,  alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F, axes = F)+theme(panel.grid = element_blank())


genes_to_plot<-c("IL32","GPNMB")
colours_for_genes<-c("red", "green")
#plot1<-ImageDimPlot(merged_CosMX,fov="Leuven1lower",group.by = "Cluster", molecules=genes_to_plot, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot2<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom1",group.by = "Cluster", molecules=genes_to_plot, size=1, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot3<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom2",group.by = "Cluster", molecules=genes_to_plot,size=1.5, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot4<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom3",group.by = "Cluster", molecules=genes_to_plot,size=2, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.2, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())

genes_to_plot<-c("HLA-DRA","GPNMB", "FABP5","LPL")
colours_for_genes<-c("red", "green", "blue", "orange")
#plot1<-ImageDimPlot(merged_CosMX,fov="Leuven1lower",group.by = "Cluster", molecules=genes_to_plot, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot2<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom1",group.by = "Cluster", molecules=genes_to_plot,size=1, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot3<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom2",group.by = "Cluster", molecules=genes_to_plot,size=1.5, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot4<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom3",group.by = "Cluster", molecules=genes_to_plot,size=2, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.1, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())


genes<-c("IL32", "GPNMB","FABP5","LPL","MSR1","CD163","HLA-DRA") 
for(i in 1:length(genes)){
  plot1<- ImageFeaturePlot(merged_CosMX, fov="Leuven1zoom3",genes[i], combine=T,  coord.fixed = T, size=2,crop=T, axes=F,border.size = NA, min.cutoff=0, max.cutoff=5)+coord_flip()+theme(panel.grid = element_blank())
  ggsave(plot=plot1, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_zoom3_CosMX_featureplot_",genes[i],".pdf", sep = "") ,height=4.5, width=9, units="in", dpi=320)
  
}







#Supplemental
dotplot_Liver_genes<-c("LPL", "GPNMB", "FABP5",  "LYZ",  "ACP5", "SPP1")
dotplot_Liver<-DotPlot(Liver_clean, group.by= "Cluster",features = rev(dotplot_Liver_genes) , assay = "RNA", cols = c("steelblue2", "red"))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) +xlab("")+ylab("")+coord_flip()

DefaultAssay(Myeloid_postclean)<-"RNA"
genes<-c("LPL", "GPNMB", "FABP5","TREM2", "CD9", "LYZ", "HLA-DRA", "CD163", "SPP1","MSR1", "ACP5")
for(i in genes){
  feature_cite<-FeaturePlot(Myeloid_postclean, label = F, features = i, order = T, pt.size = 1.2, max.cutoff =4, min.cutoff = 2 , cols = c("grey", "brown4"))+theme( 
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(),axis.line = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), axis.title = element_blank())+NoLegend()+ggtitle("")
  ggsave(plot=feature_cite, filename=paste(getwd(),"/Figures/Olivier/",i,"_Feature_Myeloid.pdf", sep = "") ,height=3, width=3, units="in", dpi=320)
}

genes<-c("LPL", "GPNMB", "FABP5", "TREM2", "CD9", "LYZ", "HLA-DRA", "CD163", "SPP1","MSR1", "ACP5")
for(i in genes){
  vlnplot<-VlnPlot(subset(Myeloid_downstream, subset=new_subcluster %in% c("KC", "GPNMB Mac","TransMac","preMac",  "Monocyte" )), assay = "RNA", pt.size = 0, features = i, group.by = "new_subcluster", split.by = "disease", cols = colours_disease)+NoLegend()+theme( axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize), axis.text.y = element_text(size=fontsize))+NoLegend()+xlab("")
  ggsave(plot=vlnplot, filename=paste(getwd(),"/Figures/Olivier/",i,"_VlnPlot_Myeloid_disease.pdf", sep = "") ,height=3, width=3, units="in", dpi=320)
}


#Get proportion of cells expressing certain markers per disease

#combination and percentage of >0 counts
get_percentage<-function(disease_group, gene_vector){
  sub<-subset(Myeloid_postclean, subset=disease_MASH %in% disease_group)
  df<-sub@assays$RNA@counts
  # Filter the dataframe for only the genes in the vector
  filtered_df <- df[rownames(df) %in% gene_vector, ]
  
  # Function to calculate proportions of cells with counts > 0 for a gene combination
  calculate_proportion <- function(gene_combination, data) {
    # Subset the data based on the gene combination
    subset_data <- data[gene_combination, , drop = FALSE]  # drop = FALSE ensures it's still a dataframe
    
    # For single gene combinations, we check directly for > 0 counts
    if (length(gene_combination) == 1) {
      # Calculate the proportion of cells with counts > 0
      count_above_zero <- sum(subset_data > 0)
      total_cells <- ncol(subset_data)
      proportion <- ifelse(total_cells > 0, count_above_zero / total_cells, 0)
      
    } else {
      # For multiple genes, check if all conditions are > 0
      condition_met <- apply(subset_data, 2, function(col) all(col > 0))
      total_cells <- ncol(subset_data)
      proportion <- ifelse(total_cells > 0, sum(condition_met) / total_cells, 0)
    }
    
    return(proportion)  # Return the calculated proportion
  }
  
  # Initialize an empty data frame to store the results
  results <- data.frame(genes = character(), proportion = numeric(), stringsAsFactors = FALSE)
  
  # Loop through combinations of 1, 2, and 3 genes
  for (i in 1:3) {
    gene_combinations <- combn(gene_vector, i, simplify = FALSE)
    
    for (combination in gene_combinations) {
      proportion <- calculate_proportion(combination, filtered_df)
      gene_comb_name <- paste(combination, collapse = " & ")
      
      # Append the result to the results data frame
      results <- rbind(results, data.frame(genes = gene_comb_name, proportion = proportion, stringsAsFactors = FALSE))
    }
  }
  
  # Print the final results
  results$proportion<-   results$proportion*100 
  colnames(results)<-  c("genes",paste("proportion_",disease_group,"(%)",sep=""))                    
  return(results)}

Myeloid_postclean$disease_MASH<-Myeloid_postclean$disease
Idents(Myeloid_postclean)<-Myeloid_postclean$disease_MASH
Myeloid_postclean<-RenameIdents(Myeloid_postclean,"MASL"="No MASH","lean"="No MASH","obese"="No MASH")
Myeloid_postclean$disease_MASH<-Idents(Myeloid_postclean)

genes<-c("TREM2", 'CD9', 'TIMD4', 'S100A4', 'LPL', 'GPNMB', 'SPP1', 'CD163', 'MARCO', 'CD5L', 'CCL19')
MASH_proportion<-get_percentage(gene_vector=genes, disease_group="MASH")                          
no_MASH_proportion<-get_percentage(gene_vector=genes, disease_group="No MASH") 

merged<-merge(MASH_proportion, no_MASH_proportion) 
