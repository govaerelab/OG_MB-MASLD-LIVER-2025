#Analysis
set.seed(1234567)
setwd("") ##change to output directory
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
library(DropletUtils)
library(SoupX)

#run Soupx for ambient correction and write new matrix
ambient_removal<-function(folder, sample_name){
  filt.matrix <- Read10X (data.dir = paste0("filtered_feature_bc_matrix/", folder, sep=""))
  raw.matrix  <- Read10X (data.dir = paste0("raw_feature_bc_matrix/", folder, sep=""))
  srat    <-CreateSeuratObject(counts = filt.matrix , min.cells = 0, min.features  = 0, project = "liver", assay = "RNA")
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  srat    <- NormalizeData(object = srat, normalization.method = "LogNormalize", scale.factor = 10000)
  srat    <- FindVariableFeatures(object = srat, verbose=F)
  srat    <- ScaleData(object = srat)
  srat    <- RunPCA(object = srat, npcs = 30)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel  <- setDR(soup.channel, umap)
  soup.channel  <- autoEstCont(soup.channel)
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
  DropletUtils:::write10xCounts(paste0("SoupX/", sample_name, "_soupx_counts", sep=""), adj.matrix)
}

dir.create("./SoupX")

ambient_removal(folder = "Healthy controls/Pt 12", sample_name = "Healthy_Pt12" )
ambient_removal(folder = "Healthy controls/Pt 14", sample_name = "Healthy_Pt14" )
ambient_removal(folder = "Healthy controls/Pt 16", sample_name = "Healthy_Pt16" )
ambient_removal(folder = "Healthy controls/Pt 45", sample_name = "Healthy_Pt45" )

ambient_removal (folder = "obese controls/LB222",sample_name = "obese_LB222" )
ambient_removal (folder = "obese controls/LB226",sample_name = "obese_LB226")
ambient_removal (folder = "obese controls/LB228",sample_name = "obese_LB228")
ambient_removal (folder = "obese controls/LB231",sample_name = "obese_LB231")
ambient_removal (folder = "obese controls/LB235",sample_name = "obese_LB235")
ambient_removal (folder = "NAFL/21B6674",sample_name = "NAFL_21B6674")
ambient_removal (folder = "NAFL/AB54",sample_name = "NAFL_AB54")
ambient_removal (folder = "NAFL/LB219",sample_name = "NAFL_LB219")
ambient_removal (folder = "NAFL/LB227",sample_name = "NAFL_LB227")
ambient_removal (folder = "NASH/20B6320",sample_name = "NASH_20B6320")
ambient_removal (folder = "NASH/21B948",sample_name = "NASH_21B948")
ambient_removal (folder = "NASH/AB61",sample_name = "NASH_AB61")
ambient_removal (folder = "NASH/LB224",sample_name = "NASH_LB224")
ambient_removal (folder = "NASH/LB234",sample_name = "NASH_LB234")

#load in samples and add metadata
p12.data <- Read10X (data.dir = "SoupX/Healthy_Pt12_soupx_counts")
p14.data <- Read10X (data.dir = "SoupX/Healthy_Pt14_soupx_counts")
p16.data <- Read10X (data.dir = "SoupX/Healthy_Pt16_soupx_counts")
p45.data<- Read10X (data.dir = "SoupX/Healthy_Pt45_soupx_counts")

LB222.data <- Read10X (data.dir = "SoupX/obese_LB222_soupx_counts")
LB226.data <- Read10X (data.dir = "SoupX/obese_LB226_soupx_counts")
LB228.data <- Read10X (data.dir = "SoupX/obese_LB228_soupx_counts")
LB231.data <- Read10X (data.dir = "SoupX/obese_LB231_soupx_counts")
LB235.data <- Read10X (data.dir = "SoupX/obese_LB235_soupx_counts")

AB21B6674.data <- Read10X (data.dir = "SoupX/NAFL_21B6674_soupx_counts")
AB54.data <- Read10X (data.dir = "SoupX/NAFL_AB54_soupx_counts")
LB219.data <- Read10X (data.dir = "SoupX/NAFL_LB219_soupx_counts")
LB227.data <- Read10X (data.dir = "SoupX/NAFL_LB227_soupx_counts")

AB20B6320.data <- Read10X (data.dir = "SoupX/NASH_20B6320_soupx_counts")
AB21B948.data <- Read10X (data.dir = "SoupX/NASH_21B948_soupx_counts")
AB61.data <- Read10X (data.dir = "SoupX/NASH_AB61_soupx_counts")
LB224.data <- Read10X (data.dir = "SoupX/NASH_LB224_soupx_counts")
LB234.data <- Read10X (data.dir = "SoupX/NASH_LB234_soupx_counts")


samples<-c("p12", "p14", "p16", "p45", "LB222", "LB226", "LB228", "LB231", "LB235", "AB21B6674",  "AB54",
           "LB219", "LB227", "AB20B6320", "AB21B948", "AB61", "LB224", "LB234")

for (i in 1:length(samples)) {
  assign(samples[i], CreateSeuratObject(counts = get(paste(samples[i], ".data", sep="")) , min.cells = 3, min.features  = 200, project = samples[i], assay = "RNA"))
}

rm(list=ls(pattern = ".data"))
##Add metadata

All.list<-samples
list.lean<-c("p12", "p14", "p16", "p45")
list.obese<-c( "LB222", "LB226", "LB228", "LB231", "LB235")
list.MASL<-c("AB21B6674", "AB54","LB219", "LB227")
list.MASH<-c( "AB20B6320", "AB21B948", "AB61", "LB224", "LB234")

for (i in list.lean){
  .GlobalEnv[[i]] <- AddMetaData(object = .GlobalEnv[[i]],metadata = "lean", col.name = "disease")}
for (i in list.obese){
  .GlobalEnv[[i]] <- AddMetaData(object = .GlobalEnv[[i]],metadata = "obese", col.name = "disease")}
for (i in list.MASL){
  .GlobalEnv[[i]] <- AddMetaData(object = .GlobalEnv[[i]],metadata = "MASL", col.name = "disease")}
for (i in list.MASH){
  .GlobalEnv[[i]] <- AddMetaData(object = .GlobalEnv[[i]],metadata = "MASH", col.name = "disease")}

for (i in list.lean){
  .GlobalEnv[[i]] <- AddMetaData(object = .GlobalEnv[[i]],metadata = "lean", col.name = "bmi_range")}
for (i in c(list.obese, list.MASL, list.MASH)){
  .GlobalEnv[[i]] <- AddMetaData(object = .GlobalEnv[[i]],metadata = "obese", col.name = "bmi_range")}


#Add mitochondria and cellcyle info
for (i in All.list){
  .GlobalEnv[[i]][["percent.mito"]]<-PercentageFeatureSet(.GlobalEnv[[i]], pattern = "^MT-")
}
for (i in All.list){
  .GlobalEnv[[i]] <- CellCycleScoring(object = .GlobalEnv[[i]], g2m.features = cc.genes$g2m.genes,
                                      s.features = cc.genes$s.genes)}

#decide cutoffs
mitovsfeature.list<-list()
for (i in All.list){
  sample<-.GlobalEnv[[i]]
  mitovsfeature.list[[i]]<-FeatureScatter(object = sample, feature1 = "nFeature_RNA", feature2 = "percent.mito")+ggtitle(i)
}
mitovsfeature<-wrap_plots(mitovsfeature.list)

#run QC
for (i in All.list){
  .GlobalEnv[[i]]<-subset(x = .GlobalEnv[[i]], subset = nFeature_RNA > 700 & nFeature_RNA < 9000 & percent.mito >  -Inf & percent.mito < 30 )
  
  ##Normalize
  .GlobalEnv[[i]]<-NormalizeData(object = .GlobalEnv[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  
  ##Find variable genes across the single cells (~2000genes, otherwise remove nfeatures and look with length(nSAH@assays$RNA@var.features) )
  .GlobalEnv[[i]]<-FindVariableFeatures(object = .GlobalEnv[[i]], verbose=F)
}

#remove doublets based on scDblFinder
remove_doublets<-function(file){
  sub<-get(file)
  DefaultAssay(sub)<-"RNA"
  sce<-as.SingleCellExperiment(sub)
  sce <- scDblFinder(sce, dims = 50)
  sub@meta.data<-as.data.frame(sce@colData)
  sub2<-subset(sub,subset = scDblFinder.class %in% c("singlet"))
  sub2<-FindVariableFeatures(object = sub2, verbose=F)
  return(sub2)
}

library(scDblFinder)
for (i in All.list){
  .GlobalEnv[[i]]<-remove_doublets(file=i)
}

#integration
objectlist<-mget(samples)
disease.anchors <- FindIntegrationAnchors(object.list = objectlist,dims = 1:50)
snNAFLDLiver <- IntegrateData(anchorset =disease.anchors, dims = 1:50)

#PCA, UMAP and clustering
DefaultAssay(snNAFLDLiver)<-"RNA"
snNAFLDLiver<-FindVariableFeatures(object = snNAFLDLiver, verbose=F)
DefaultAssay(snNAFLDLiver)<-"integrated" ##use only for "clustering" not differential analysis
snNAFLDLiver <- ScaleData(object = snNAFLDLiver,
                          vars.to.regress = c("nCount_RNA", "percent.mito"))
snNAFLDLiver <- RunPCA(object = snNAFLDLiver, npcs = 50)
snNAFLDLiver <- RunUMAP(snNAFLDLiver, dims = 1:50)
snNAFLDLiver <- FindNeighbors(snNAFLDLiver, reduction = "pca", dims = 1:50)
snNAFLDLiver <- FindClusters(snNAFLDLiver, resolution =1.5, algorithm = 1)

#markers for clustering
Hepatocyte.genes<-c("ACSS2","AKR1C1","ALB","ALDH6A1","BCHE","CYP2A6","CYP2A7","CYP3A7","G6PC","GHR","GSTA2","HAMP","HMGCS1","HPR","HSD11B1","MASP2","PCK1","RND3","RPP25L","SCD","SEC16B","SLBP","TM7SF2")
Cholangiocyte.genes<-c("EPCAM","KRT19","KRT7","SOX9","FXYD2","CLDN4","CLDN10","MMP7","CXCL1","CFTR","CD24", "BICC1")
Tcell.genes<-c("CCL5","CD2","CD3D","GNLY","GZMK","GZMB","HMGB2","PTGDS","STMN1","TRAC","TRDC","TYMS")
NK.genes<-c("CD7","NKG7","KLRD1","NCR1","NCAM1","CMC1","XCL2")
Bcell.genes<-c("CD37","CD79B","IGHG1","IGHG2","IGKC","IGLC2","LTB","MS4A1")
Endothelial.genes<-c("CCL14","CLEC14A","CLEC1B","DNASE1L3","FCN2","INMT","LIFR","MGP","RAMP2","RAMP3","S100A13","SPARCL1","TM4SF1")
Stellate.genes<-c("ACTA2","COL1A1","RBP1","TAGLN","COL1A2","COL3A1","SPARC","DCN","MYL9")
Myeloid.genes<-c("CD68","CD5L","HLA-DPB1","LYZ","MARCO","S100A8","S100A9","VSIG4","CPVL","CD14")
Plasma.genes<-c("FCRL5", "JSRP1")
pDC.genes<-c("LILRA4", "PTCRA", "CLEC4C", "LRRC26")
marker.genes<-c(Hepatocyte.genes,Cholangiocyte.genes,Stellate.genes, Myeloid.genes, Endothelial.genes, Tcell.genes, NK.genes, Bcell.genes, Plasma.genes, pDC.genes )


##Rename
Idents(snNAFLDLiver)<-snNAFLDLiver$integrated_snn_res.1.5
snNAFLDLiver<-RenameIdents(object = snNAFLDLiver, "0"="Hepatocyte",
                           "1"="Hepatocyte",
                           "2"="Hepatocyte",
                           "3"="Hepatocyte",
                           "4"="Hepatocyte",
                           "5"="Hepatocyte",
                           "6"="Lymphocyte",
                           "7"="Hepatocyte",
                           "8"="Cholangiocyte",
                           "9"="Hepatocyte",
                           "10"="Endothelial",
                           "11"="Myeloid",
                           "12"="Mesenchymal",
                           "13"="Endothelial",
                           "14"="Lymphocyte",
                           "15"="Endothelial",
                           "16"="Myeloid",
                           "17"="Myeloid",
                           "18"="Lymphocyte",
                           "19"="Endothelial",
                           "20"="Endothelial",
                           "21"="Hepatocyte",
                           "22"="Mesenchymal",
                           "23"="Mesenchymal",
                           "24"="Endothelial",
                           "25"="Lymphocyte",
                           "26"="Hepatocyte",
                           "27"="Mesenchymal",
                           "28"="Cholangiocyte",
                           "29"="Myeloid",
                           "30"="Endothelial",
                           "31"="lowQ",
                           "32"="Endothelial",
                           "33"="Endothelial",
                           "34"="Endothelial",
                           "35"="lowQ",
                           "36"="Hepatocyte",
                           "37"="B cell",
                           "38"="Cholangiocyte",
                           "39"="Lymphocyte",
                           "40"="Endothelial",
                           "41"="B cell",
                           "42"="Endothelial",
                           "43"="Myeloid",
                           "44"="Mesenchymal",
                           "45"="Hepatocyte",
                           "46"="B cell",
                           "47"="Myeloid",
                           "48"="Cholangiocyte",
                           "49"="Hepatocyte")
snNAFLDLiver[["Clusters"]] <- Idents(object = snNAFLDLiver)
Idents(snNAFLDLiver)<-snNAFLDLiver$Clusters

#subcluster and continue analysis for each Cluster

subclustering<-function(subsets, dims, group){
  sub<-subset(snNAFLDLiver, invert=F, subset = Clusters %in% subsets)
  sub@meta.data$Clusters <- droplevels(sub@meta.data$Clusters)
  DefaultAssay(sub)<-"RNA"
  sub<-FindVariableFeatures(object = sub, verbose=F)
  DefaultAssay(sub)<-"integrated" ##use only for "clustering" not differential analysis
  sub <- ScaleData(object = sub,
                   vars.to.regress = c("nCount_RNA", "percent.mito"))
  sub <- RunPCA(object = sub, npcs = dims)
  sub <- RunUMAP(sub, dims = 1:dims, return.model = T)
  sub <- FindNeighbors(sub, reduction = "pca", dims = 1:dims)
  sub <- FindClusters(sub, resolution = 2, algorithm = 1)
  DefaultAssay(sub)<-"RNA"
  saveRDS(sub, file = paste(getwd(),"/data/",group,".rds", sep = ""))
  return(sub)
}

Hepato_Chol<-subclustering(subset=c("Hepatocyte", "Cholangiocyte"), dims = 50, group = "Hepato_Chol")
Endothelial<-subclustering(subset="Endothelial", dims = 30, group="Endothelial")
Mesenchymal<-subclustering(subset="Mesenchymal", dims = 30, group="Mesenchymal")
Myeloid<-subclustering(subset="Myeloid", dims = 30, group="Myeloid")
Lymphocyte<-subclustering(subset=c("Lymphocyte"), dims = 30, group="Lymphocyte")
Bcell<-subclustering(subset=c("B cell"), dims = 15, group="B_cell")

#after cleaning of the individual Clusters, add labels to the metadata and remove all manually identifed lowQ and doublet cells

Liver_clean<-subset(snNAFLDLiver, invert=T, subset=Cluster %in% "Contamination" ) #removal of lowQ and doublets together
Liver_clean$Cluster<-factor(x=Liver_clean$Cluster, levels = c("Hepatocytes",  "Cholangiocytes",    "Endothelial", "Mesenchymal", "Lymphocytes","Myeloid-like cells", "B cells", "pDC"))
Idents(Liver_clean)<-Liver_clean$Cluster
DefaultAssay(Liver_clean)<-"RNA"
Liver_clean<-FindVariableFeatures(object = Liver_clean, verbose=F)
DefaultAssay(Liver_clean)<-"integrated" ##use only for "clustering" not differential analysis
Liver_clean <- ScaleData(object = Liver_clean,
                         vars.to.regress = c("nCount_RNA", "percent.mito"))
Liver_clean <- RunPCA(object = Liver_clean, npcs = 50)
Liver_clean <- RunUMAP(Liver_clean, dims = 1:50)

DefaultAssay(Liver_clean)<-"RNA"
