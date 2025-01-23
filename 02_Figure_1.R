#Figure 1 plots
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

Liver_clean <- readRDS("data/Liver_clean.rds")
Liver_clean$disease<-factor(x=Liver_clean$disease, levels = c("lean", "obese", "MASL","MASH"))
Liver_clean$Cluster<-factor(x=Liver_clean$Cluster, levels = c("Hepatocytes",  "Cholangiocytes",    "Endothelial", "Mesenchymal", "Lymphocytes","Myeloid-like cells", "B cells", "pDC"))
DefaultAssay(Liver_clean)<-"RNA"

#Figure 1
Umap_Liver<-DimPlot(Liver_clean, reduction = "umap",  group.by = "Cluster", repel = F, label = F, pt.size = 0.5, cols=colours_Liver, order = F, raster = F)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend() 

#script for DoMultiBarHeatmap.R
Idents(Myeloid_postclean)<-Myeloid_postclean$new_subcluster
Myeloid_markers<-FindAllMarkers(object = Myeloid_postclean, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
wiltopTrans <- Myeloid_markers[Myeloid_markers$p_val_adj<1e-10,]
wiltopTrans2<-wiltopTrans %>% group_by(cluster) %>% top_n(10, avg_log2FC)
Myeloid_postclean <- ScaleData(Myeloid_postclean, features = as.character(unique(wiltopTrans2$gene)), assay = "RNA")

MultiHeatmap_Myeloid_Cluster<-DoMultiBarHeatmap(subset(Myeloid_postclean, downsample=1000), features = as.character(unique(wiltopTrans2$gene)), group.by = "new_subcluster", additional.group.by = "disease", 
                                                additional.group.sort.by = "disease",assay = "RNA", disp.min = -2, disp.max = 2, cols.use = list(new_subcluster=colours_Myeloid_new,disease=colours_disease), size=6)+ scale_fill_gradientn(colors = c("steelblue2", "white", "red2"))+guides(colour="none")+FontSize(x.text = fontsize, y.text = 12, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))

Umap_Myeloid_postclean<-DimPlot(Myeloid_postclean, reduction = "umap",  group.by = "new_subcluster", repel = F, label = F, pt.size = 0.5, cols=colours_Myeloid_new, order = F, 
                                raster = F)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), 
                                                                 legend.title=element_text(size=fontsize ))+NoLegend()
Umap_Myeloid_disease<-DimPlot(Myeloid_downstream, reduction = "umap",  group.by = "new_subcluster", split.by = "disease", repel = F, label = F, pt.size = 0.5, cols=colours_Myeloid_new, order = F, raster = F)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend() 

##boxplot percentage patients
Idents(Myeloid_postclean)<-Myeloid_postclean@meta.data$orig.ident
table<-table(Myeloid_postclean@meta.data$new_subcluster, Myeloid_postclean@active.ident)
table<-as.data.frame.matrix(table)
table<-table[order(row.names(table)),]
table <- table[ order(as.numeric(row.names(table))), ]
table.perc<-apply(table[],2,function (x){(x/sum(x))*100})
table.perc.round<-as.data.frame(t(round(table.perc, digits = 2))) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-c("lean", "lean", "lean", "lean","obese", "obese","obese", "obese","obese","MASL",
                            "MASL","MASL","MASL","MASH","MASH","MASH","MASH","MASH")
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("lean", "obese",  "MASL","MASH"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("KC", "GPNMB Mac", "TransMac","preMac",  "Monocyte",  "cDC1", "cDC2", "migDC"),  names_to = "Cluster", values_to = "Percentage")

box_clusters_Myeloid_postclean<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", palette=colours_disease, ylim=c(0,150))+theme(axis.text.x = element_text(angle = 45, hjust=1), 
                                  axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

#GSVA
library(GSVA)
library(msigdbr)
library(tidyverse)
library(cerebroApp)
library(GSEABase)
library(limma)

genesetsH <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select("gs_name","gene_symbol") %>% as.data.frame()
genesetsRea <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME") %>% dplyr::select("gs_name","gene_symbol") %>% as.data.frame()
genesetsKEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>% dplyr::select("gs_name","gene_symbol") %>% as.data.frame()
genesetsGO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select("gs_name","gene_symbol") %>% as.data.frame()

geneset<-rbind(genesetsH,genesetsRea,genesetsKEGG,genesetsGO)
genesets <- split(geneset$gene_symbol, geneset$gs_name)

Idents(Myeloid_postclean)<-Myeloid_postclean$new_subcluster
sub<-subset(Myeloid_postclean, subset=new_subcluster %in% c("KC", "GPNMB Mac","TransMac","preMac",  "Monocyte"))
Idents(sub)<-sub$new_subcluster
sub$new_subcluster<-droplevels(sub$new_subcluster)
sub$new_subcluster<-as.character(sub$new_subcluster)
sub$str_subcluster<-str_c(sub$new_subcluster,"_", sub$disease)
sub$str_subcluster<- gsub(" ", "_", sub$str_subcluster)
sub$str_subcluster<-as.factor(sub$str_subcluster)
Idents(sub)<-sub$str_subcluster
df.data <- GetAssayData(object = sub, slot = "data")
df.group <- data.frame(umi = names(Idents(sub)), 
                       cluster = as.character(sub@meta.data$str_subcluster), 
                       stringsAsFactors = F)
writeClipboard(unique(paste(shQuote(readClipboard()), collapse=", ")))
gsva_terms=c("REACTOME_INFLAMMASOMES","GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
             "GOBP_INNATE_IMMUNE_RESPONSE",
             "GOBP_MACROPHAGE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
             "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE",
             "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
             "GOBP_CHEMOKINE_PRODUCTION",
             "REACTOME_INTERLEUKIN_10_SIGNALING",
             "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",
             "REACTOME_SIGNALING_BY_VEGF",
             "HALLMARK_FATTY_ACID_METABOLISM",
             "HALLMARK_GLYCOLYSIS",
             "HALLMARK_APOPTOSIS",
             "GOBP_PHAGOCYTOSIS",
             "REACTOME_AUTOPHAGY",
             "HALLMARK_HYPOXIA")
picked_gsva_terms<-genesets[gsva_terms]
gsvascore_all  = gsva(data.matrix(df.data), picked_gsva_terms, method="ssgsea", parallel.sz=10)

clusters<-levels(sub)
mean_gsva<-list()
for(i in 1:length(clusters)){
  mean_gsva[[i]]<-as.data.frame(rowMeans(gsvascore_all[,subset(df.group, cluster == clusters[i], select = umi)$umi]))
  names(mean_gsva[[i]])<-clusters[i]
}
mean_gsva<-do.call(cbind, mean_gsva)

rownames(mean_gsva)<- gsub("HALLMARK", "", rownames(mean_gsva))
rownames(mean_gsva)<- gsub("GOBP", "", rownames(mean_gsva))
rownames(mean_gsva)<- gsub("REACTOME", "", rownames(mean_gsva))
rownames(mean_gsva)<- gsub("_", " ", rownames(mean_gsva))

Myeloid_annotation_bar<-data.frame(celltype=rep(c("KC", "GPNMB_Mac", "TransMac","preMac",  "Monocyte"),c(4,4,4,4,4)),disease=unlist(rep(list(c("lean", "obese", "MASL","MASH")),5)))
row.names(Myeloid_annotation_bar) <- str_c(Myeloid_annotation_bar$celltype,"_", Myeloid_annotation_bar$disease)
Myeloid_annotation_bar$disease <- factor(Myeloid_annotation_bar$disease, levels = c("lean", "obese","MASL","MASH"))
Myeloid_annotation_bar<-Myeloid_annotation_bar[,c(2,1)]

mean_gsva<-mean_gsva[,row.names(Myeloid_annotation_bar)]

my_colour = list(
  celltype = c(KC = "slategray3", GPNMB_Mac="#00B6EB",TransMac = "navy", preMac = "darkgreen",Monocyte="darkmagenta"),
  disease = c(lean = "palegreen4", obese = "steelblue3", MASL="sienna3", MASH="red4"))

heatmap_gsva2<-pheatmap(mean_gsva,scale="row", cluster_cols = F, cluster_rows = F,annotation_colors = my_colour ,annotation_col = Myeloid_annotation_bar, show_rownames = T,show_colnames = F,annotation_legend = F, fontsize = fontsize, treeheight_row = 1,treeheight_col = 2,gaps_col = c(4,8,12,16), color =colorRampPalette( c("steelblue2", "white", "red2"))(50), angle_col = 90 )

#Figure S1
Umap_Liver_orig<-DimPlot(Liver_clean, reduction = "umap",  group.by = "orig.ident", repel = F, label = F, pt.size = 0.5,  raster = T, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))
Umap_Liver_disease<-DimPlot(Liver_clean, reduction = "umap",  group.by = "disease", repel = F, label = F, pt.size = 0.5, cols=colours_disease, raster = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend() 

##DotPlot
dotplot_Liver_genes<-c("CYP3A4", "SDS", "ASGR1", #Hep
                       "BICC1", "PKHD1", "DCDC2", #Chol
                       "PTPRB",  "FLT1", "NRG3",#Endo
                       "LAMA2", "PTH1R", "C7", #mesen
                       "THEMIS", "CD96", "NCAM1",#Lymph
                       "CD163",  "MARCO", "FCGR3A",#Myeloid
                       "BANK1",  "MS4A1", "CD79A", #Bcells 
                       "CLEC4C", "P2RY14" #pDC
)
dotplot_Liver<-DotPlot(Liver_clean, group.by= "Cluster",features = rev(dotplot_Liver_genes) , assay = "RNA", cols = c("steelblue2", "red"))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) +xlab("")+ylab("")+coord_flip()

Umap_Liver_disease_split<-DimPlot(Liver_clean, reduction = "umap",  group.by = "Cluster", split.by = "disease", repel = F, label = F, pt.size = 0.5, cols=colours_Liver, raster = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend() 

#miloR
#miloR
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
sub<-Myeloid_postclean
Idents(sub)<-sub$new_subcluster
dims=30
DefaultAssay(sub)<-"RNA"
MILO  <- as.SingleCellExperiment(sub)
reducedDim(MILO, "PCA", withDimnames=TRUE) <- sub[['pca']]@cell.embeddings
reducedDim(MILO, "UMAP", withDimnames=TRUE) <- sub[['umap']]@cell.embeddings
MILO_obj <- Milo(MILO)
MILO_obj <- buildGraph(MILO_obj, k = 30, d = dims, reduced.dim = "PCA")
MILO_obj <- makeNhoods(MILO_obj, prop = 0.2, k = 30, d=dims, refined = TRUE, reduced_dims = "PCA", refinement_scheme="graph")
plotNhoodSizeHist(MILO_obj) ##need distribution peak betwwen 50 and 100, otherwise up k and lower prop
MILO_obj <- countCells(MILO_obj, meta.data = data.frame(colData(MILO_obj)), sample="orig.ident")
#MILO_obj <- calcNhoodDistance(MILO_obj, d=dims, reduced.dim = "PCA") #~5min

milo.design <- as.data.frame(xtabs(~disease + orig.ident, data=data.frame(colData(MILO_obj))))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$orig.ident
milo.design <- milo.design[colnames(nhoodCounts(MILO_obj)),]
milo.res <- testNhoods(MILO_obj, design=~disease, design.df=milo.design, fdr.weighting="graph-overlap")
MILO_obj <- buildNhoodGraph(MILO_obj)
milo.res <- annotateNhoods(MILO_obj, milo.res, coldata_col = "new_subcluster")##change to metadata
milo.res$new_subcluster<-factor(x=milo.res$new_subcluster, levels = levels(sub))

plot<-plotNhoodGraphDA(MILO_obj, milo.res, alpha=0.05)+ggtitle("")

#miloR barplot
df<-list()
for(i in levels(Myeloid_downstream$new_subcluster)){
  df[[i]]<-mean(milo.res[milo.res$new_subcluster %in% i,]$logFC)
  names(df[[i]])<-i
}
milo_logfc<-t(do.call(cbind, df))
colnames(milo_logfc)<-"logFC"
milo_logfc<-as.data.frame(milo_logfc)
milo_logfc$name<-rownames(milo_logfc)
milo_logfc$group<-ifelse(milo_logfc$logFC <0, "low", "high")
barplot_miloR<-ggbarplot(milo_logfc, x = "name", y = "logFC",
                         fill="group",  # change fill color by mpg_level
                         color = "white",            # Set bar border colors to white
                         palette = c("red","blue"),            # jco journal color palett. see ?ggpar
                         sort.val = "none",           # Sort the value in ascending order
                         sort.by.groups = FALSE,     # Don't sort inside each group
                         x.text.angle = 90,          # Rotate vertically x axis texts
                         ylab = "logFC",
                         xlab = "",
                         legend.title = "Abundance"
)+coord_flip()+theme(axis.text.x = element_text(angle = 0, hjust=0.5))+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, 
                                                                                y.title = fontsize,legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend() 

#Quality Plots
Vlnplotstats<-c("percent.mito", "nFeature_RNA", "nCount_RNA", "S.Score","G2M.Score")
for(i in Vlnplotstats){ 
  Figure<-VlnPlot(Liver_clean, group.by= "Cluster",features = i , assay = "RNA", cols = colours_Liver, pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))+NoLegend()
  ggsave(plot=Figure, filename=paste(getwd(),"/Extra/",i,"_Liver_stats_Cluster_Violin.svg", sep = "") ,height=4, width=4, units="in", dpi=320)
}

#for cellphoneDB
Liver_clean$subcluster_cellphonedb<-Liver_clean$Cluster
Idents(Myeloid_postclean)<-Myeloid_postclean$new_subcluster

to_KC<-WhichCells(Myeloid_postclean, idents = c("KC"))
to_GPNMB<-WhichCells(Myeloid_postclean, idents = c("GPNMB Mac"))
to_TransMac<-WhichCells(Myeloid_postclean, idents = c("TransMac"))
to_preMac<-WhichCells(Myeloid_postclean, idents = c("preMac"))
to_mono<-WhichCells(Myeloid_postclean, idents = c("Monocyte"))
to_DC<-WhichCells(Myeloid_postclean, idents = c("cDC1","cDC2","migDC"))

Idents(Endothelial_postclean)<-Endothelial_postclean$subcluster
to_LSEC<-WhichCells(Endothelial_postclean, idents = c("LSEC"))
to_EC<-WhichCells(Endothelial_postclean, idents = c("hepatic artery EC","portal vein EC","central vein EC", "scar EC", "lymphatic EC"))

Idents(Mesenchymal_postclean)<-Mesenchymal_postclean$subcluster
to_HSC<-WhichCells(Mesenchymal_postclean, idents = c("HSC","Fibroblast"))
to_VSMC<-WhichCells(Mesenchymal_postclean, idents = c("VSMC"))

Idents(Lymphocyte_postclean)<-Lymphocyte_postclean$subcluster
to_MAIT<-WhichCells(Lymphocyte_postclean, idents = c("MAIT"))
to_ILC<-WhichCells(Lymphocyte_postclean, idents = c("ILC"))
to_NK<-WhichCells(Lymphocyte_postclean, idents = c("NK_resident","NK_cyto"))
to_CD4<-WhichCells(Lymphocyte_postclean, idents = c("CD4_TN","CD4_Treg","CD4_TRM","CD4_Th17"))
to_CD8<-WhichCells(Lymphocyte_postclean, idents = c("CD8_TN","CD8_TRM","CD8_Tcyto"))

meta.data<-Liver_clean@meta.data
meta.data$subcluster_cellphonedb<-as.character(meta.data$subcluster_cellphonedb)
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_GPNMB]<-"GPNMB Mac"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_KC]<-"KC"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_TransMac]<-"TransMac"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_preMac]<-"preMac"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_mono]<-"Monocyte"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_DC]<-"DCs"

meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_LSEC]<-"LSEC"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_EC]<-"EC"

meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_HSC]<-"HSC"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_VSMC]<-"VSMC"

meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_MAIT]<-"MAIT"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_ILC]<-"ILC"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_NK]<-"NK"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_CD4]<-"CD4_T"
meta.data$subcluster_cellphonedb[rownames(meta.data) %in% to_CD8]<-"CD8_T"

meta.data$subcluster_cellphonedb<-as.factor(meta.data$subcluster_cellphonedb)
Liver_clean@meta.data<-meta.data
Idents(Liver_clean)<-Liver_clean$subcluster_cellphonedb

#write files for cellphonedb
cellphone_for_seurat <- function(seu, name,disease_cluster){
  Seurat_obj<-subset(seu, subset=disease %in% disease_cluster)
  Seurat_obj<-NormalizeData(object = Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  counts <- as.data.frame(
    as.matrix(
      Seurat_obj@assays$RNA@data)
  )
  
  #colnames(counts) <- paste('d-pos_', colnames(counts), sep = '')
  
  library("biomaRt")
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  genes  <-  getBM(filters='hgnc_symbol',
                   attributes = c('ensembl_gene_id','hgnc_symbol'),
                   values = rownames(counts),
                   mart = ensembl)
  
  counts <- counts[rownames(counts) %in% genes$hgnc_symbol,]
  
  counts <- tibble::rownames_to_column(
    as.data.frame(counts), var = 'hgnc_symbol')
  
  counts <- plyr::join(counts, genes)
  
  counts$hgnc_symbol <- NULL
  
  counts <- cbind(counts[,which(colnames(counts) == 'ensembl_gene_id')], counts)
  
  colnames(counts)[1] <- 'Gene'
  counts$ensembl_gene_id <- NULL
  ## change @ seurat_cluster to ident
  metadata <- data.frame(Cell = rownames(Seurat_obj@meta.data),
                         cell_type = Seurat_obj$subcluster_cellphonedb
  )
  
  #metadata$Cell <- paste('d-pos_', metadata$Cell, sep = '')
  
  data.table::fwrite(counts,
                     file = paste0("CellPhoneDB/",disease_cluster, "/counts_",name,"_",disease_cluster, ".txt", sep=""),
                     quote = F,
                     col.names = T,
                     row.names = F,
                     sep = '\t')
  
  data.table::fwrite(metadata,
                     file = paste0("CellPhoneDB/",disease_cluster, "/metadata_",name,"_",disease_cluster,".txt", sep=""),
                     quote = F,
                     col.names = T,
                     row.names = F,
                     sep = '\t')
  
  #system('cellphonedb method statistical_analysis metadata.txt counts.txt --iterations=10 --threads=2')
  
  #system('cellphonedb plot dot_plot')
  
  #system('cellphonedb plot heatmap_plot metadata.txt')
}

cellphone_for_seurat(seu=Liver_clean,name = "Liver", disease_cluster = "lean" )
cellphone_for_seurat(seu=Liver_clean,name = "Liver", disease_cluster = "obese" )
cellphone_for_seurat(seu=Liver_clean,name = "Liver", disease_cluster = "MASL" )
cellphone_for_seurat(seu=Liver_clean,name = "Liver", disease_cluster = "MASH" )

#run default cellphonedb
#plot cellphonedb
heatmap_cellphonedb_diseases<-function(diseases, ordering){
  mypvals <- read.table(paste0("CellPhoneDB/",diseases,"/pvalues.txt", sep=""),
                        header = T,sep = "\t",stringsAsFactors = F)
  mymeans <- read.table(paste0("CellPhoneDB/",diseases,"/means.txt", sep=""),
                        header = T,sep = "\t",stringsAsFactors = F) 
  mypvals<-mypvals %>% distinct(interacting_pair, .keep_all = T) ##to make unique and get correct counting
  
  colnames(mypvals)<-gsub("B.cells", "B_cells", colnames(mypvals))
  colnames(mymeans)<-gsub("B.cells", "B_cells", colnames(mymeans))
  
  colnames(mypvals)<-gsub("GPNMB.Mac", "GPNMB_Mac", colnames(mypvals))
  colnames(mymeans)<-gsub("GPNMB.Mac", "GPNMB_Mac", colnames(mymeans))
  
  #colnames(mypvals)[12:ncol(mypvals)] 
  sm = as.data.frame(
    do.call(rbind,
            lapply( 12:ncol(mypvals) , function(i){
              return(c( strsplit(colnames(mypvals)[i],'\\.')[[1]],
                        sum(mypvals[,i] <0.05)))
            }))
  )
  #head(sm)
  colnames(sm)=c('SOURCE' ,'TARGET' ,'count')
  sm$count = as.numeric( sm$count )
  #sm
  #write.table(sm,file = 'count_network.txt',
  #            sep = '\t',
  #            quote = F,row.names = F)
  
  
  library(reshape2)
  sm_df =dcast(as.data.frame(sm),SOURCE~TARGET )
  sm_df[is.na(sm_df)]=0
  #sm_df
  rownames(sm_df) = sm_df[,1]
  sm_df = sm_df[,-1]
  sm_df2<-sm_df[,ordering]
  sm_df2<-sm_df2[ordering,]
  p1<-pheatmap::pheatmap(sm_df2, show_rownames = T, show_colnames = T, scale="none", cluster_cols = F,
                         border_color="white", cluster_rows = F,fontsize=fontsize, fontsize_row = fontsize, fontsize_col = fontsize,
                         main = diseases, treeheight_row = 0, family = "Arial",color = colorRampPalette(c("dodgerblue4","peachpuff","deeppink4" ))(25), breaks = seq(0, 25, by = 1) ,
                         treeheight_col = 0, display_numbers = F, number_color = "white")
  #p1=pheatmap::pheatmap(sm_df2,display_numbers = T, cluster_cols = F, cluster_rows = F, main=diseases,
  #                     color =colorRampPalette( c("steelblue2", "white", "red2"))(30), breaks = seq(0, 30, by = 1))
  ggsave(plot=p1, filename=paste0("CellPhoneDB/Results/",diseases,"_Heatmap.png", sep="") ,height=7, width=7, units="in", dpi=320)
  ggsave(plot=p1, filename=paste0("CellPhoneDB/Results/",diseases,"_Heatmap.svg", sep="") ,height=7, width=7, units="in", dpi=320)
  
  Under005<-mypvals[,12:ncol(mypvals)]<0.05
  ANY <- apply(Under005, 1, any)
  dfUnder005<-mypvals[ANY,]
  write.xlsx(dfUnder005, file=paste0("CellPhoneDB/Results/",diseases,"_Dotplot_genes_sig.xlsx", sep=""),
             col.names=TRUE, row.names=TRUE)
  
  GPNMB_dot<- dfUnder005[,c(1:12,grep("GPNMB", colnames(dfUnder005)))]
  write.xlsx(GPNMB_dot, file=paste0("CellPhoneDB/Results/",diseases,"GPNMB_Macs_Dotplot_genes_sig.xlsx", sep=""),
             col.names=TRUE, row.names=TRUE)
}

diseases_levels<-c("lean", "obese", "MASL", "MASH")
order<-c("Hepatocytes","Cholangiocytes","LSEC","EC","HSC","VSMC","CD4_T","CD8_T","MAIT","ILC","NK","B_cells","pDC", "KC", "GPNMB_Mac", "TransMac","preMac",  "Monocyte", "DCs")
for(i in diseases_levels){
  heatmap_cellphonedb_diseases(diseases=i,  ordering=order)}





