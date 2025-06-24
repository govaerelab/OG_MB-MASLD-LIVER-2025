#CosMx analysis
library(Seurat) #v5
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize= 8912896000)
setwd("C:/Users/u0125188/Desktop/Liver")
fontsize=14.5

Leuven1<-LoadNanostring("CosMX/Leuven_1/", fov = "Leuven1")
Leuven1$orig.ident<-"Leuven1"
Leuven2<-LoadNanostring("CosMX/Leuven_2/", fov = "Leuven2")
Leuven2$orig.ident<-"Leuven2"
Leuven3<-LoadNanostring("CosMX/Leuven_3/", fov = "Leuven3")
Leuven3$orig.ident<-"Leuven3"
Leuven4<-LoadNanostring("CosMX/Leuven_4/", fov = "Leuven4")
Leuven4$orig.ident<-"Leuven4"

#integrate
names<-c("Leuven1","Leuven2","Leuven3","Leuven4")
merged_CosMX <- merge(Leuven1, y = list(Leuven2,Leuven3,Leuven4), add.cell.ids =names, project = "CosMX_Liver")
DefaultAssay(merged_CosMX)<-"Nanostring"
merged_CosMX <-SCTransform(merged_CosMX, assay = "Nanostring",clip.range = c(-10, 10))
merged_CosMX <- RunPCA(merged_CosMX, npcs = 30)
merged_CosMX <- RunUMAP(merged_CosMX, dims = 1:30, min.dist = 0.1)
merged_CosMX <- FindNeighbors(merged_CosMX, reduction = "pca", dims=1:30, k.param = 30)
merged_CosMX <- FindClusters(merged_CosMX, resolution = 0.3) #Louvain Clustering

#check Quality and if you need to set cutoffs before, then subset before merging with nCount_Nanostring cutoff
hist(merged_CosMX$nCount_Nanostring, breaks=500, xlim=c(0,500))

#Cluster labelling
DefaultAssay(merged_CosMX)<-"Nanostring"
merged_CosMX<-NormalizeData(object = merged_CosMX, normalization.method = "LogNormalize", scale.factor = 10000, assay = "Nanostring")
merged_CosMX<-JoinLayers(merged_CosMX, assay  = "Nanostring")
#Check cluster, Liver_genes could also be based on FindAllMarkers of the snRNA-seq data
DotPlot(merged_CosMX, features = Liver_genes , assay = "Nanostring", cols = c("steelblue2", "red"), col.min = 0, dot.min = 0.1)+coord_flip()

#Label clusters
Idents(merged_CosMX)<-merged_CosMX$SCT_snn_res.0.3
new.cluster.ids <- c("Hepatocytes", "Mesenchymal", "Myeloid-like cells", "Hepatocytes",
                     "Cholangiocytes", "Endothelial", "Hepatocytes", "Lymphocytes", "B cells", 
                     "Mesenchymal", "Mesenchymal", "B cells", "Mesenchymal", "Hepatocytes",
                     "B cells", "Hepatocytes", "B cells", "Lymphocytes", "Hepatocytes", 
                     "Hepatocytes", "B cells", "Lymphocytes", "Mesenchymal"
)
names(new.cluster.ids) <- seq(0,22)
merged_CosMX <- RenameIdents(merged_CosMX, new.cluster.ids)
merged_CosMX$Cluster<-Idents(merged_CosMX)
DimPlot(merged_CosMX, pt.size = 1, label = T)

#Afterwards we changed the name of Lymphocytes and loaded in data again
merged_CosMX<-readRDS("data/merged_CosMX.rds")
colours_CosMX<-c("lightpink3" , "slategray3","mediumpurple2","navy","#00B6EB","turquoise4", "burlywood")
Idents(merged_CosMX)<-merged_CosMX$Cluster
merged_CosMX<-RenameIdents(merged_CosMX,"Lymphocytes"="Lympho-/Granulocytes")
merged_CosMX$Cluster<-Idents(merged_CosMX)
merged_CosMX$Cluster<-factor(x=merged_CosMX$Cluster, levels = c("Hepatocytes",  "Cholangiocytes",    "Endothelial", "Mesenchymal", "Lympho-/Granulocytes","Myeloid-like cells", "B cells"))
Idents(merged_CosMX)<-merged_CosMX$Cluster
DefaultAssay(merged_CosMX)<-"Nanostring"
DimPlot(merged_CosMX, reduction = "umap",  group.by = "Cluster", repel = F, label = T,  cols=colours_CosMX, order = T)


#make FeaturePlots on UMAP
CosMX_genes<-c("CD14", "CD163", "CD5L","GPNMB", "LYZ", "MARCO", "S100A4", "SPP1", "TREM2")
for(i in CosMX_genes){
  feature_cite<-FeaturePlot(merged_CosMX, label = F, features = i, order = T, pt.size = 1.2, max.cutoff =4, min.cutoff=0, cols = c("grey", "brown4"))+theme( 
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(),axis.line = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), axis.title = element_blank())+NoLegend()+ggtitle("")
  ggsave(plot=feature_cite, filename=paste0("/data/leuven/343/vsc34335/CosMX/",i,"_CosMX_FeaturePlot.pdf", sep = "") ,height=3, width=3, units="in", dpi=320)
}

#if you want to use segmentation (wasn't good here)
#border.size = NULL and border.color = NA is important for segmentation to not have boarders arounds cells
names<-c("Leuven1", "Leuven2","Leuven3","Leuven4")
for (i in names){
  DefaultBoundary(merged_CosMX[[i]]) <- "segmentation"
}
for(i in names){ #to image
  plot<-ImageDimPlot(merged_CosMX,fov=i,group.by = "Cluster",  alpha=1,border.size = NULL,border.color = NA,dark.background = T,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = T)+scale_y_reverse()+theme(panel.grid = element_blank())
  ggsave(plot=plot, filename=paste0("/data/leuven/343/vsc34335/CosMX/",i,"_CosMX_Clustering.pdf", sep = "") ,height=20, width=20, units="in", dpi=320)
}

#use centroids, border.size=NA is important so that there are no lines around the centroids
names<-c("Leuven1", "Leuven2","Leuven3","Leuven4")
for (i in names){
  DefaultBoundary(merged_CosMX[[i]]) <- "centroids"
}
for(i in names){
  plot<-ImageDimPlot(merged_CosMX,fov=i,group.by = "Cluster",  alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F)+theme(panel.grid = element_blank())
  ggsave(plot=plot, filename=paste0("/data/leuven/343/vsc34335/CosMX/",i,"_CosMX_Clustering_centroids.pdf", sep = "") ,height=20, width=20, units="in", dpi=320)
}

#Crop ROIs for Leuven1
cropped.coords <- Crop(merged_CosMX[["Leuven1"]], y = c(40000, 125000), x = c(10000, 72000), coords = "plot")
merged_CosMX[["Leuven1lower"]] <- cropped.coords
cropped.coords <- Crop(merged_CosMX[["Leuven1"]], y = c(40000, 75000), x = c(50000, 72000), coords = "plot")
merged_CosMX[["Leuven1zoom1"]] <- cropped.coords
cropped.coords <- Crop(merged_CosMX[["Leuven1"]], y = c(60000, 75000), x = c(57000, 72000), coords = "plot")
merged_CosMX[["Leuven1zoom2"]] <- cropped.coords


#plot ROIs and adjust dot size if more zoomed in
plot1<-ImageDimPlot(merged_CosMX,fov="Leuven1lower",group.by = "Cluster", size=0.5,  alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F, axes = F)+theme(panel.grid = element_blank())
plot2<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom1",group.by = "Cluster", size=1, alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F, axes = F)+theme(panel.grid = element_blank())
plot3<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom2",group.by = "Cluster",size=1.5,  alpha=1,dark.background = T,border.size=NA,  cols =colours_CosMX,crop = T,nmols = 10000,flip_xy = F, axes = F)+theme(panel.grid = element_blank())

ggsave(plot=plot1, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_lower_CosMX_Clustering_centroids.pdf", sep = "") ,height=10, width=10, units="in", dpi=320)
ggsave(plot=plot2, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_zoom1_CosMX_Clustering_centroids_bigger.pdf", sep = "") ,height=10, width=10, units="in", dpi=320)
ggsave(plot=plot3, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_zoom2_CosMX_Clustering_centroids_bigger.pdf", sep = "") ,height=10, width=10, units="in", dpi=320)

#If you want to plot some molecules extra
genes_to_plot<-c("IL32","GPNMB")
colours_for_genes<-c("red", "green")
plot1<-ImageDimPlot(merged_CosMX,fov="Leuven1lower",group.by = "Cluster", molecules=genes_to_plot, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot2<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom1",group.by = "Cluster", molecules=genes_to_plot, size=1, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot3<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom2",group.by = "Cluster", molecules=genes_to_plot,size=1.5, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.08, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())
plot4<-ImageDimPlot(merged_CosMX,fov="Leuven1zoom3",group.by = "Cluster", molecules=genes_to_plot,size=2, alpha=0.8, mols.cols = colours_for_genes, mols.size = 0.2, cols =colours_CosMX,crop = T, border.size = NA,nmols = 10000, flip_xy = F)+theme(panel.grid = element_blank())

ggsave(plot=plot1, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_lower_CosMX_Clustering_centroids_IL32_GPNMB.pdf", sep = "") ,height=10, width=10, units="in", dpi=320)
ggsave(plot=plot2, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_zoom1_CosMX_Clustering_centroids_IL32_GPNMB.pdf", sep = "") ,height=10, width=10, units="in", dpi=320)
ggsave(plot=plot3, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_zoom2_CosMX_Clustering_centroids_IL32_GPNMB.pdf", sep = "") ,height=10, width=10, units="in", dpi=320)
ggsave(plot=plot4, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_zoom3_CosMX_Clustering_centroids_IL32_GPNMB.pdf", sep = "") ,height=10, width=10, units="in", dpi=320)

#FeaturePlot on the image, not UMAP, on ROI
genes<-c("IL32", "GPNMB","FABP5","LPL","MSR1","CD163","HLA-DRA") 
for(i in 1:length(genes)){
  plot1<- ImageFeaturePlot(merged_CosMX, fov="Leuven1zoom3",genes[i], combine=T,  coord.fixed = T, size=2,crop=T, axes=F,border.size = NA, min.cutoff=0, max.cutoff=5)+coord_flip()+theme(panel.grid = element_blank())
  ggsave(plot=plot1, filename=paste0("/data/leuven/343/vsc34335/CosMX/Leuven1_zoom3_CosMX_featureplot_",genes[i],".pdf", sep = "") ,height=4.5, width=9, units="in", dpi=320)
}

#make a DotPlot of the markers for clusters
Idents(merged_CosMX)<-merged_CosMX$Cluster
downed_cosmx<-subset(merged_CosMX, downsample=1000)
genes<-FindAllMarkers(downed_cosmx,  only.pos = T , verbose=T)
wiltopTrans <- genes[genes$p_val_adj<1e-10,]
wiltopTrans2<-wiltopTrans %>% group_by(cluster) %>% top_n(-7, p_val_adj)
plot<-DotPlot(merged_CosMX, features = as.character(unique(wiltopTrans2$gene)) ,group.by="Cluster", assay = "Nanostring", cols = c("steelblue2", "red"), col.min = 0, dot.min = 0.1)+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) +xlab("")+ylab("")
)
ggsave(plot=plot, filename="/data/leuven/343/vsc34335/CosMX/DotPlot_markers_CosMX.pdf" ,height=10, width=10, units="in", dpi=320)



