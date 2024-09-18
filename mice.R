library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd("/Users/sophiezhuang/Downloads")

Read10X_h5("lt3") -> mice_10x

mice_seurat <- CreateSeuratObject(counts = mice_10x, 
                           project = "lt3", 
                           min.cells = 5, 
                           min.features = 400)

mice_seurat[["percent.mt"]] <- PercentageFeatureSet(mice_seurat, pattern = "^mt-")
mice_seurat[["percent.Ig"]]<- PercentageFeatureSet(mice_seurat, pattern = "^Ighv|^Iglv|Igkv")
mice_seurat[["percent.Tcr"]]<- PercentageFeatureSet(mice_seurat, pattern = "^Trb[vjdc]|^Tra[vjdc]")
mice_seurat[["percent.Ribos"]] <- PercentageFeatureSet(mice_seurat, pattern = "^Rpl|^Rps")

mice_seurat$percent_lg <- apply(
  mice_seurat@assays$RNA@counts, # Inputting counts
  2, # Sum columns
  function(x)(100*max(x))/sum(x) # Find maximum expression
)

VlnPlot(mice_seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent_lg"))
quantile(mice_seurat$nCount_RNA, 0.25)




mice_seurat <- subset(mice_seurat, subset = nCount_RNA > 1000 & nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 10 & percent_lg < 20)
mice_seurat<- SCTransform(mice_seurat, verbose = FALSE)

n=50
mice_seurat<- RunPCA(mice_seurat, verbose = T)
mice_seurat <- RunUMAP(mice_seurat, dims = 1:n)
mice_seurat<- FindNeighbors(mice_seurat, reduction = "pca", dims = 1:n)
mice_seurat <- FindClusters(mice_seurat, resolution = 0.4)

DimPlot(mice_seurat)
#remove clusters 3 and 5
#mice_seurat <- subset(mice_seurat, seurat_clusters%in% c(3,5), invert=T)

DefaultAssay(mice_seurat) <- "RNA"
mice_seurat <- NormalizeData(mice_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
OA_mark <- FindAllMarkers(mice_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data", test.use = "wilcox", return.thresh = 0.25)
top10RNA <- OA_mark%>%group_by(cluster)%>%subset(avg_log2FC>0.5 & p_val_adj<0.05) %>%top_n(n = 10, wt = avg_log2FC)
top10RNA$gene
mice_seurat <- ScaleData(mice_seurat, features = top10RNA$gene)

pdf("lt3_markers rna method1.pdf", width = 18, height = 15)
DoHeatmap( mice_seurat, features = c(top10RNA$gene), slot = "scale.data", raster = F, draw.lines = T)
dev.off()

FeaturePlot(object = mice_seurat, pt.size = 2,ncol = 3,raster = T,
            features = c("Acan","Malat1","Neat1","Icam1","mt-Co3","Prg4","Lars2"),
            order = TRUE,
            min.cutoff = 'q20', max.cutoff = "q99",
            label = F,cols = c("gray","darkred"),
            repel = TRUE
)
DefaultAssay(mice_seurat)<-"SCT"
list_features<-VariableFeatures(mice_seurat) # to get the variable genes
newlist.features<-unique(c(list_features[-grep("-as1$|^Hla-|^Ig[jkl]|^Igh[vjdc]|^Trb[vjdc]|^Tra[vjdc]|^Rna|^mt|^Rpl|^Rps|^Malat|^Neat|^Acyb|^Linc|G0s2",list_features)]))

mice_seurat<- RunPCA(mice_seurat, features = newlist.features)
mice_seurat <- RunUMAP(mice_seurat, dims = 1:n)
mice_seurat<- FindNeighbors(mice_seurat, reduction = "pca", dims = 1:n)
mice_seurat <- FindClusters(mice_seurat, resolution = 0.4)
DimPlot(mice_seurat)

#heatmap
FeaturePlot(object = mice_seurat, pt.size = 2,ncol = 2,raster = T,
            features = c("nCount_RNA","nFeature_RNA","percent.mt","percent_lg"),
            order = TRUE,
            min.cutoff = 0, max.cutoff = "q99",
            label = F,cols = c("gray","darkred"),
            repel = TRUE
)

mice_seurat <- subset(mice_seurat, subset = nCount_RNA >2000 & nFeature_RNA > 500 & 
                 nFeature_RNA < 6000 & percent.mt < 10 & percent_lg < 20)
#mice_seurat <- subset(mice_seurat, seurat_clusters%in% c(3,5), invert=T)
mice_seurat<- SCTransform(mice_seurat, verbose = FALSE)

newlist.features<-unique(c(list_features[-grep("^Hsp|^Gm|-as1$|^Hla-|^Ig[jkl]|^Igh[vjdc]|^Trb[vjdc]|^Tra[vjdc]|^Rna|^mt|^Rpl|^Rps|^Malat|^Neat|^Acyb|^Linc|G0s2",list_features)]))

n=40
mice_seurat <- RunPCA(mice_seurat, features = newlist.features)
mice_seurat <- RunUMAP(mice_seurat, dims = 1:n)
mice_seurat <- FindNeighbors(mice_seurat, reduction = "pca", dims = 1:n)
mice_seurat <- FindClusters(mice_seurat, resolution = 0.4)


DimPlot(mice_seurat)
#remove clusters 3 and 5
mice_seurat<-subset(mice_seurat,seurat_clusters%in% c(3,5), invert=T)
#===== find markers
DefaultAssay(mice_seurat)<-"RNA"
mice_seurat <- NormalizeData(mice_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

OA_mark <- FindAllMarkers(mice_seurat, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25,
                          slot = "data",
                          test.use ="wilcox",
                          return.thresh = 0.25)

top10RNA <- OA_mark%>%group_by(cluster)%>%subset(avg_log2FC>0.5 & p_val_adj<0.05) %>%top_n(n = 10, wt = avg_log2FC)
top10RNA$gene


mice_seurat <- ScaleData(mice_seurat, features = top10RNA$gene)

pdf("wt3_markers rna method3.pdf",width = 18,height = 15)
DoHeatmap(mice_seurat, features = c(top10RNA$gene),slot = "scale.data",raster = F,draw.lines = T)
dev.off()


#===== find markers
DefaultAssay(mice_seurat)<-"RNA"
mice_seurat <- NormalizeData(mice_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

OA_mark <- FindAllMarkers(mice_seurat, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25,
                          slot = "data",
                          test.use ="wilcox",
                          return.thresh = 0.25)

top10RNA <- OA_mark%>%group_by(cluster)%>%subset(avg_log2FC>0.75 & p_val_adj<0.05) %>%top_n(n = 5, wt = avg_log2FC)
top10RNA$gene

mice_seurat <- ScaleData(mice_seurat, features = top10RNA$gene)

pdf("lt3_markers rna method3.pdf",width = 18,height = 15)
DoHeatmap(mice_seurat, features = c(top10RNA$gene),slot = "scale.data",raster = F,draw.lines = T)
dev.off()

saveRDS(mice_seurat, "ICAM1_goodQC.rds")
icam1_QC<-mice_seurat

saveRDS(mice_seurat, "WT1_goodQC.rds")
wt1_QC<-mice_seurat

saveRDS(mice_seurat, "LT1_goodQC.rds")
lt1_QC<-mice_seurat

saveRDS(mice_seurat, "ICAM2_goodQC.rds")
icam2_QC <- mice_seurat

saveRDS(mice_seurat, "WT2_goodQC.rds")
wt2_QC <- mice_seurat

saveRDS(mice_seurat, "ICAM3_goodQC.rds")
icam3_QC <- mice_seurat

saveRDS(mice_seurat, "WT3_goodQC.rds")
wt3_QC <- mice_seurat

saveRDS(mice_seurat, "LT2_goodQC.rds")
LT2_QC <- mice_seurat

saveRDS(mice_seurat, "LT3_goodQC.rds")
LT3_QC <- mice_seurat

#integration
DefaultAssay(icam1_QC) <- "SCT"
DefaultAssay(wt1_QC) <- "SCT"
DefaultAssay(lt1_QC) <- "SCT"
DefaultAssay(icam2_QC) <- "SCT"
#make list
all.list <- list(icam1_QC, wt1_QC)
names(all.list)<-c("icam1","wt1")
integ_features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 3000)

newlist.features<-unique(c(integ_features[-grep("^Hsp|^Gm|-as1$|^Hla-|^Ig[jkl]|^Igh[vjdc]|^Trb[vjdc]|^Tra[vjdc]|^Rna|^mt|^Rpl|^Rps|^Malat|^Neat|^Acyb|^Linc|G0s2",integ_features)]))

int.list <- PrepSCTIntegration(object.list = all.list, anchor.features = newlist.features, 
                               verbose = T)
integ_anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", anchor.features = newlist.features)

#Integrate across datasets
oa_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
n = 40
oa_integrated <- RunPCA(object = oa_integrated)
oa_integrated <- RunUMAP(oa_integrated, dims = 1:40)
oa_integrated <- FindNeighbors(oa_integrated, reduction = "pca", dims = 1:n)
oa_integrated <- FindClusters(oa_integrated, resolution = 0.4)
#plot data
DimPlot(oa_integrated, split.by = "orig.ident")
head(oa_integrated)


