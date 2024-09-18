
#======== alternative method for integration ==============================


# make a list
all.list<-list(icam1_QC,wt1_QC,wt1_QC)
names(all.list)<-c("icam1","wt1","wt1r")

# merge
merged_Seurat <- merge(all.list[[1]], y = c(all.list[c(-1)]),project = "mice")

tail(merged_Seurat[[]])

merged_Seurat

DefaultAssay(merged_Seurat)<-"RNA"

#  Log Normalization
merged_Seurat <- NormalizeData(merged_Seurat, normalization.method = "LogNormalize", 
                               scale.factor = 10000)

merged_Seurat <- FindVariableFeatures(merged_Seurat, selection.method = "vst", 
                                      nfeatures = 3000)

# filter variable genes

list.features<-VariableFeatures(merged_Seurat)
newlist.features<-unique(c(list.features[-grep("^Hsp|^Gm|-as1$|^Hla-|^Ig[jkl]|^Igh[vjdc]|^Trb[vjdc]|^Tra[vjdc]|^Rna|^mt|^Rpl|^Rps|^Malat|^Neat|^Acyb|^Linc|G0s2",list.features)]))

VariableFeatures(merged_Seurat)<-newlist.features

# scale data

merged_Seurat <- ScaleData(merged_Seurat)

# PCA 

merged_Seurat <- RunPCA(merged_Seurat, 
                        features = VariableFeatures(object = merged_Seurat),
                        npcs = 30)

#== cluster directlry
n=30
merged_Seurat <- RunUMAP(merged_Seurat,dims = 1:n, reduction = "pca")
merged_Seurat <- FindNeighbors(merged_Seurat, reduction = "pca", dims = 1:n)
merged_Seurat <- FindClusters(merged_Seurat, resolution = 0.4)
DimPlot(merged_Seurat,split.by = "orig.ident",cols = c25) # save pdf send to me

#==== Harmony ========

merged_Seurat <- RunHarmony(merged_Seurat, group.by.vars=c("orig.ident"))
merged_Seurat <- RunUMAP(object=merged_Seurat, reduction="harmony", 
                         dims=1:n, verbose=FALSE,n.components = 5)
merged_Seurat <- FindNeighbors(object=merged_Seurat, reduction="harmony", dims=1:n, verbose=FALSE)
merged_Seurat <- FindClusters(object=merged_Seurat, resolution=.4, verbose=T)
DimPlot(merged_Seurat,split.by = "orig.ident",cols = c25) # save PDF


FeaturePlot(object = merged_Seurat, pt.size = 2,ncol = 3,raster = T,
            features = c("Il1b","Icam1","Prg4","Mmp3","Cxcl2","Col1a1"),
            order = TRUE,
            min.cutoff = 'q20', max.cutoff = "q99",
            label = F,cols = c("gray","darkred"),
            repel = TRUE
) 

