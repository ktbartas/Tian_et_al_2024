# conda activate signac  
setwd('/media/emem/Lab/kbartas/gpe_cocsal')
library(Seurat)   
library(magrittr)
library(ggplot2)
library(data.table)

control1<-readRDS('R_outs/control_2023.rds')
cocaine1<-readRDS('R_outs/cocaine_2023.rds')
print('1 - data loaded')

DefaultAssay(control1) <- "RNA"
DefaultAssay(cocaine1) <- "RNA"
# run sctransform
control1 <- SCTransform(control1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
cocaine1<- SCTransform(cocaine1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
#SC transform replaces normalize, scale, and find variable features
#make list of seurat objects
plusmin.list <- list(control1,cocaine1)

#prepare for integration
features <- SelectIntegrationFeatures(object.list =plusmin.list, nfeatures = 3000)
plusmin.list <- PrepSCTIntegration(object.list = plusmin.list, anchor.features = features)

#find anchors and then integrate
plusmin.anchors <- FindIntegrationAnchors(object.list = plusmin.list, normalization.method = "SCT",
                                          anchor.features = features)
combined.sct <- IntegrateData(anchorset = plusmin.anchors, normalization.method = "SCT")
print('integrated data')
#set default to corrected integrated data
DefaultAssay(combined.sct) <- "integrated"

# Run the standard workflow for visualization and clustering
#no scaling data with SCtransform
combined.sct <- RunPCA(combined.sct, npcs = 30, verbose = FALSE);
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:20,n.neighbors=20)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:20)
combined.sct <- FindClusters(combined.sct, resolution = 0.8)

#view data
p1 <- DimPlot(combined.sct, reduction = "umap", group.by = "treatment");
p2 <- DimPlot(combined.sct, reduction = "umap", label = TRUE, repel = TRUE);

ggsave(file="neuron-figs/integrated_treatment_umap.svg", plot=p1);
ggsave(file="neuron-figs/integrated_clusters_umap.svg", plot=p2);
#see conditions side-by-side
p1<-DimPlot(combined.sct, reduction = "umap", split.by = "treatment")
ggsave(file="neuron-figs/integrated_treatment_umap-split.svg", plot=p1,width=12)
saveRDS(combined.sct,'R_outs/neuron-outs/gpe_combined-rna.rds')
print('complete')
