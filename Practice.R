library(Seurat)
data("pbmc_small")

NormalizeData(pbmc_small)
pbmc_small=NormalizeData(pbmc_small)
ScaleData(pbmc_small)
pbmc_small=ScaleData(pbmc_small)
RunPCA(pbmc_small)
pbmc_small=RunPCA(pbmc_small)
RunUMAP(pbmc_small)
pbmc_small=RunUMAP(pbmc_small,dims = 1:2)
UMAPPlot(pbmc_small)
FindNeighbors(pbmc_small)
pbmc_small=FindNeighbors(pbmc_small,reduction = "umap",dims = 1:2)
FindClusters(pbmc_small)
pbmc_small=FindClusters(pbmc_small,resolution = 0.8)
UMAPPlot(pbmc_small)

# More Practice from https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- Read10X(data.dir = "hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# Shifts the expression of each gene, so that the mean expression across cells is 0
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
# JackStrawPlot Function 
#'Significant' PCs will show and the solid curve is lowest p-value
# This can take a very long time to run for large data sets
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
# Elbow Plot- a ranking of principle components based on the percentage of variance  
# 'elbow' around PC9-10; majority of true signal is captured in the first 10 PCs.
ElbowPlot(pbmc)
# Clustering Cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "CaiLabRotation")
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
# plotting the top 20 markers (or all markers if less than 20) for each cluster
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Assign Cell Types!
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# library(clustermole)
#clustermole_markers("hs")
#markers <- clustermole_markers("hs")
rownames(pbmc)


de = (FindMarkers(object = pbmc, ident.1 = "CD8 T"))
gcol= ifelse(abs(de$avg_logFC)>1 & de$p_val_adj<0.05, yes = "red", no = "black")
plot(de$avg_logFC,-log10(de$p_val), col=gcol, pch=16)

#markers <- markers[markers$celltype == "B cells",]
m#arkers <- list(B=markers$gene)

library(fgsea)
markers=gmtPathways("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/markerGenes/hsaPanglaoDB.gmt")
fc=de$avg_logFC
names(fc)=rownames(de)

library(ggplot2)
e=fgsea(markers,fc,1000)
plotEnrichment(markers$`T cells`,fc) + xlab("Gene Rank") + ylab("Enrichment Score") + labs(title = "T Cells", subtitle = "FDR = 4.850000e-09", caption = "CCL5, GZMK, GZMH, NKG7, GZMA")


# ALLIE TEST
source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R")
AllieTest <- Read10X("filtered_feature_bc_matrix")
AllieTest =scQC(AllieTest)
AllieTest <- CreateSeuratObject(counts = AllieTest, project = "filtered_feature_bc_matrix")
#NormalizeData(AllieTest)
AllieTest=NormalizeData(AllieTest)
#FindVariableFeatures(AllieTest)
AllieTest=FindVariableFeatures(AllieTest)
#ScaleData(AllieTest)
AllieTest=ScaleData(AllieTest)
AllieTest = RunPCA(AllieTest)
PCAPlot(AllieTest)
AllieTest = RunTSNE(AllieTest)
TSNEPlot(AllieTest)
VG = VariableFeatures(AllieTest)
FeaturePlot(AllieTest, VG[1:10], order = TRUE)
AllieTest <- RunUMAP(AllieTest, dims = 1:10)
AllieTest <- FindNeighbors(AllieTest, reduction = "umap", dims = 1:2)
AllieTest <- FindClusters(AllieTest, resolution = 0.01)
head(Idents(AllieTest), 5)
DimPlot(AllieTest, reduction = "umap")



library(fgsea)
de = (FindMarkers(object = AllieTest, ident.1= "0"))
markers=gmtPathways("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/markerGenes/mmuPanglaoDB.gmt")
fc=de$avg_logFC
names(fc)=rownames(de)
gcol= ifelse(abs(de$avg_logFC)>1 & de$p_val_adj<0.05, yes = "red", no = "black")
plot(de$avg_logFC,-log10(de$p_val), col=gcol, pch=16)
e=fgsea(markers,fc,1000)
plotEnrichment(markers$`Fibroblasts`,fc) + xlab("Gene Rank") + ylab("Enrichment Score") + labs(title = "Fibroblasts", subtitle = "FDR = 3.020000e-09", caption = "Gsn, Mgp, Fn1, Lum, Col3a1")               
               
               
               
               