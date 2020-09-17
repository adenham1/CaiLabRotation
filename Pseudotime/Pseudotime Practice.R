# Pseudotime Practice

#Install Monacle
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("monocle")
library(Matrix)
library(monocle)
library(Seurat)
source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R")
cMatrix <- Read10X("filtered_feature_bc_matrix")
cMatrix <- cMatrix[rowSums(cMatrix) > 0,]
cMatrix =scQC(cMatrix)
cMatrix <- CreateSeuratObject(counts = cMatrix)
cMatrix=NormalizeData(cMatrix)
cMatrix=FindVariableFeatures(cMatrix)
cMatrix=ScaleData(cMatrix)
cMatrix = RunPCA(cMatrix)
PCAPlot(cMatrix)
cMatrix = RunTSNE(cMatrix)
TSNEPlot(cMatrix)
VG = VariableFeatures(cMatrix)
cMatrix <- RunUMAP(cMatrix, dims = 1:10)
UMAPPlot(cMatrix)
cMatrix <- FindNeighbors(cMatrix, reduction = "umap", dims = 1:2)
cMatrix <- FindClusters(cMatrix, resolution = 0.01)
head(Idents(cMatrix), 5)
DimPlot(cMatrix, reduction = "umap")

c2=subset(cMatrix,ident=0)
cMatrix=c2@assays$RNA@counts
source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/getPseudotime.R")
PT= computePseudoTime(cMatrix)
write.csv(PT,"Output.csv")

c2$pseudotime = PT
FeaturePlot(c2, 'pseudotime', order = TRUE)

source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')
c2= CellCycleScoring(c2,s.features = hsa2mmu_SYMBOL(cc.genes.updated.2019$s.genes), g2m.features = hsa2mmu_SYMBOL(cc.genes.updated.2019$g2m.genes))
FeaturePlot(c2, 'Phase', order = TRUE)
boxplot(c2$pseudotime~c2$Phase)
plot(log1p(pseudoTiveV$Pseudotime))
