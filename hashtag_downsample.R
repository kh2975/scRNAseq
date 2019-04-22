install.packages('devtools')
install.packages('raster')
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(raster)
#install.packages("reticulate")
#library(reticulate)

# create a new environment 
#conda_create("r-reticulate")

# install SciPy
#conda_install("r-reticulate", "scipy")

# import SciPy (it will be automatically discovered in "r-reticulate")
#scipy <- import("scipy")

#use_python("/usr/local/bin/python")
#py_install("umap-learn")

umi.data <- Read10X(data.dir = "/Users/krystenharvey/Documents/Aifantis_Lab/filtered_gene_matrix/")

#umi.data <- CreateSeuratObject(counts = umi.data, project = "aml", min.cells = 3, min.features = 200)
umi.data
hto.data <- read.table(file = "/Users/krystenharvey/Documents/Aifantis_Lab/HTO_file/Hashtag_results.tsv", as.is = TRUE)
hto.data

joint.bcs <- intersect(x = colnames(x = umi.data), y = colnames(x = hto.data))

umi.data <- (x= umi.data[, joint.bcs])
hto.data <- as.matrix(x = hto.data[, joint.bcs])

rownames(x = hto.data)
rownames(x = umi.data)


# Setup Seurat object
hashtags <- CreateSeuratObject(counts = umi.data)


# Normalize RNA data with log normalization
hashtags <- NormalizeData(object = hashtags)
# Find and scale variable features
hashtags <- FindVariableFeatures(object = hashtags, selection.method = "mean.var.plot")
hashtags <- ScaleData(object = hashtags, features = VariableFeatures(object = hashtags))


hashtags[["HTO"]] <- CreateAssayObject(counts = hto.data)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
hashtags <- NormalizeData(object = hashtags, assay = "HTO", normalization.method = "CLR")

hashtags <- HTODemux(hashtags, assay = "HTO", positive.quantile = 0.99)

table(hashtags$HTO_classification.global)


FeatureScatter(object = hashtags, feature1 = "Hashtag2-TGATGGCCTATTGGG", feature2 = "Hashtag3-TTCCGCCTCTCTTTG")

Idents(object = hashtags) <- "HTO_classification.global"
VlnPlot(object = hashtags, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


hashtag.subset <- subset(x = hashtags, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(x = dist(x = t(x = GetAssayData(object = hashtag.subset, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
pbmc.hashtag.subset <- RunTSNE(object = hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(object = pbmc.hashtag.subset)



HTOHeatmap(object = hashtags, assay = "HTO")

singlets <- subset(x = hashtags, idents = "Singlet")

# Select the top 1000 most variable features
singlets <- FindVariableFeatures(object = singlets, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
singlets <- ScaleData(object = singlets, features = VariableFeatures(object = singlets))

# Run PCA
singlets <- RunPCA(object = singlets, features = VariableFeatures(object = singlets))

singlets <- FindNeighbors(object = singlets, reduction = "pca", dims = 1:10)
singlets <- FindClusters(object = singlets, resolution = 0.6)


singlets <- RunTSNE(object = singlets, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
DimPlot(object = singlets, group.by = "HTO_classification")


cntrl <- singlets[,singlets$HTO_classification=="Hashtag2-TGATGGCCTATTGGG"]
aml<- singlets[,singlets$HTO_classification=="Hashtag3-TTCCGCCTCTCTTTG"]

cntrl
aml
install.packages('cowplot')
library(cowplot)
ctrl <- cntrl
aml <- aml

ctrl$stim <- "CTRL"
ctrl <- subset(x = ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(object = ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
aml$stim <- "AML"
aml <- subset(x = aml, subset = nFeature_RNA > 500)


immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, aml), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(object = immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(object = immune.combined, verbose = FALSE)
immune.combined <- RunPCA(object = immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(object = immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(object = immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)


#find markers
markers0 <- FindMarkers(immune.combined,0)
markers0$avg_logFC <- abs(markers0$avg_logFC)
markers0 <- markers0[order(markers0$avg_logFC),] 
l <- length(markers0$avg_logFC)
fin_markers0 <- markers0[l:(l-20),]

markers1 <- FindMarkers(immune.combined, 1)
markers1$avg_logFC <- abs(markers1$avg_logFC)
markers1 <- markers1[order(markers1$avg_logFC),] 
l <- length(markers1$avg_logFC)
fin_markers1 <- markers1[l:(l-20),]


markers2 <- FindMarkers(immune.combined, 2)
markers2$avg_logFC <- abs(markers2$avg_logFC)
markers2 <- markers2[order(markers2$avg_logFC),] 
l <- length(markers2$avg_logFC)
fin_markers2 <- markers2[l:(l-20),]

markers3 <- FindMarkers(immune.combined, 3)
markers3$avg_logFC <- abs(markers3$avg_logFC)
markers3 <- markers3[order(markers3$avg_logFC),] 
l <- length(markers3$avg_logFC)
fin_markers3 <- markers3[l:(l-20),]

markers4 <- FindMarkers(immune.combined, 4)
markers4$avg_logFC <- abs(markers4$avg_logFC)
markers4 <- markers4[order(markers4$avg_logFC),] 
l <- length(markers4$avg_logFC)
fin_markers4 <- markers4[l:(l-20),]

markers5 <- FindMarkers(immune.combined, 5)
markers5$avg_logFC <- abs(markers5$avg_logFC)
markers5 <- markers5[order(markers5$avg_logFC),] 
l <- length(markers5$avg_logFC)
fin_markers5 <- markers5[l:(l-20),]


markers6 <- FindMarkers(immune.combined, 6)
markers6$avg_logFC <- abs(markers6$avg_logFC)
markers6 <- markers6[order(markers6$avg_logFC),] 
l <- length(markers6$avg_logFC)
fin_markers6 <- markers6[l:(l-20),]


fin_markers0
fin_markers1
fin_markers2
fin_markers3
fin_markers4
fin_markers5
fin_markers6


#plotting

p1 <- DimPlot(object = immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(object = immune.combined, reduction = "umap", split.by = "stim")

FeaturePlot(object = immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A","FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9")

immune.combined <- RenameIdents(object = immune.combined, `0` = "B", `1` = "T", `2` = "T",`3` = "T", `4` = "T", `5` = "T", `6` = "?1", `7` = "CD16 Mono", `8` = "?2", `9` = "B", `10` = "?3")

DimPlot(object = immune.combined, label = TRUE)

DimPlot(object = immune.combined, label = TRUE,split.by = "stim")

#markers

cntrl <- immune.combined[,immune.combined$stim=="CTRL"]
aml <- immune.combined[,immune.combined$stim=="AML"]



