install.packages('devtools')
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)

install.packages('cowplot')
library(cowplot)

#Sample-Group A

umi.data <- Read10X(data.dir = "/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/A/filtered_gene_matrix/")



hto.data <- read.table(file = "/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/Hashtagsample1.tsv", as.is = TRUE)

joint.bcs <- intersect(x = colnames(x = umi.data), y = colnames(x = hto.data))

umi.data <- (x= umi.data[, joint.bcs])
hto.data <- as.matrix(x = hto.data[, joint.bcs])

length(colnames(x = hto.data))
length(colnames(x = umi.data))


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


FeatureScatter(object = hashtags, feature1 = "Hashtag1-ACCCACCAGTAAGAC", feature2 = "Hashtag2-GGTCGAGAGCATTCA", feature3="Hashtag4-AAAGCATTCTTCACG")

Idents(object = hashtags) <- "HTO_classification.global"
VlnPlot(object = hashtags, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


hashtag.subset <- subset(x = hashtags, idents = "Negative", invert=TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(x = dist(x = t(x = GetAssayData(object = hashtag.subset, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
#pbmc.hashtag.subset <- RunTSNE(object = hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
#DimPlot(object = pbmc.hashtag.subset)

HTOHeatmap(object = hashtags, assay = "HTO")

#filter singlets

singlets <- subset(x = hashtags, idents = "Singlet")

HTOHeatmap(object = singlets, assay = "HTO")


mouse_56_f1 <- singlets[,singlets$HTO_classification=="Hashtag1-ACCCACCAGTAAGAC"]
mouse_56_f2 <- singlets[,singlets$HTO_classification=="Hashtag2-GGTCGAGAGCATTCA"]
mouse_57_f1 <- singlets[,singlets$HTO_classification=="Hashtag4-AAAGCATTCTTCACG"]


mouse_56_f1 <- subset(x = mouse_56_f1, subset = nFeature_RNA > 500)
mouse_56_f2 <- subset(x = mouse_56_f2, subset = nFeature_RNA > 500)
mouse_57_f1 <- subset(x = mouse_57_f1, subset = nFeature_RNA > 500)

mouse_56_f1$stim <- "mouse_56_f1"
mouse_56_f2$stim <-"mouse_56_f2"
mouse_57_f1$stim <-"mouse_57_f1"


anchors <- FindIntegrationAnchors(object.list = list(mouse_56_f1, mouse_56_f2,mouse_57_f1), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(object = combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

p1 <- DimPlot(object = combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

plot_dir <-"/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/A/"
prostate_A_ridge_plot <- RidgePlot(hashtags, assay = "HTO", features = colnames(hashtags[["HTO"]]), ncol = 2)
ggsave(prostate_A_ridge_plot, filename = paste(plot_dir, 
                                               "prostate_A_ridge_plot.png",
                                               sep= ""), height = 7, width =7 , units = "in")

# heatmap
prostate_A_heatmap <- HTOHeatmap(hashtags, assay = "HTO", ncells = 5000)
ggsave(prostate_A_heatmap, filename =  paste(plot_dir, 
                                             "prostate_A_heatmap.png",
                                             sep= ""), height = 7, width =7 , units = "in")

# scatter plots
HTO_seurat <- as.data.frame(t(as.data.frame(hashtags@assays$HTO@data)))

install.packages('GGally')
library('GGally')
hto_pairs <- ggpairs(HTO_seurat)
hto_pairs<-pairs(HTO_seurat)
ggsave(hto_pairs, filename =  paste(plot_dir,
                                    "_hto_seurat_pairs.png",
                                    sep= ""), height = 7, width =7 , units = "in")

# plot the doublet trends
doublet_trend <- ggplot(as.data.frame(hashtags$HTO_classification), aes(hashtags$HTO_classification)) +
  geom_bar() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(doublet_trend, filename =  paste(plot_dir,
                                        "_hto_seurat_doubet.png",
                                        sep= ""), height = 7, width =7 , units = "in")




#Sample B

umi.data2 <- Read10X(data.dir = "/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/B/filtered_gene_matrix/")

hto.data2 <- read.table(file = "/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/HashtagsampleB.tsv", as.is = TRUE)

joint.bcs2 <- intersect(x = colnames(x = umi.data2), y = colnames(x = hto.data2))

umi.data2 <- (x= umi.data2[, joint.bcs2])
hto.data2 <- as.matrix(x = hto.data2[, joint.bcs2])

length(colnames(x = hto.data2))
length(colnames(x = umi.data2))


# Setup Seurat object
hashtags2 <- CreateSeuratObject(counts = umi.data2)


# Normalize RNA data with log normalization
hashtags2 <- NormalizeData(object = hashtags2)
# Find and scale variable features
hashtags2 <- FindVariableFeatures(object = hashtags2, selection.method = "mean.var.plot")
hashtags2 <- ScaleData(object = hashtags2, features = VariableFeatures(object = hashtags2))


hashtags2[["HTO"]] <- CreateAssayObject(counts = hto.data2)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
hashtags2 <- NormalizeData(object = hashtags2, assay = "HTO", normalization.method = "CLR")
hashtags2 <- HTODemux(hashtags2, assay = "HTO", positive.quantile = 0.95)

table(hashtags2$HTO_classification.global)


#FeatureScatter(object = hashtags2, feature1 = "Hashtag4-AAAGCATTCTTCACG", feature2 = "Hashtag5-CTTTGTCTTTGTGAG", feature3="Hashtag6-TATGCTGCCACGGTA")

Idents(object = hashtags2) <- "HTO_classification.global"
#VlnPlot(object = hashtags2, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


hashtag2.subset <- subset(x = hashtags2, idents = "Negative", invert=TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx2 <- as.matrix(x = dist(x = t(x = GetAssayData(object = hashtag2.subset, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
#pbmc.hashtag.subset2 <- RunTSNE(object = hashtag2.subset, distance.matrix = hto.dist.mtx2, perplexity = 100)
#DimPlot(object = pbmc.hashtag.subset2)

#HTOHeatmap(object = hashtags2, assay = "HTO")

#filter singlets

singlets2 <- subset(x = hashtags2, idents = "Singlet")

HTOHeatmap(object = singlets2, assay = "HTO")


Nd52_f1 <- singlets2[,singlets2$HTO_classification=="Hashtag4-AAAGCATTCTTCACG"]
Nd52_f2 <- singlets2[,singlets2$HTO_classification=="Hashtag5-CTTTGTCTTTGTGAG"]
Nd54_f1 <- singlets2[,singlets2$HTO_classification=="Hashtag6-TATGCTGCCACGGTA"]

saveRDS(Nd52_f1,file="/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/Nd52_f1.rds")
saveRDS(Nd52_f2,file="/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/Nd52_f2.rds")
saveRDS(Nd54_f1,file="/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/Nd54_f1.rds")

#Nd52_f1 <- subset(x = Nd52_f1, subset = nFeature_RNA > 500)
#Nd52_f2 <- subset(x = Nd52_f2, subset = nFeature_RNA > 500)
#Nd54_f1 <- subset(x = Nd54_f1, subset = nFeature_RNA > 500)

Nd52_f1$stim <- "Nd52_f1"
Nd52_f2$stim <-"Nd52_f2"
Nd54_f1$stim <-"Nd54_f1"


anchors2 <- FindIntegrationAnchors(object.list = list(Nd52_f1, Nd52_f2,Nd54_f1), dims = 1:20)
combined2 <- IntegrateData(anchorset = anchors2, dims = 1:20)

DefaultAssay(object = combined2) <- "integrated"

# Run the standard workflow for visualization and clustering
combined2 <- ScaleData(object = combined2, verbose = FALSE)
combined2 <- RunPCA(object = combined2, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
combined2 <- RunUMAP(object = combined2, reduction = "pca", dims = 1:20)
combined2 <- FindNeighbors(object = combined2, reduction = "pca", dims = 1:20)
combined2 <- FindClusters(combined2, resolution = 0.5)

p1 <- DimPlot(object = combined2, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = combined2, reduction = "umap", label = TRUE)
plot_grid(p1, p2)





plot_dir <-"/Users/krystenharvey/Documents/Aifantis_Lab/AML_healthy/prostate/B/"
prostate_A_ridge_plot <- RidgePlot(hashtags2, assay = "HTO", features = rownames(hashtags2[["HTO"]]), ncol = 2)
ggsave(prostate_A_ridge_plot, filename = paste(plot_dir, 
                                               "prostate_A_ridge_plot.png",
                                               sep= ""), height = 7, width =7 , units = "in")

# heatmap
prostate_A_heatmap <- HTOHeatmap(hashtags2, assay = "HTO", ncells = 5000)
ggsave(prostate_A_heatmap, filename =  paste(plot_dir, 
                                             "prostate_A_heatmap.png",
                                             sep= ""), height = 7, width =7 , units = "in")

# scatter plots
HTO_seurat <- as.data.frame(t(as.data.frame(hashtags2@assays$HTO@data)))


hto_pairs <- ggpairs(HTO_seurat)
hto_pairs<-pairs(HTO_seurat)
hto_pairs
ggsave(hto_pairs, filename =  paste(plot_dir,
                                    "_hto_seurat_pairs.png",
                                    sep= ""), height = 7, width =7 , units = "in")

# plot the doublet trends
doublet_trend <- ggplot(as.data.frame(hashtags2$HTO_classification), aes(hashtags2$HTO_classification)) +
  geom_bar() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(doublet_trend, filename =  paste(plot_dir,
                                        "_hto_seurat_doubet.png",
                                        sep= ""), height = 7, width =7 , units = "in")






#COMBINED

anchorstotal <- FindIntegrationAnchors(object.list = list(combined,combined2), dims = 1:20)
total_combined <- IntegrateData(anchorset = anchorstotal, dims = 1:20)

total_combined <- NormalizeData(object = total_combined)
total_combined <- FindVariableFeatures(object = total_combined, selection.method = "mean.var.plot")


# Run the standard workflow for visualization and clustering
total_combined <- ScaleData(object = total_combined, verbose = FALSE)
total_combined <- RunPCA(object = total_combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
total_combined <- RunUMAP(object = total_combined, reduction = "pca", dims = 1:20)
total_combined <- FindNeighbors(object = total_combined, reduction = "pca", dims = 1:20)
total_combined <- FindClusters(total_combined, resolution = 0.5)

p3 <- DimPlot(object = total_combined, reduction = "umap", group.by = "stim")
p4 <- DimPlot(object = total_combined, reduction = "umap", label = TRUE)
plot_grid(p3, p4)
