library(dplyr)
library(parallel)
library(Seurat)
library(scater)
library(scran)
library(purrr)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(data.table)
library(stringr)

path <- "~/"

# Sample names
my_samples <- c("CMVpos2", "CMVpos3", "CMVpos4", "CMVneg1", "CMVneg2")

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
clusterpal <- brewer.pal(8, "Dark2")
clusterpal <- c(clusterpal[1:3],clusterpal[5:8])

viriscale <- viridis(9)

##### Load the datasets #####

CMVpos2.data <- Read10X(data.dir = paste0(path, "CMVpos2"), strip.suffix = T)
CMVpos3.data <- Read10X(data.dir = paste0(path, "CMVpos3"), strip.suffix = T)
CMVpos4.data <- Read10X(data.dir = paste0(path, "CMVpos4"), strip.suffix = T)
CMVneg1.data <- Read10X(data.dir = paste0(path, "CMVneg1"), strip.suffix = T)
CMVneg2.data <- Read10X(data.dir = paste0(path, "CMVneg2"), strip.suffix = T)

# Initialize the Seurat objects with the raw (non-normalized data).
# Keep all features expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected features
CMVpos2 <- CreateSeuratObject(counts = CMVpos2.data, min.cells = 3, min.features = 200, project = "CMVpos2")
CMVpos3 <- CreateSeuratObject(counts = CMVpos3.data, min.cells = 3, min.features = 200, project = "CMVpos3")
CMVpos4 <- CreateSeuratObject(counts = CMVpos4.data, min.cells = 3, min.features = 200, project = "CMVpos4")
CMVneg1 <- CreateSeuratObject(counts = CMVneg1.data, min.cells = 3, min.features = 200, project = "CMVneg1")
CMVneg2 <- CreateSeuratObject(counts = CMVneg2.data, min.cells = 3, min.features = 200, project = "CMVneg2")

CMVpos2
CMVpos3
CMVpos4
CMVneg1
CMVneg2

##### QC #####
# store mitochondrial percentage in object meta data
CMVpos2 <- PercentageFeatureSet(CMVpos2, pattern = "^MT-", col.name = "percent.mt")
CMVpos3 <- PercentageFeatureSet(CMVpos3, pattern = "^MT-", col.name = "percent.mt")
CMVpos4 <- PercentageFeatureSet(CMVpos4, pattern = "^MT-", col.name = "percent.mt")
CMVneg1 <- PercentageFeatureSet(CMVneg1, pattern = "^MT-", col.name = "percent.mt")
CMVneg2 <- PercentageFeatureSet(CMVneg2, pattern = "^MT-", col.name = "percent.mt")

# Filter based on QC metrics
CMVpos2 <- subset(CMVpos2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 6)
CMVpos3 <- subset(CMVpos3, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVpos4 <- subset(CMVpos4, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg1 <- subset(CMVneg1, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg2 <- subset(CMVneg2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)


# Visualize markers of contaminating cells

VlnPlot(object = CMVpos2, features = c("HBA1", "HBA2", "IGJ"), slot = "counts")
VlnPlot(object = CMVpos3, features = c("HBA1", "HBA2", "IGJ"), slot = "counts")
VlnPlot(object = CMVpos4, features = c("HBA1", "HBA2", "IGJ"), slot = "counts")
VlnPlot(object = CMVneg1, features = c("HBA1", "HBA2", "IGJ"), slot = "counts")
VlnPlot(object = CMVneg2, features = c("HBA1", "HBA2", "IGJ"), slot = "counts")

#Normalize objects in the list

CMVpos2 <- NormalizeData(CMVpos2)
CMVpos3 <- NormalizeData(CMVpos3)
CMVpos4 <- NormalizeData(CMVpos4)
CMVneg1 <- NormalizeData(CMVneg1)
CMVneg2 <- NormalizeData(CMVneg2)


#Filter out RBCs from each Seurat object
CMVpos2 <- subset(CMVpos2, subset = HBA1 == 0 & HBA2 == 0 & IGJ == 0)
CMVpos3 <- subset(CMVpos3, subset = HBA1 == 0 & HBA2 == 0 & IGJ == 0)
CMVpos4 <- subset(CMVpos4, subset = HBA1 == 0 & HBA2 == 0 & IGJ == 0)
CMVneg1 <- subset(CMVneg1, subset = HBA1 == 0 & HBA2 == 0 & IGJ == 0)
CMVneg2 <- subset(CMVneg2, subset = HBA1 == 0 & HBA2 == 0 & IGJ == 0)




##### Variable Features, Scaling #####

# Find variable features

CMVpos2 <- FindVariableFeatures(CMVpos2, selection.method = "vst", nfeatures = 2000)
CMVpos3 <- FindVariableFeatures(CMVpos3, selection.method = "vst", nfeatures = 2000)
CMVpos4 <- FindVariableFeatures(CMVpos4, selection.method = "vst", nfeatures = 2000)
CMVneg1 <- FindVariableFeatures(CMVneg1, selection.method = "vst", nfeatures = 2000)
CMVneg2 <- FindVariableFeatures(CMVneg2, selection.method = "vst", nfeatures = 2000)

#ScaleData

CMVpos2 <- ScaleData(CMVpos2, features = rownames(CMVpos2))
CMVpos3 <- ScaleData(CMVpos3, features = rownames(CMVpos3))
CMVpos4 <- ScaleData(CMVpos4, features = rownames(CMVpos4))
CMVneg1 <- ScaleData(CMVneg1, features = rownames(CMVneg1))
CMVneg2 <- ScaleData(CMVneg2, features = rownames(CMVneg2))

#RunPCA

CMVpos2 <- RunPCA(CMVpos2, verbose = FALSE)
CMVpos3 <- RunPCA(CMVpos3, verbose = FALSE)
CMVpos4 <- RunPCA(CMVpos4, verbose = FALSE)
CMVneg1 <- RunPCA(CMVneg1, verbose = FALSE)
CMVneg2 <- RunPCA(CMVneg2, verbose = FALSE)






# Inspect PCs in Elbowplot and DimHeatmaps
pdf(file = "PCA_results.pdf", title = "PCA Results", paper = "a4")

DimHeatmap(CMVpos2, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVpos2, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVpos2, dims = 19:27, cells = 500, balanced = TRUE, nfeatures = 22)

DimHeatmap(CMVpos3, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVpos3, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVpos3, dims = 19:27, cells = 500, balanced = TRUE, nfeatures = 22)

DimHeatmap(CMVpos4, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVpos4, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVpos4, dims = 19:27, cells = 500, balanced = TRUE, nfeatures = 22)

DimHeatmap(CMVneg1, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVneg1, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVneg1, dims = 19:27, cells = 500, balanced = TRUE, nfeatures = 22)

DimHeatmap(CMVneg2, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVneg2, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
DimHeatmap(CMVneg2, dims = 19:27, cells = 500, balanced = TRUE, nfeatures = 22)

ElbowPlot(CMVpos2)
ElbowPlot(CMVpos3)
ElbowPlot(CMVpos4)
ElbowPlot(CMVneg1)
ElbowPlot(CMVneg2)

dev.off()




# RunUMAP
CMVpos2 <- RunUMAP(CMVpos2, dims = 1:13, verbose = FALSE)
CMVpos3 <- RunUMAP(CMVpos3, dims = 1:14, verbose = FALSE)
CMVpos4 <- RunUMAP(CMVpos4, dims = 1:13, verbose = FALSE)
CMVneg1 <- RunUMAP(CMVneg1, dims = 1:10, verbose = FALSE)
CMVneg2 <- RunUMAP(CMVneg2, dims = 1:10, verbose = FALSE)

##### Clustering #####
# High resolution in some donors (i.e. overclustering)
# to identify contaminating ILCs (are in proximity to bright)
CMVpos2 <- FindNeighbors(CMVpos2, dims = 1:13, verbose = FALSE)
CMVpos2 <- FindClusters(CMVpos2, verbose = FALSE, resolution = 2, n.start = 100)

CMVpos3 <- FindNeighbors(CMVpos3, dims = 1:14, verbose = FALSE)
CMVpos3 <- FindClusters(CMVpos3, verbose = FALSE, resolution = 1, n.start = 100)

CMVpos4 <- FindNeighbors(CMVpos4, dims = 1:13, verbose = FALSE)
CMVpos4 <- FindClusters(CMVpos4, verbose = FALSE, resolution = 0.6, n.start = 100)

CMVneg1 <- FindNeighbors(CMVneg1, dims = 1:10, verbose = FALSE)
CMVneg1 <- FindClusters(CMVneg1, verbose = FALSE, resolution = 0.6, n.start = 100)

CMVneg2 <- FindNeighbors(CMVneg2, dims = 1:10, verbose = FALSE)
CMVneg2 <- FindClusters(CMVneg2, verbose = FALSE, resolution = 0.4, n.start = 100)


# Plot UMAPs with clustering and expression of ILC/NK markers

pdf(file = "UMAP_CMVpos2_ILCs.pdf", paper = "a4")
DimPlot(CMVpos2, label = T)
FeaturePlot(CMVpos2, features = c("CD56", "IL7R", "GATA3", "IL2RA", "CD40LG"), cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE)

pdf(file = "UMAP_CMVpos3_ILCs.pdf", paper = "a4")
DimPlot(CMVpos3, label = T)
FeaturePlot(CMVpos3, features = c("CD56", "IL7R", "GATA3", "IL2RA", "CD40LG"), cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE)

pdf(file = "UMAP_CMVpos4_ILCs.pdf", paper = "a4")
DimPlot(CMVpos4, label = T)
FeaturePlot(CMVpos4, features = c("CD56", "IL7R", "GATA3", "IL2RA", "CD40LG"), cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE)

pdf(file = "UMAP_CMVneg1_ILCs.pdf", paper = "a4")
DimPlot(CMVneg1, label = T)
FeaturePlot(CMVneg1, features = c("CD56", "IL7R", "GATA3", "IL2RA", "CD40LG"), cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE)

pdf(file = "UMAP_CMVneg2_ILCs.pdf", paper = "a4")
DimPlot(CMVneg2, label = T)
FeaturePlot(CMVneg2, features = c("CD56", "IL7R", "GATA3", "IL2RA", "CD40LG"), cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE)

dev.off()





##### Removal of contaminating cells (mostly ILC2s: High IL7R, GATA3, IL2RA, CD40LG, negative for CD56) #####
CMVpos2 <- subset(CMVpos2, idents = "13", invert = TRUE) # Cluster 13 are ILC2
CMVpos3 <- subset(CMVpos3, idents = "6", invert = TRUE) # Cluster 6 are ILC2
CMVpos4 <- subset(CMVpos4, idents = c("1", "4", "5", "6"), invert = TRUE) # Cluster 1, 4, 5, 6 are ILC2
CMVneg1 <- subset(CMVneg1, idents = "0", invert = TRUE) # Cluster 0 are ILC2
CMVneg2 <- subset(CMVneg2, idents = "3", invert = TRUE) # Cluster 3 and possibly 1 are ILC2


##### Re-run Dimred and Clustering #####

CMVpos2 <- RunPCA(CMVpos2, verbose = FALSE)
CMVpos3 <- RunPCA(CMVpos3, verbose = FALSE)
CMVpos4 <- RunPCA(CMVpos4, verbose = FALSE)
CMVneg1 <- RunPCA(CMVneg1, verbose = FALSE)
CMVneg2 <- RunPCA(CMVneg2, verbose = FALSE)

# RunUMAP
CMVpos2 <- RunUMAP(CMVpos2, dims = 1:11, verbose = FALSE, n.neighbors = 20)
CMVpos3 <- RunUMAP(CMVpos3, dims = 1:14, verbose = FALSE, n.neighbors = 20)
CMVpos4 <- RunUMAP(CMVpos4, dims = 1:15, verbose = FALSE, n.neighbors = 20)
CMVneg1 <- RunUMAP(CMVneg1, dims = 1:10, verbose = FALSE, n.neighbors = 20)
CMVneg2 <- RunUMAP(CMVneg2, dims = 1:13, verbose = FALSE, n.neighbors = 20)


# Add donor annotation
CMVpos2$donor <- "CMVpos2"
CMVpos3$donor <- "CMVpos3"
CMVpos4$donor <- "CMVpos4"
CMVneg1$donor <- "CMVneg1"
CMVneg2$donor <- "CMVneg2"


# Assign serostatus
CMVpos2$serostatus <- "CMVpos"
CMVpos3$serostatus <- "CMVpos"
CMVpos4$serostatus <- "CMVpos"
CMVneg1$serostatus <- "CMVneg"
CMVneg2$serostatus <- "CMVneg"


##### Integration of datasets #####

# Find overlapping variable features of CMVpos
shared_variable_features_CMVpos <- Reduce(intersect, list(VariableFeatures(CMVpos2),
                                                          VariableFeatures(CMVpos3),
                                                          VariableFeatures(CMVpos4)
))

shared_features <- Reduce(intersect, list(rownames(CMVpos2),
                                          rownames(CMVpos3),
                                          rownames(CMVpos4),
                                          rownames(CMVneg1),
                                          rownames(CMVneg2)
))

all_features <- Reduce(unique, list(rownames(CMVpos2),
                                    rownames(CMVpos3),
                                    rownames(CMVpos4),
                                    rownames(CMVneg1),
                                    rownames(CMVneg2)
))


##############Integrate CMVpos###########

integration_anchors_CMVpos <- FindIntegrationAnchors(c(CMVpos2, CMVpos3, CMVpos4),
                                                     dims = 1:10, anchor.features = shared_variable_features_CMVpos)
CMVpos <- IntegrateData(integration_anchors_CMVpos, dims = 1:10, features.to.integrate = all_features)
DefaultAssay(CMVpos) <- "integrated"


#Join Layers
CMVpos[["RNA"]] <- JoinLayers(CMVpos[["RNA"]])

#Add CD161+/- labels based on gene expression instead of hashtags
grep('KLRB1', rownames(CMVpos@assays$RNA$counts))
length(which(CMVpos@assays$RNA$counts[rownames(CMVpos@assays$RNA$counts)[9330], ] != 0)) 
CMVpos@meta.data$KLRB1 <- 'Neg'
CMVpos@meta.data$KLRB1[which(CMVpos@assays$RNA$counts[rownames(CMVpos@assays$RNA$counts)[9330], ] != 0)] <- 'Pos'



# Process integrated dataset
CMVpos <- ScaleData(CMVpos, vars.to.regress = "nCount_RNA")

CMVpos <- RunPCA(CMVpos)

CMVpos <- RunUMAP(CMVpos, dims = 1:10, n.neighbors = 20, min.dist = 0.2)
CMVpos <- FindNeighbors(CMVpos, dims = 1:10)
CMVpos <- FindClusters(CMVpos, resolution = 0.3)


#Visualize clustering in objects:

DimPlot(CMVpos)

DimPlot(CMVpos, group.by = "KLRB1", cols = clusterpal, shuffle = TRUE)



# Integration of CMVneg (use the same anchor features as for CMVpos to enable merging afterwards)

integration_anchors_CMVneg <- FindIntegrationAnchors(c(CMVneg1, CMVneg2),
                                                     dims = 1:10, anchor.features = shared_variable_features_CMVpos)
CMVneg <- IntegrateData(integration_anchors_CMVneg, dims = 1:10, features.to.integrate = all_features)

CMVneg[["RNA"]] <- JoinLayers(CMVneg[["RNA"]])

#Add CD161 expression to object:
grep('KLRB1', rownames(CMVneg@assays$RNA$counts))
length(which(CMVneg@assays$RNA$counts[rownames(CMVneg@assays$RNA$counts)[8943], ] != 0)) 
CMVneg@meta.data$KLRB1 <- 'Neg'
CMVneg@meta.data$KLRB1[which(CMVneg@assays$RNA$counts[rownames(CMVneg@assays$RNA$counts)[8943], ] != 0)] <- 'Pos'


# Process integrated dataset
DefaultAssay(CMVneg) <- "integrated"
CMVneg <- ScaleData(CMVneg, vars.to.regress = "nCount_RNA")
CMVneg <- RunPCA(CMVneg)
CMVneg <- RunUMAP(CMVneg, dims = 1:10)
DimPlot(CMVneg)

CMVneg <- FindNeighbors(CMVneg, dims = 1:10)
CMVneg <- FindClusters(CMVneg, resolution = 0.2)


#visualize CD161+/- CMVneg objects before labelling clusters:
DimPlot(CMVneg)

DimPlot(CMVneg, group.by = "KLRB1", cols = clusterpal, shuffle = TRUE)

# Normalize RNA assay for full dataset for downstream analysis
DefaultAssay(CMVpos) <- "RNA"
DefaultAssay(CMVneg) <- "RNA"

CMVpos <- NormalizeData(CMVpos)
CMVneg <- NormalizeData(CMVneg)

CMVpos <- ScaleData(CMVpos)
CMVneg <- ScaleData(CMVneg)

# Export integrated datasets
saveRDS(CMVpos, file = "CMVpos_scRNA")
saveRDS(CMVneg, file = "CMVneg_scRNA")
