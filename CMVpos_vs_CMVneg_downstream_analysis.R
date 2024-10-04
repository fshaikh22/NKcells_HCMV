#load packages:
library(dplyr)
library(parallel)
library(Seurat)
library(SeuratObject)
library(scater)
library(scran)
library(purrr)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(data.table)
library(stringr)
library(tibble)
library(ggrepel)

#import R.data files:

path = path <- "~/"

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
clusterpal <- brewer.pal(8, "Dark2")
clusterpal <- c(clusterpal[1:3],clusterpal[5:8])

# From GreenleafLab/ArchR
blueyellow <- c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60')
solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')


viriscale <- viridis(9)

#import seurat objects, re-normalize to regress out depth effects:
CMVpos_scRNA <- readRDS(paste0(path, "CMVpos_scRNA"))
CMVpos_scRNA <- NormalizeData(CMVpos_scRNA)
CMVpos_scRNA <- ScaleData(CMVpos_scRNA)

CMVneg_scRNA <- readRDS(paste0(path, "CMVneg_scRNA"))
CMVneg_scRNA <- NormalizeData(CMVneg_scRNA)
CMVneg_scRNA <- ScaleData(CMVneg_scRNA)

#visualize seurat objects with clustering
DimPlot(CMVpos_scRNA)
DimPlot(CMVneg_scRNA)

#visualize based on CD161 expression:
DimPlot(CMVpos_scRNA, group.by = "KLRB1", cols = clusterpal, shuffle = TRUE)
DimPlot(CMVneg_scRNA, group.by = "KLRB1", cols = clusterpal, shuffle = TRUE)

#Set idents to CD161+/- for downstream analysis
Idents(CMVpos_scRNA)
Idents(CMVpos_scRNA) <- CMVpos_scRNA$KLRB1

Idents(CMVneg_scRNA)
Idents(CMVneg_scRNA) <- CMVneg_scRNA$KLRB1

#change clustering via dims or res if required:
CMVpos_scRNA <- FindNeighbors(object = CMVpos_scRNA, dims = 1:10)
CMVpos_scRNA@graphs
CMVpos_scRNA <- FindClusters(object = CMVpos_scRNA, verbose = FALSE, graph.name = "integrated_snn", resolution = 0.5)
CMVpos_scRNA <- RunUMAP(object = CMVpos_scRNA, dims = 1:10)


CMVneg_scRNA <- FindNeighbors(object = CMVneg_scRNA, dims = 1:10)
CMVneg_scRNA@graphs
CMVneg_scRNA <- FindClusters(object = CMVneg_scRNA, verbose = FALSE, graph.name = "integrated_snn", resolution = 0.5)
CMVneg_scRNA <- RunUMAP(object = CMVneg_scRNA, dims = 1:10)


View(CMVpos_scRNA@meta.data)
View(CMVneg_scRNA@meta.data)

#further QC and filtering if required:
CMVpos_scRNA <- subset(CMVpos_scRNA, subset = nCount_RNA > 800 &
                         nFeature_RNA > 500 &
                         percent.mt < 5)

CMVneg_scRNA <- subset(CMVneg_scRNA, subset = nCount_RNA > 800 &
                         nFeature_RNA > 500 &
                         percent.mt < 5)


#genes for feature plots:
featureplot_genes = c("KLRC2", "B3GAT1", "FCER1G", "CD3E", "IL32", "SYK", "SH2D1B", "CD52")

#feature plots of cluster-wise gene expression:
FeaturePlot(CMVpos_scRNA, features = featureplot_genes, coord.fixed = T, pt.size = 0.2,
            cols =  colorscale, slot = "data", min.cutoff = "q1", max.cutoff = "q95")&NoLegend()&NoAxes()

FeaturePlot(CMVneg_scRNA, features = featureplot_genes, coord.fixed = T, pt.size = 0.2,
            cols =  colorscale, slot = "data", min.cutoff = "q1", max.cutoff = "q95")&NoLegend()&NoAxes()


# Extract the counts matrices:
counts_CMVP <- GetAssayData(CMVpos_scRNA, assay = "RNA", layer = "counts")
counts_CMVN <- GetAssayData(CMVneg_scRNA, assay = "RNA", layer = "counts")

# Get the metadata that includes the grouping variable
metadata_CMVP <- CMVpos_scRNA@meta.data
metadata_CMVN <- CMVneg_scRNA@meta.data


# Visualize DEGs between KLRB1+/-
# Then downsample to equal numbers to enable fair comparison
table(CMVpos_scRNA$KLRB1)
table(CMVneg_scRNA$KLRB1)

downsampled_CMVpos <- c(sample(names(which(CMVpos_scRNA$KLRB1 == "Pos")), size = 7500),
                        sample(names(which(CMVpos_scRNA$KLRB1 == "Neg")), size = 7500))

downsampled_CMVneg <- c(sample(names(which(CMVneg_scRNA$KLRB1 == "Pos")), size = 1800),
                        sample(names(which(CMVneg_scRNA$KLRB1 == "Neg")), size = 1800))


# Export downsampled cells to new object for reproducibility
CMVpos_down <- subset(CMVpos_scRNA, cells = downsampled_CMVpos)
CMVneg_down <- subset(CMVneg_scRNA, cells = downsampled_CMVneg)

#Run FindMarkers() for generating DEGs (default test.use parameter is Wilcoxon rank sum test):
markers_CD161_CMVpos <- FindMarkers(CMVpos_down, ident.1 = "Pos", ident.2 = "Neg", logfc.threshold = 0)
markers_CD161_CMVneg <- FindMarkers(CMVneg_down, ident.1 = "Pos", ident.2 = "Neg", logfc.threshold = 0)

#save markers as rds files:
saveRDS(markers_CD161_CMVpos, file = "CMVpos_markers")
saveRDS(markers_CD161_CMVneg, file = "CMVneg_markers")

#export to csv files:
write.csv(markers_CD161_CMVpos, "CMVpos_markers.csv", row.names=TRUE)
write.csv(markers_CD161_CMVneg, "markers_CD161_CMVneg.csv", row.names=TRUE)

# Add Diffexpr column
markers_CD161_CMVpos$diffexp <- "NO"
markers_CD161_CMVpos$diffexp[markers_CD161_CMVpos$p_val_adj < 0.05 & markers_CD161_CMVpos$avg_log2FC > 0.25] <- "UP"
markers_CD161_CMVpos$diffexp[markers_CD161_CMVpos$p_val_adj < 0.05 & markers_CD161_CMVpos$avg_log2FC < -0.25] <- "DOWN"

markers_CD161_CMVneg$diffexp <- "NO"
markers_CD161_CMVneg$diffexp[markers_CD161_CMVneg$p_val_adj < 0.05 & markers_CD161_CMVneg$avg_log2FC > 0.25] <- "UP"
markers_CD161_CMVneg$diffexp[markers_CD161_CMVneg$p_val_adj < 0.05 & markers_CD161_CMVneg$avg_log2FC < -0.25] <- "DOWN"


# Layer so that significant genes are on top
layer2_CMVpos <- markers_CD161_CMVpos[markers_CD161_CMVpos$diffexp%in%c("UP", "DOWN"),]
layer1_CMVpos <- markers_CD161_CMVpos[markers_CD161_CMVpos$diffexp%in%c("NO"),]

layer2_CMVneg <- markers_CD161_CMVneg[markers_CD161_CMVneg$diffexp%in%c("UP", "DOWN"),]
layer1_CMVneg <- markers_CD161_CMVneg[markers_CD161_CMVneg$diffexp%in%c("NO"),]

#Add a "gene" column
markers_CD161_CMVpos$gene <- rownames(markers_CD161_CMVpos)
markers_CD161_CMVneg$gene <- rownames(markers_CD161_CMVneg)


################# Plots for total NK cells in CMVpos#################

#remove all mitochondrial genes:
markers_CD161_CMVpos <- markers_CD161_CMVpos[!grepl('MT-', markers_CD161_CMVpos$gene),]

#remove all ribosomal genes:
markers_CD161_CMVpos <- markers_CD161_CMVpos[!grepl('^RP', markers_CD161_CMVpos$gene),]


# Create a new column "DElabel" to CMVpos, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
markers_CD161_CMVpos$DElabel_CMVpos <- ifelse(markers_CD161_CMVpos$gene %in% head(markers_CD161_CMVpos[order(markers_CD161_CMVpos$p_val_adj), "gene"], 30), markers_CD161_CMVpos$gene, NA)

#export DElabel column into a vector on top 30 DEGs:
top30_CMVpos <- as.vector(markers_CD161_CMVpos$DElabel_CMVpos)



#Volcano plot:
ggplot()+
  geom_point(data = layer1_CMVpos, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  geom_point(data = layer2_CMVpos, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = markers_CD161_CMVpos[markers_CD161_CMVpos$gene%in%c(top30_CMVpos)&!markers_CD161_CMVpos$diffexp=="NO",],
                  size = 5, box.padding = 0.5, max.overlaps = Inf, fontface = 'italic')+
  scale_y_continuous(limits = c(0, 350))+
  scale_x_continuous(limits = c(-2.65, 2.65))+
  xlab(label = "CD161+ vs CD161- in CMVpos patients (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()



#DotPlots:

#remove KLRB1 from the DEGs list for DotPlot:
top30_CMVpos_DotPlot <- top30_CMVpos[!(top30_CMVpos == "KLRB1")]

levels(CMVpos_scRNA) <- rev(levels(CMVpos_scRNA))
DotPlot(CMVpos_scRNA, features = top30_CMVpos_DotPlot,
        dot.min = 0.01, idents = c("Pos", "Neg"))+scale_color_gradientn(colours = solarExtra)+theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, face = 'italic'))


#rearrange genes in dotplot with upregulated ones in the beginning:
DotPlot(CMVpos_scRNA, features = c("FCER1G", "CLIC3", "AKR1C3", "CD160",  "CD7", "TXNIP", "PRF1", "TMIGD2", "JAK1", 
                                   "XCL2",  "ALOX5AP", "GSN", "CD247", "UBB", "GIMAP7", "PLAC8", "NKG7", "APMAP", "IL2RB", 
                                   "KLRC2", "IL32", "CD3E", "CD52","CCL5","GNLY", "CD3D", "GZMH", "TMSB4X", "PATL2"),
        dot.min = 0.01, dot.scale = 7, col.min = -6, col.max = 6) +scale_color_gradientn(colours = solarExtra)+theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, face = 'italic'))


##########Plots for total NK cells in CMVneg:#################

#remove all mitochondrial genes:
markers_CD161_CMVneg <- markers_CD161_CMVneg[!grepl('MT-', markers_CD161_CMVneg$gene),]

#remove all ribosomal genes:
markers_CD161_CMVneg <- markers_CD161_CMVneg[!grepl('^RP', markers_CD161_CMVneg$gene),]

# Create a new column "DElabel" to CMVneg, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)

markers_CD161_CMVneg$DElabel_CMVneg <- ifelse(markers_CD161_CMVneg$gene %in% head(markers_CD161_CMVneg[order(markers_CD161_CMVneg$p_val_adj), "gene"], 30), markers_CD161_CMVneg$gene, NA)

#export DElabel column into a vector on top 30 DEGs:
top30_CMVneg <- as.vector(markers_CD161_CMVneg$DElabel_CMVneg)

top10_CMVneg <- top30_CMVneg[1:10]

#volcano plot:

#Top 10 DEGs in CMVneg:

ggplot()+
  geom_point(data = layer1_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  geom_point(data = layer2_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = markers_CD161_CMVneg[markers_CD161_CMVneg$gene%in%c(top10_CMVneg)&!markers_CD161_CMVneg$diffexp=="NO",],
                  size = 5, box.padding = 0.7, max.overlaps = Inf, fontface = 'italic')+
  scale_y_continuous(limits = c(0, 350))+
  scale_x_continuous(limits = c(-2.65, 2.65))+
  xlab(label = "CD161+ vs CD161- in CMVneg individuals (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

#Top 30 DEGs in CMVneg:
ggplot()+
  geom_point(data = layer1_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  geom_point(data = layer2_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = markers_CD161_CMVneg[markers_CD161_CMVneg$gene%in%c(top30_CMVneg)&!markers_CD161_CMVneg$diffexp=="NO",],
                  size = 5, box.padding = 0.7, max.overlaps = Inf, fontface = 'italic')+
  scale_y_continuous(limits = c(0, 350))+
  scale_x_continuous(limits = c(-2.65, 2.65))+
  xlab(label = "CD161+ vs CD161- in CMVneg individuals (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()


#generate dotplot for clusters

top30_CMVneg_DotPlot <- top30_CMVneg[!(top30_CMVneg == "KLRB1")]
top10_CMVneg_DotPlot <- top30_CMVneg_DotPlot[1:10]

levels(CMVneg_scRNA) <- rev(levels(CMVneg_scRNA))
DotPlot(CMVneg_scRNA, features = top30_CMVneg_DotPlot,
        dot.min = 0.01, idents = c("Pos", "Neg"))+scale_color_gradientn(colours = solarExtra)+theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, face = 'italic'))

DotPlot(CMVneg_scRNA, features = top10_CMVneg_DotPlot,
        dot.min = 0.01, idents = c("Pos", "Neg"))+scale_color_gradientn(colours = solarExtra)+theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1, face = 'italic'))

