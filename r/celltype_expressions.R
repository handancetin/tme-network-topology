# Dependecies
library(dplyr)
library(readxl)
library(Seurat) #Seurat_5.0.1   
library(harmony)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

########## Set Environment ##########
set.seed(123)
setwd("/Users/hcetin/Documents/CSBL/github-projects/tme-cfm-gems/r/scripts")

########## Seurat analysis ##########
## Load data from Qi et. al. (Nat Commun. 2022)
metadata <- read_excel("../input/41467_2022_29366_MOESM7_ESM.xlsx")
counts <- readr::read_rds("../input/41467_2022_29366_MOESM6_ESM.gz")
seuobj <- CreateSeuratObject(counts = counts , project = "scrna-crc-r", meta.data = metadata)

## Investigate & check if qc is performed
Idents(seuobj) <- "MainTypes"
#VlnPlot(seuobj, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), pt.size = 0)  + theme(legend.position = 'none')
#FeatureScatter(seuobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Get maintype identities of interest (i.e. remove T cells, B cells, endothelials etc.) 

seuobj <- subset(seuobj, idents = c("Epithelial", "Stroma", "Myeloid"))
Idents(seuobj) <- "Cell Types"
VlnPlot(seuobj, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), pt.size = 0)  + theme(legend.position = 'none')
table(seuobj@meta.data$`Cell Types`, seuobj@meta.data$Tissues)

# Filter by cell types and investigate
levels(Idents(seuobj))
selectedTypes <- c("Malignant cells",
                  "FGFR2+ fibroblasts","FAP+ fibroblasts", "CD73+ fibroblasts", 
                  "DES+ myofibroblasts",  "MFAP5+ myofibroblasts", "ICAM1+ telocytes", "ICAM1- telocytes", 
                  "THBS1+ Macrophage", "VCAN+ Monocyte", "MARCO+ Macrophage")
seuobj.filtered <- subset(seuobj, idents = selectedTypes)
FeatureScatter(seuobj.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #0.92
table(seuobj.filtered@meta.data$`Cell Types`, seuobj.filtered@meta.data$Tissues) 
# Note: MARCO+ Macrophage from Normal tissue, and MFAP5+ myofibroblasts from Tumor tissue 
# will be removed later to obtain UMAP figure as in Qi et. al.

# Run seurat pipeline
seuobj.filtered <- NormalizeData(seuobj.filtered) 
seuobj.filtered <- FindVariableFeatures(seuobj.filtered, mean.function = ExpMean, dispersion.function = LogVMR)
seuobj.filtered <- ScaleData(seuobj.filtered)
seuobj.filtered <- RunPCA(seuobj.filtered, npcs = 30, verbose = FALSE) 
#ElbowPlot(object = seuobj.filtered)
seuobj.filtered <- RunUMAP(seuobj.filtered, reduction = "pca", dims = 1:10,  n.neighbors = 30L,   n.components = 2L) 

# Visualize with UMAP
selectedTypesColors <- c("Malignant cells" = "#4D4D4D",
                         "FGFR2+ fibroblasts" = "#1B76DE",
                         "FAP+ fibroblasts" = "#DD1C1A", 
                         "CD73+ fibroblasts" = "#7CB466", 
                         "DES+ myofibroblasts" = "#E8A621",  
                         "MFAP5+ myofibroblasts" = "#F8DD6C", 
                         "ICAM1+ telocytes" = "#086788", 
                         "ICAM1- telocytes" = "#42AA95", 
                         "THBS1+ Macrophage" = "#AB82FF", 
                         "VCAN+ Monocyte" =  "#DB7092", 
                         "MARCO+ Macrophage" = "#FFB5C5")
umap1 <- DimPlot(seuobj.filtered, reduction = "umap", label = TRUE, cols = selectedTypesColors) + NoLegend()
umap2 <- DimPlot(seuobj.filtered, reduction = "umap", split.by = "Tissues", cols = selectedTypesColors)
umap1 + umap2 + plot_layout(widths = c(1, 2))

# Extract UMAP coordinates
umap_coords <- Embeddings(seuobj.filtered, reduction = "umap")
metadata <- seuobj.filtered@meta.data
umap_data <- cbind(umap_coords, metadata)
write.csv(umap_data, "../output/umap_pos_data.csv", row.names = FALSE)


########## Generate pseudobulk matrices ##########
metgenes <- read.delim2("/Users/hcetin/Documents/CSBL/github-projects/third-party/Human-GEM-1.18.0/model/genes.tsv")
metgenes <- data.frame( genes = metgenes$genes,  geneSymbols = metgenes$geneSymbols) 

pseudobulk <- PseudobulkExpression(seuobj.filtered, 
                                   features = metgenes$geneSymbols,  
                                   method = "average",
                                   normalization.method = "LogNormalize",  
                                   group.by = c("ident", "Tissues"),
                                   verbose = TRUE,
                                   return.seurat = TRUE)
pseudobulk <- GetAssayData(pseudobulk)
pseudobulk <- pseudobulk * 10^6 / sum(pseudobulk)
boxplot(pseudobulk,range=0)
head(pseudobulk)

# Change geneSymbols back to ENSGs to match human-gem
row.names(pseudobulk) <- metgenes$genes[match(row.names(pseudobulk), metgenes$geneSymbols)]

# Save the matrix for downstream
write.csv(pseudobulk,"../output/averaged_lognorm_expressions_cpm.csv")


########## Perform differential expression analysis by clusters ##########

# # Get the list of metabolic genes
seuobj.metfeats <- subset(seuobj.filtered, features = metgenes$geneSymbols)

Idents(seuobj.metfeats) <- "Tissues"

# Eliminate cells  n<50
seuobj.metfeats.normal <- subset(seuobj.metfeats, idents = "N")
Idents(seuobj.metfeats.normal) <- "Cell Types"
Idents(seuobj.metfeats.normal) <- factor(Idents(seuobj.metfeats.normal), levels=rev(selectedTypes))
seuobj.metfeats.normal <- subset(seuobj.metfeats.normal, idents = "MARCO+ Macrophage", invert = TRUE)

seuobj.metfeats.tumor <- subset(seuobj.metfeats, idents = "T")
Idents(seuobj.metfeats.tumor) <- "Cell Types"
Idents(seuobj.metfeats.tumor) <- factor(Idents(seuobj.metfeats.tumor), levels=rev(selectedTypes))
seuobj.metfeats.tumor <- subset(seuobj.metfeats.tumor, idents = "MFAP5+ myofibroblasts", invert = TRUE)

# Find what genes derive clusters   : NORMAL
markers.met.normal.wilcox <- FindAllMarkers(seuobj.metfeats.normal, test.use = "wilcox") 
markers.met.normal.top <- markers.met.normal.wilcox %>% filter(pct.1 > 0.5) %>%  filter(pct.2 < 0.2) %>% group_by(cluster) %>% top_n(5, avg_log2FC)
markers.met.normal.top

# Find what genes derive clusters   : TUMOR
markers.met.tumor.wilcox <- FindAllMarkers(seuobj.metfeats.tumor, test.use = "wilcox") 
markers.met.tumor.top <- markers.met.tumor.wilcox %>% filter(pct.1 > 0.5) %>%  filter(pct.2 < 0.2) %>% group_by(cluster) %>% top_n(5, avg_log2FC)
markers.met.tumor.top

p1 <- DotPlot(seuobj.metfeats.normal, features = unique(markers.met.normal.top$gene)) + RotatedAxis() +  xlab('Normal')+ ylab('')
p2 <- DotPlot(seuobj.metfeats.tumor, features = unique(markers.met.tumor.top$gene)) + RotatedAxis()  +  xlab('Tumor') + ylab('')
p1 / p2 

# Save data for python
write.csv(p1$data, file = "../output/dotplot_data_metmarkers_normal.csv", row.names = FALSE)
write.csv(p2$data, file = "../output/dotplot_data_metmarkers_tumor.csv", row.names = FALSE)
