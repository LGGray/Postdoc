library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(scDEED)
options(future.globals.maxSize = 8 * 1024^3)

adult.file <- '9w/outs/cellbender_output_seurat.h5'
data.data <- Read10X_h5(filename = adult.file, use.names = TRUE)

adult <- CreateSeuratObject(counts = data.data)

# Preprocessing and Normalisation
adult <- PercentageFeatureSet(adult, pattern = "^MT-", col.name = "percent.mt")
adult <- subset(adult, subset = nCount_RNA > 500 & nCount_RNA < 20000 & nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
adult <- subset(adult, features = rownames(adult)[Matrix::rowSums(GetAssayData(adult, layer = "counts") > 0) > 10])
adult <- SCTransform(adult, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimension reduction
adult <- RunPCA(adult)

pdf('9w/elbow_plot.pdf')
ElbowPlot(adult)
dev.off()

adult <- RunUMAP(adult, dims = 1:10)

# Doublet detection
sweep.res.list_adult <- paramSweep(adult, PCs = 1:10, sct = TRUE)
sweep.stats_adult <- summarizeSweep(sweep.res.list_adult, GT = FALSE)
bcmvn_adult <- find.pK(sweep.stats_adult)

    bcmvn_adult <- bcmvn_adult %>%
    mutate(pK_num = as.numeric(as.character(pK)))

pdf('9w/pK.pdf')
ggplot(bcmvn_adult, aes(x = pK_num, y = BCmetric)) + 
geom_line() +
geom_point() +
labs(x = "pK", y = "BCmvn (BCmetric)", title = "DoubletFinder pK sweep") +
theme_minimal()
dev.off()

## pK = 0.05
annotations <- adult@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(adult@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
adult <- doubletFinder(adult, PCs = 1:10, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)

# Remove doublets
adult <- subset(adult, subset = DF.classifications_0.25_0.05_865 == "Singlet")

cell_index <- as.data.frame(colnames(adult))
rownames(cell_index) <- cell_index[,1]
write.table(cell_index, file = '9w/cell_index.txt', col.names = FALSE, row.names = TRUE, quote = FALSE)

saveRDS(adult, file = '9w/adult.RDS')

###################
aged.file <- '78w/outs/cellbender_output_seurat.h5'
aged.data <- Read10X_h5(filename = aged.file, use.names = TRUE)

aged <- CreateSeuratObject(counts = aged.data)

# Preprocessing and normalisation
aged <- PercentageFeatureSet(aged, pattern = "^MT-", col.name = "percent.mt")
aged <- subset(aged, subset = nCount_RNA > 500 & nCount_RNA < 20000 & nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
aged <- subset(aged, features = rownames(aged)[Matrix::rowSums(GetAssayData(aged, layer = "counts") > 0) > 10])
aged <- SCTransform(aged, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimension reduction
aged <- RunPCA(aged)

pdf('78w/elbow_plot.pdf')
ElbowPlot(aged)
dev.off()

aged <- RunUMAP(aged, dims = 1:10)

# Doublet detection
sweep.res.list_aged <- paramSweep(aged, PCs = 1:10, sct = TRUE)
sweep.stats_aged <- summarizeSweep(sweep.res.list_aged, GT = FALSE)
bcmvn_aged <- find.pK(sweep.stats_aged)

bcmvn_aged <- bcmvn_aged %>%
  mutate(pK_num = as.numeric(as.character(pK)))

pdf('78w/pK.pdf')
ggplot(bcmvn_aged, aes(x = pK_num, y = BCmetric)) + 
  geom_line() +
  geom_point() +
  labs(x = "pK", y = "BCmvn (BCmetric)", title = "DoubletFinder pK sweep") +
  theme_minimal()
dev.off()

## pK = 0.19
annotations <- aged@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(aged@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
aged <- doubletFinder(aged, PCs = 1:10, pN = 0.25, pK = 0.19, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)

# Remove doublets
aged <- subset(aged, subset = DF.classifications_0.25_0.19_804 == "Singlet")

cell_index <- as.data.frame(colnames(aged))
rownames(cell_index) <- cell_index[,1]
write.table(cell_index, file = '78w/cell_index.txt', col.names = FALSE, row.names = TRUE, quote = FALSE)

saveRDS(aged, file = '78w/aged.RDS')

# Integration - https://satijalab.org/seurat/articles/seurat5_integration
adult <- readRDS('9w/adult.RDS')
aged <- readRDS('78w/aged.RDS')

adult$age <- 'adult'
aged$age  <- 'aged'

merged <- merge(adult, aged, add.cell.ids = c("adult","aged"), project = "combined")
DefaultAssay(merged) <- "RNA"
merged <- SCTransform(merged, verbose = FALSE)

merged <- RunPCA(merged)
merged[["SCT"]] <- split(merged[["SCT"]], f = merged$age)
merged <- IntegrateLayers(object = merged, method = RPCAIntegration, assay='SCT')
merged <- JoinLayers(merged)


pdf('elbow_plot.pdf')
ElbowPlot(merged, reduction = "pca", ndims=30)
dev.off()


merged <- FindNeighbors(merged, reduction = "integrated.dr", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.1, cluster.name = "seurat_clusters")
Idents(merged) <- 'seurat_clusters'

merged <- RunUMAP(merged, reduction = "integrated.dr", dims = 1:30)

pdf('adult_aged_umap.pdf')
DimPlot(merged, reduction = "umap", label = TRUE)
dev.off()


# ### scDEED for integrated data - https://github.com/JSB-UCLA/scDEED
# K=
# pre_embedding = 'integrated.dr'
# rpca_space = Embeddings(merged, pre_embedding)
# permuted_rpca_space = rpca_space
# ## We will now permute the integrated RPCA space.
# ##This is the same as the Permuted function, but we have swapped rows/columns because embeddings are cell x features while the scale.data slot is features x cells
# set.seed(1000)
# for (i in 1:dim(permuted_rpca_space)[2])
# {
#   row = pracma::randperm(dim(permuted_rpca_space)[1])
#   permuted_rpca_space[,i]=rpca_space[row,i]

# }
# ##Now, we will add our permuted RPCA space as a reduced dimension object. 
# data.permuted = merged
# data.permuted[[pre_embedding]] <- CreateDimReducObject(embeddings =permuted_rpca_space , key = "integratedrpca_", assay = DefaultAssay(data.permuted))

# ##Now we can call scDEED and 1) provide our own permuted object 2) specify that the pre_embedding space is integrated.cca
# result = scDEED(merged, K = K, pre_embedding = pre_embedding, permuted = data.permuted, reduction.method = 'umap')
# result$num_dubious

# # Minimum number of dubious cells at n_neighbors, min_dist
# m <- 0.1
# n <- 5

# merged <- RunUMAP(merged, dims = 1:5, min.dist = m, n.neighbors = n, seed.use = 100)

# ### Find neighbours and clusters
# merged <- FindNeighbors(merged, reduction = "pca", dims = 1:K)
# merged <- FindClusters(merged, resolution = 0.1, cluster.name = "seurat_clusters")



pdf('adult_aged_umap.pdf')
DimPlot(merged, reduction = "umap", label = TRUE)
dev.off()

saveRDS(merged, file = 'adult_aged_merged.RDS')

# Load in marker genes
marker_genes <- read.table('marker.genes.txt', header = FALSE)

# Plot marker genes
pdf('adult_aged_marker_genes.pdf', width = 10, height = 8)
DotPlot(merged, features = marker_genes$V1) + RotatedAxis() +
theme(axis.text.x = element_text(angle = 90, vjust=1, size=8))
dev.off()

pdf('adult_aged_marker_genes_heatmap.pdf', width = 10, height = 8)
DoHeatmap(merged, features = marker_genes$V1, assay = "SCT") + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


merged$celltype <- plyr::mapvalues(Idents(merged),
from = 0:11,
to = c("FB","ventCM1","ventCM2","EC2","EC1","SMC",
"Macro","atrialCM","LEC","Epicard","Adipo","Glial"))

Idents(merged) <- 'celltype'

pdf('adult_aged_umap_celltype.pdf')
DimPlot(merged, reduction = "umap")
dev.off()

saveRDS(merged, file = 'adult_aged_merged.RDS')

pdf('batch_umap.pdf')
DimPlot(merged, reduction = "umap", group.by = "age")
dev.off()

cell_ids <- data.frame(cell_id = colnames(merged), cell_type = merged$celltype, age = merged$age)
cell_ids$cell_id <- gsub('adult_|aged_', '', cell_ids$cell_id)

adult <- cell_ids[cell_ids$age == 'adult', 1:2]
aged <- cell_ids[cell_ids$age == 'aged', 1:2]

write.table(adult, file = '9w/cell_index.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(aged, file = '78w/cell_index.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
