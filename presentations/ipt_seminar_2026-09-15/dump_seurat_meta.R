# Run ON THE CLUSTER (seurat_env), from the OCM directory. Takes seconds there;
# the 2.4 GB object does not load on a 16 GB laptop over the mount.
#
#   cd $BASE/OCM && conda activate seurat_env && Rscript ~/Postdoc/presentations/ipt_seminar_2026-09-15/dump_seurat_meta.R
#
# Writes OCM/seurat_metadata_umap.tsv: one row per nucleus (all four samples) with
# the QC metrics, cell type and UMAP coordinates. Copy it to
# ~/Downloads/IPT_Seminar_2026-09-15/data/ on the laptop and rerun make_figures.R
# + build_deck.py; slide 6 then uses all nuclei instead of the >= 30-read subset.
suppressPackageStartupMessages(library(Seurat))
h <- readRDS("heart_seurat_object_SCT.rds")
md <- h@meta.data
md$celltype <- as.character(Idents(h))
um <- Embeddings(h, "umap")
md$UMAP_1 <- um[rownames(md), 1]; md$UMAP_2 <- um[rownames(md), 2]
md$cell_barcode <- rownames(md)
write.table(md, "seurat_metadata_umap.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(md), "nuclei\n")
