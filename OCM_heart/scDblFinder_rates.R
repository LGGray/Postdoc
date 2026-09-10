# ---------------------------------------------------------------------------
# scDblFinder doublet rates on the already-built merged heart object.
#
# Run from the OCM/ directory:  Rscript $POSTDOC_ROOT/OCM_heart/scDblFinder_rates.R
#
# This is the read-only companion to the scDblFinder block that now lives in
# Seurat_preprocessing.R. It reports the per-sample doublet rate from the
# saved object without re-running the pipeline, so the existing rds and the
# sinto_ID.txt barcode lists are left untouched. Re-run
# Seurat_preprocessing.R when you actually want the doublets dropped and the
# barcode lists regenerated.
#
# Writes: Allelic_ratio_results/scDblFinder_rates_per_sample.txt
#         Allelic_ratio_results/scDblFinder_per_cell.txt
#         Allelic_ratio_results/scDblFinder_QC.pdf
# ---------------------------------------------------------------------------
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)
library(ggplot2)

options(future.globals.maxSize = 8 * 1024^3)

DBL_THREADS <- as.integer(Sys.getenv("DBL_THREADS", "4"))
OUT_DIR     <- "Allelic_ratio_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

heart <- readRDS("heart_seurat_object_SCT.rds")
if (!"celltype" %in% colnames(heart@meta.data)) heart$celltype <- Idents(heart)

# scDblFinder wants raw counts, not SCT residuals.
DefaultAssay(heart) <- "RNA"
heart[["RNA"]] <- JoinLayers(heart[["RNA"]])

dbl_sce <- SingleCellExperiment(
  assays = list(counts = GetAssayData(heart, assay = "RNA", layer = "counts"))
)
dbl_sce$sample <- heart$sample

# samples= keeps artificial doublets within a sample; clusters=TRUE builds them
# between cell types rather than at random, which is the right choice in heart.
set.seed(42)
dbl_sce <- scDblFinder(
  dbl_sce,
  samples  = "sample",
  clusters = TRUE,
  BPPARAM  = MulticoreParam(DBL_THREADS)
)

stopifnot(identical(colnames(dbl_sce), colnames(heart)))
heart$scDblFinder.class <- factor(as.character(dbl_sce$scDblFinder.class),
                                  levels = c("singlet", "doublet"))
heart$scDblFinder.score <- dbl_sce$scDblFinder.score

# Per-sample rate
dbl_tab <- table(heart$sample, heart$scDblFinder.class)
dbl_rates <- data.frame(
  sample      = rownames(dbl_tab),
  n_cells     = as.integer(rowSums(dbl_tab)),
  singlet     = as.integer(dbl_tab[, "singlet"]),
  doublet     = as.integer(dbl_tab[, "doublet"]),
  pct_doublet = round(100 * dbl_tab[, "doublet"] / rowSums(dbl_tab), 2),
  row.names   = NULL
)
cat("\n--- scDblFinder doublet rate per sample ---\n")
print(dbl_rates)
cat(sprintf("\nOverall: %d / %d cells = %.2f%% doublets\n",
            sum(dbl_rates$doublet), sum(dbl_rates$n_cells),
            100 * sum(dbl_rates$doublet) / sum(dbl_rates$n_cells)))

# Which cell types the doublets land in - a doublet call concentrated in one
# annotated cluster usually means that cluster IS the doublet population.
ct_tab <- table(heart$celltype, heart$scDblFinder.class)
ct_rates <- data.frame(
  celltype    = rownames(ct_tab),
  n_cells     = as.integer(rowSums(ct_tab)),
  doublet     = as.integer(ct_tab[, "doublet"]),
  pct_doublet = round(100 * ct_tab[, "doublet"] / rowSums(ct_tab), 2),
  row.names   = NULL
)
ct_rates <- ct_rates[order(-ct_rates$pct_doublet), ]
cat("\n--- doublet rate per annotated cell type ---\n")
print(ct_rates)

write.table(dbl_rates, file.path(OUT_DIR, "scDblFinder_rates_per_sample.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ct_rates, file.path(OUT_DIR, "scDblFinder_rates_per_celltype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(
    cell     = colnames(heart),
    sample   = heart$sample,
    celltype = heart$celltype,
    class    = heart$scDblFinder.class,
    score    = heart$scDblFinder.score
  ),
  file.path(OUT_DIR, "scDblFinder_per_cell.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pdf(file.path(OUT_DIR, "scDblFinder_QC.pdf"), width = 11, height = 5)
print(DimPlot(heart, reduction = "umap", group.by = "scDblFinder.class",
              order = "doublet", raster = FALSE))
print(FeaturePlot(heart, features = "scDblFinder.score", raster = FALSE))
print(VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA"), group.by = "sample",
              split.by = "scDblFinder.class", pt.size = 0, ncol = 2))
print(
  ggplot(ct_rates, aes(x = reorder(celltype, pct_doublet), y = pct_doublet)) +
    geom_col(fill = "grey30") + coord_flip() +
    labs(x = NULL, y = "% doublet", title = "scDblFinder rate per cell type") +
    theme_minimal()
)
dev.off()
