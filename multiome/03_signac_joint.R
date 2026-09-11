# ---------------------------------------------------------------------------
# 03 - Signac: common peak set, TRUE cross-sample joint WNN, coverage tracks.
#
# What this adds over 02_cluster_umap.R, and why it needed Signac:
#
#   02 could only do WNN PER SAMPLE. cellranger called peaks independently for
#   each sample (52,231 at 9w, 71,694 at 78w) and computed each LSI on its own
#   peak set, so 9w component 3 and 78w component 3 are not the same axis.
#   Concatenating them would give a joint embedding whose structure is an
#   artefact of that mismatch.
#
#   Here the peaks are unified first, both samples are requantified against the
#   SAME features from their fragment files, and LSI is computed once on the
#   merged matrix. That is a shared space, so the joint WNN spans both ages
#   legitimately.
#
# ANNOTATION COMES FROM THE ARC REFERENCE'S OWN GTF, NOT EnsDb. Every Signac
# mouse tutorial uses EnsDb.Mmusculus.v79, which is mm10. Wiring that to GRCm39
# data gives silently wrong coordinates - the same class of failure as reusing
# -q 255 on a BWA BAM. genes/genes.gtf.gz ships with the reference and is
# GRCm39 by construction, so there is no version to match.
#
# BARCODE SPACES, THIRD INSTANCE. This project has three files keyed three
# ways, and they do not agree:
#   filtered_feature_bc_matrix / GEX BAM   -> gex_barcode (== canonical)
#   atac_fragments.tsv.gz                  -> canonical barcode  (verified:
#                                             20/20 sampled match gex_barcode,
#                                             0/20 match atac_barcode)
#   ATAC BAM (CB tag)                      -> canonical barcode, NOT
#                                             atac_barcode (verified: 25/25
#                                             distinct CB tags over three loci
#                                             match gex_barcode, 0/25 match
#                                             atac_barcode; CR holds the raw
#                                             ATAC barcode and cellranger-arc
#                                             has already translated it)
# So all three are canonical and NOTHING needs translating. This entry read
# `-> atac_barcode` until it was measured, and the sinto file exported on that
# basis selected zero reads. Getting it backwards produces an empty matrix, or
# an empty BAM, rather than an error - which is why the ATAC Allelome job
# (slurm/multiome_allelome_atac.slurm) refuses to start without a live check.
#
# Run in seurat_env, after installing Signac:
#   conda activate seurat_env
#   Rscript multiome/03_signac_joint.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat); library(Signac)
  library(GenomicRanges); library(rtracklayer)
  library(ggplot2); library(dplyr); library(readr); library(tibble); library(patchwork)
})
options(future.globals.maxSize = 32 * 1024^3)
set.seed(1)

BASE <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2"
WORK <- file.path(BASE, "adult_aged_multiome")
QC   <- file.path(WORK, "qc")
OUT  <- file.path(WORK, "figures_signac")
REF  <- file.path(BASE, "refdata-cellranger-arc-GRCm39-2024-A")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SAMPLES    <- c("9w", "78w")
KEEP_CHR   <- c(paste0("chr", 1:19), "chrX")   # female, so chrY dropped; chrM meaningless for peaks
RESOLUTION <- 0.1
N_PCS      <- 30
LSI_DIMS   <- 2:30    # LSI component 1 tracks sequencing depth; conventionally dropped
PEAK_WIDTH <- c(20, 10000)

# Core escape genes as defined in OCM_heart/core_escape_SNPs.R, plus Xist.
ESCAPE_GENES <- c("Kdm5c","Kdm6a","Ddx3x","Eif2s3x","Ftx","Jpx",
                  "Pbdc1","Utp14a","Akap17a","Sts")
PANELS <- list(
  "Ventricular Cardiomyocytes"      = c("Tnnt2","Myh6","Actc1","Ttn","Ryr2","Pln","Myl2","Tnni3","Actn2","Trdn"),
  "Cardiomyocytes (stressed)"       = c("Nppa","Nppb","Myh7","Ankrd1","Acta1","Xirp2"),
  "Fibroblasts"                     = c("Col1a1","Col3a1","Dcn","Pdgfra","Gsn","Postn","Lum"),
  "Endothelial cells"               = c("Pecam1","Cdh5","Egfl7","Kdr","Emcn","Ly6c1","Plvap","Tie1"),
  "Endocardium"                     = c("Npr3","Nfatc1","Vwf","Cytl1"),
  "Lymphatic endothelial"           = c("Lyve1","Prox1","Mmrn1","Flt4","Ccl21a"),
  "Pericytes - Smooth muscle cells" = c("Rgs5","Abcc9","Kcnj8","Myh11","Acta2","Tagln","Notch3"),
  "Macrophages"                     = c("Ptprc","Cd68","Adgre1","Csf1r","Lyz2","Mrc1","Fcgr1"),
  "T cells"                         = c("Cd3e","Cd3d","Cd8a","Themis","Skap1"),
  "B cells"                         = c("Cd79a","Cd79b","Ms4a1","Ighm"),
  "Epicardial - Mesothelial cells"  = c("Wt1","Msln","Upk3b","Krt19")
)

source(file.path(BASE, "Postdoc", "multiome", "00_helpers.R"))
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("Signac %s, Seurat %s", packageVersion("Signac"), packageVersion("Seurat"))

# ---------------------------------------------------------------------------
# annotation, GRCm39, from the reference itself
# ---------------------------------------------------------------------------
say("--- annotation from %s/genes/genes.gtf.gz ---", basename(REF))
gtf <- rtracklayer::import(file.path(REF, "genes", "genes.gtf.gz"))
gtf <- gtf[as.character(seqnames(gtf)) %in% KEEP_CHR]
# Signac wants tx_id / gene_name / gene_biotype / type on the annotation.
if (!"gene_biotype" %in% names(mcols(gtf)) && "gene_type" %in% names(mcols(gtf)))
  gtf$gene_biotype <- gtf$gene_type
if (!"tx_id" %in% names(mcols(gtf)) && "transcript_id" %in% names(mcols(gtf)))
  gtf$tx_id <- gtf$transcript_id
genome(gtf) <- "GRCm39"
say("annotation: %d ranges, %d gene names", length(gtf), length(unique(gtf$gene_name)))

# ---------------------------------------------------------------------------
# common peak set
# ---------------------------------------------------------------------------
say("--- unifying peaks ---")
read_peaks <- function(id) {
  p <- read.table(file.path(WORK, id, "outs", "atac_peaks.bed"),
                  col.names = c("chr","start","end"), comment.char = "#")
  makeGRangesFromDataFrame(p)
}
pk <- lapply(SAMPLES, read_peaks); names(pk) <- SAMPLES
for (id in SAMPLES) say("  %s: %d peaks", id, length(pk[[id]]))

combined <- reduce(do.call(c, unname(pk)))
combined <- combined[as.character(seqnames(combined)) %in% KEEP_CHR]
w <- width(combined)
combined <- combined[w > PEAK_WIDTH[1] & w < PEAK_WIDTH[2]]
seqlevels(combined) <- KEEP_CHR
say("  common peak set: %d peaks (union, reduced, %s only, width %d-%d)",
    length(combined), "chr1-19+chrX", PEAK_WIDTH[1], PEAK_WIDTH[2])
export.bed(combined, file.path(OUT, "common_peaks.bed"))

# ---------------------------------------------------------------------------
# per sample: RNA assay + ATAC requantified on the common peaks
# ---------------------------------------------------------------------------
build_one <- function(id) {
  say("--- %s ---", id)
  pass <- read_csv(file.path(QC, sprintf("qc_pass_%s.csv", id)), show_col_types = FALSE)
  h5   <- file.path(WORK, id, "outs", "filtered_feature_bc_matrix.h5")
  frag <- file.path(WORK, id, "outs", "atac_fragments.tsv.gz")

  mat <- Read10X_h5(h5); if (is.list(mat)) mat <- mat[["Gene Expression"]]
  keep <- intersect(colnames(mat), pass$barcode)
  say("  %d QC-pass nuclei", length(keep))

  # Prefix the cell names HERE rather than letting merge(add.cell.ids=) do it.
  # A Fragment object stores its own barcode->cellname map, and merge does not
  # rewrite it, so prefixing at merge time leaves the ChromatinAssay's cells
  # unprefixed while meta.data is prefixed. That mismatch is what produced
  # "No cell overlap between new meta data and Seurat object" from
  # AddModuleScore. Naming once, up front, keeps every assay consistent.
  newnames <- paste0(id, "_", keep)

  o <- CreateSeuratObject(mat[, keep], project = id, assay = "RNA")
  o <- RenameCells(o, new.names = newnames)
  o$sample     <- id
  o$percent.mt <- PercentageFeatureSet(o, pattern = "^mt-")
  meta <- pass %>% filter(barcode %in% keep) %>%
    select(barcode, atac_fragments, frip, tss_frac) %>%
    mutate(barcode = paste0(id, "_", barcode)) %>% column_to_rownames("barcode")
  o <- AddMetaData(o, meta[colnames(o), , drop = FALSE])

  # Fragments are keyed by the CANONICAL barcode (verified: 20/20 sampled match
  # gex_barcode, 0/20 atac_barcode). Signac's `cells` map takes the desired
  # object cell names as NAMES and the in-file barcodes as VALUES, so this is
  # where the prefix is declared to the fragment object.
  say("  FeatureMatrix over %d common peaks", length(combined))
  fr <- CreateFragmentObject(path = frag, cells = setNames(keep, newnames))
  cm <- FeatureMatrix(fragments = fr, features = combined, process_n = 5000)
  say("  FeatureMatrix: %d peaks x %d cells", nrow(cm), ncol(cm))

  shared <- intersect(colnames(o), colnames(cm))
  if (!length(shared)) {
    stop("no shared cells between RNA assay and FeatureMatrix for ", id,
         "\n  RNA e.g.: ", paste(head(colnames(o), 2), collapse = ", "),
         "\n  ATAC e.g.: ", paste(head(colnames(cm), 2), collapse = ", "))
  }
  say("  %d cells with both modalities", length(shared))
  o <- subset(o, cells = shared)
  o[["ATAC"]] <- CreateChromatinAssay(counts = cm[, shared, drop = FALSE],
                                      fragments = fr, annotation = gtf, genome = "GRCm39")
  o
}
objs <- lapply(SAMPLES, build_one); names(objs) <- SAMPLES

# ---------------------------------------------------------------------------
# merge, then normalise each modality ONCE over the shared features
# ---------------------------------------------------------------------------
say("--- merged: RNA PCA and ATAC LSI in a shared space ---")
# No add.cell.ids: build_one already prefixed, see the note there.
obj <- merge(objs[[1]], y = objs[[2]])

DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE)      # percent.mt NOT regressed; see 02
obj <- RunPCA(obj, npcs = N_PCS, verbose = FALSE)

DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj, verbose = FALSE)
obj <- FindTopFeatures(obj, min.cutoff = "q5")
obj <- RunSVD(obj, verbose = FALSE)
dev_open(file.path(OUT, "lsi_depth_correlation.pdf"), width = 7, height = 5)
print(DepthCor(obj) + ggtitle("LSI component vs sequencing depth"))
dev.off()
say("  DepthCor written - confirm component 1 is the depth-correlated one before")
say("  trusting LSI_DIMS = %d:%d", min(LSI_DIMS), max(LSI_DIMS))

lsi_dims <- intersect(LSI_DIMS, seq_len(ncol(Embeddings(obj, "lsi"))))

# ---------------------------------------------------------------------------
# joint WNN across BOTH samples - the thing 02 could not do
# ---------------------------------------------------------------------------
obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                               dims.list = list(1:N_PCS, lsi_dims), verbose = FALSE)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", verbose = FALSE)
obj <- FindClusters(obj, graph.name = "wsnn", resolution = RESOLUTION, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:N_PCS,
               reduction.name = "rna.umap", reduction.key = "rnaUMAP_", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "lsi", dims = lsi_dims,
               reduction.name = "atac.umap", reduction.key = "atacUMAP_", verbose = FALSE)
say("  joint clusters: %d", length(levels(obj)))

# ---- provisional cell typing on the joint clusters ----
DefaultAssay(obj) <- "SCT"
obj <- assign_celltypes(obj, PANELS, assay = "SCT", layer = "data")
obj <- set_meta(obj, "celltype_short", short_labels(obj$celltype_provisional))
CT_LEVELS <- sort(unique(obj$celltype_short))
SC_COL <- celltype_scale(CT_LEVELS, "colour")
Idents(obj) <- "celltype_provisional"
write.csv(obj@misc[["panel_cluster_means"]],
          file.path(OUT, "cluster_panel_scores.csv"))
present <- lapply(PANELS, function(g) intersect(g, rownames(obj)))
present <- present[lengths(present) > 0]

# Checkpoint here, not only at the end. FeatureMatrix plus the merged LSI is
# the multi-hour part of this script, and the last two failures both landed
# after it -- so it is written out before anything cosmetic can throw.
saveRDS(obj, file.path(WORK, "multiome_signac_joint.rds"))
say("checkpoint written: multiome_signac_joint.rds")

comp <- obj@meta.data %>% count(seurat_clusters, sample) %>%
  tidyr::pivot_wider(names_from = sample, values_from = n, values_fill = 0)
say("per-cluster sample composition - the point of a SHARED space is that these")
say("should now be mixed rather than one-sample-per-cluster:")
print(as.data.frame(comp), row.names = FALSE)
write_csv(comp, file.path(OUT, "cluster_sample_composition.csv"))

# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
# A SHARED LEGEND, not direct labels. There are 7 cell types here; direct
# labelling works up to about 4 series, and past that the labels collide with
# each other and with the data - which is exactly what happened on the first
# version of this figure, worst in the ATAC panel. One legend for three panels
# also makes the point that the three embeddings share an identity assignment.
panel <- function(red, ttl) {
  p <- DimPlot(obj, reduction = red, group.by = "celltype_short",
               pt.size = 0.3, shuffle = TRUE) +
    ggtitle(ttl) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 9))
  if (!is.null(SC_COL)) p <- p + SC_COL
  p
}
dev_open(file.path(OUT, "umap_joint.pdf"), width = 15, height = 5.5)
print(
  (panel("rna.umap", "RNA") | panel("atac.umap", "ATAC (common peaks, merged LSI)") |
   panel("wnn.umap", "WNN joint")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 9))
)
dev.off()

# The presentation panel: one large WNN plot, where 7 direct labels do fit and
# a reader does not have to move between legend and cluster.
dev_open(file.path(OUT, "umap_wnn_labelled.pdf"), width = 9, height = 7.5)
p <- DimPlot(obj, reduction = "wnn.umap", group.by = "celltype_short",
             label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.4,
             shuffle = TRUE) +
  NoLegend() + ggtitle("Joint RNA + ATAC (WNN), 9w and 78w") +
  theme(plot.title = element_text(size = 13, face = "bold"))
if (!is.null(SC_COL)) p <- p + SC_COL
print(p)
dev.off()

dev_open(file.path(OUT, "umap_joint_by_sample.pdf"), width = 11, height = 5)
print(
  (DimPlot(obj, reduction = "wnn.umap", group.by = "sample") + ggtitle("WNN, by sample")) |
  (DimPlot(obj, reduction = "wnn.umap", split.by = "sample", group.by = "celltype_short") +
     ggtitle("n=1 per age - descriptive only"))
)
dev.off()

# FindMultiModalNeighbors names its weight columns after the assays behind the
# reductions, which has varied across Seurat versions (SCT.weight/ATAC.weight
# vs pca.weight/lsi.weight). Detect them rather than hardcode: a wrong name
# here would throw at the very end of a multi-hour job.
wcols <- grep("\\.weight$", colnames(obj@meta.data), value = TRUE)
say("WNN weight columns found: %s", if (length(wcols)) paste(wcols, collapse = ", ") else "none")
if (length(wcols)) {
  dev_open(file.path(OUT, "wnn_modality_weights.pdf"), width = 9, height = 6)
  print(VlnPlot(obj, features = wcols, group.by = "celltype_short",
                pt.size = 0, ncol = 1) &
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                axis.title.x = element_blank()))
  dev.off()
}

# ---- gene activity: accessibility over gene body + promoter ----
say("--- GeneActivity ---")
DefaultAssay(obj) <- "ATAC"
ga <- try(GeneActivity(obj, features = unique(c(ESCAPE_GENES, "Xist", unlist(present)))))
if (!inherits(ga, "try-error")) {
  obj[["GeneActivity"]] <- CreateAssayObject(counts = ga)
  obj <- NormalizeData(obj, assay = "GeneActivity",
                       normalization.method = "LogNormalize", verbose = FALSE)
  say("  GeneActivity assay: %d features", nrow(ga))
  esc <- intersect(ESCAPE_GENES, rownames(obj[["GeneActivity"]]))
  if (length(esc)) {
    dev_open(file.path(OUT, "escape_gene_activity.pdf"), width = 12, height = 8)
    DefaultAssay(obj) <- "GeneActivity"
    print(FeaturePlot(obj, features = esc, reduction = "wnn.umap", ncol = 4) &
            theme(plot.title = element_text(size = 9)))
    print(DotPlot(obj, features = esc, group.by = "celltype_short") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggtitle("Core escape genes: ATAC gene activity"))
    dev.off()
  }
} else say("  GeneActivity failed, skipping - see traceback above")

# ---- coverage tracks at the escape loci ----
say("--- CoveragePlot at core escape loci ---")
DefaultAssay(obj) <- "ATAC"
avail <- intersect(c("Xist", ESCAPE_GENES), gtf$gene_name)
dev_open(file.path(OUT, "coverage_escape_loci.pdf"), width = 10, height = 7)
for (g in avail) {
  p <- try(CoveragePlot(obj, region = g, features = g, expression.assay = "SCT",
                        extend.upstream = 5000, extend.downstream = 5000), silent = TRUE)
  if (!inherits(p, "try-error")) print(p + patchwork::plot_annotation(title = g))
  else say("  CoveragePlot failed for %s", g)
}
dev.off()
say("  plotted %d loci", length(avail))

# ---------------------------------------------------------------------------
# save + export the group maps for the per-celltype allelic route
# ---------------------------------------------------------------------------
saveRDS(obj, file.path(WORK, "multiome_signac_joint.rds"))
exp <- obj@meta.data %>% rownames_to_column("cell") %>%
  mutate(barcode = sub("^(9w|78w)_", "", cell)) %>%
  select(sample, barcode, seurat_clusters, celltype_provisional)
write_csv(exp, file.path(OUT, "nucleus_celltype_joint.csv"))

for (id in SAMPLES) {
  pass <- read_csv(file.path(QC, sprintf("qc_pass_%s.csv", id)), show_col_types = FALSE)
  e <- exp %>% filter(sample == id) %>%
    inner_join(pass %>% select(barcode, gex_barcode, atac_barcode), by = "barcode") %>%
    mutate(group = gsub("[^A-Za-z0-9]+", "_", celltype_provisional))
  # BOTH BAMs are canonical-barcode-keyed - see the barcode note in the header.
  # The two files are therefore byte-identical, and that is the point: the _atac
  # one used to be written from `atac_barcode` and selected zero reads. Keeping
  # the filename so nothing that already references it breaks. Note the ATAC
  # Allelome job still reads the _rna file by name: the _atac file sitting on
  # DSS predates this fix, and will only become correct once 03 is re-run.
  write.table(e[, c("gex_barcode","group")],
              file.path(OUT, sprintf("sinto_%s_rna_bycelltype.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(e[, c("gex_barcode","group")],
              file.path(OUT, sprintf("sinto_%s_atac_bycelltype.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  say("%s: per-celltype sinto files, %d nuclei, %d groups", id, nrow(e), length(unique(e$group)))
}
say("")
say("output under %s", OUT)
say("PROVISIONAL cell types - check against the dotplot before presenting.")
