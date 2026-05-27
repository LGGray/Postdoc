
library(DESeq2)
library(ggplot2)


# Compare gene-gene correlation rewiring between aging (aged-adult)
# and TAC injury (TAC-sham) within each cell type.

output_root <- "LINKS_study/gene_correlation"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

celltype_adult_map <- c(
  MP = "Macrophages",
  CF = "Cardiac fibroblasts",
  EC = "Endothelial cells",
  CM = "Cardiomyocytes"
)

read_counts <- function(path) {
  read.delim(path, header = TRUE, row.names = 1, check.names = FALSE)
}

load_celltype_data <- function(cell_type, adult_label) {
  adult <- read_counts(paste0("cardiac_RNAseq/DEG/", adult_label, "_GEM.txt"))
  aged <- read_counts(paste0("aged_cardiac_RNAseq/DEG/", cell_type, "_GEM.txt"))
  tac <- read_counts(paste0("TAC_cardiac_RNAseq/DEG/TAC_", cell_type, "_GEM.txt"))
  sham <- read_counts(paste0("TAC_cardiac_RNAseq/DEG/Sham_", cell_type, "_GEM.txt"))

  stopifnot(identical(rownames(adult), rownames(aged)))
  stopifnot(identical(rownames(tac), rownames(sham)))
  stopifnot(identical(rownames(adult), rownames(tac)))

  counts <- data.matrix(cbind(adult, aged, sham, tac))
  storage.mode(counts) <- "integer"

  metadata <- data.frame(
    row.names = colnames(counts),
    condition = factor(
      c(
        rep("adult", ncol(adult)),
        rep("aged", ncol(aged)),
        rep("sham", ncol(sham)),
        rep("TAC", ncol(tac))
      ),
      levels = c("adult", "aged", "sham", "TAC")
    )
  )

  list(counts = counts, metadata = metadata)
}

compute_vst <- function(counts, metadata, min_count = 10, min_samples = 3) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ 1
  )

  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds <- dds[keep, ]

  assay(vst(dds, blind = TRUE))
}

select_variable_genes <- function(vst_mat, top_n = 2000) {
  if (nrow(vst_mat) <= top_n) {
    return(rownames(vst_mat))
  }

  vars <- apply(vst_mat, 1, var)
  names(sort(vars, decreasing = TRUE))[seq_len(top_n)]
}

cor_matrix <- function(vst_mat, metadata, group_label, method = "spearman") {
  idx <- metadata$condition == group_label
  group_mat <- vst_mat[, idx, drop = FALSE]
  cor(t(group_mat), method = method, use = "pairwise.complete.obs")
}

matrix_to_edges <- function(delta_age, delta_tac) {
  idx <- which(upper.tri(delta_age), arr.ind = TRUE)

  edge_df <- data.frame(
    gene1 = rownames(delta_age)[idx[, 1]],
    gene2 = colnames(delta_age)[idx[, 2]],
    delta_age = delta_age[idx],
    delta_TAC = delta_tac[idx],
    stringsAsFactors = FALSE
  )

  edge_df$abs_delta_age <- abs(edge_df$delta_age)
  edge_df$abs_delta_TAC <- abs(edge_df$delta_TAC)
  edge_df$direction <- ifelse(sign(edge_df$delta_age) == sign(edge_df$delta_TAC), "same", "opposite")
  edge_df$joint_score <- edge_df$abs_delta_age * edge_df$abs_delta_TAC

  edge_df[order(edge_df$joint_score, decreasing = TRUE), ]
}

plot_delta_scatter <- function(edge_df, out_file, max_points = 150000) {
  plot_df <- edge_df
  if (nrow(plot_df) > max_points) {
    set.seed(42)
    plot_df <- plot_df[sample(seq_len(nrow(plot_df)), max_points), ]
  }

  p <- ggplot(plot_df, aes(x = delta_age, y = delta_TAC, color = direction)) +
    geom_point(alpha = 0.25, size = 0.7) +
    scale_color_manual(values = c(same = "#1f78b4", opposite = "#e31a1c")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(
      x = "Delta correlation (aged - adult)",
      y = "Delta correlation (TAC - sham)",
      color = "Direction"
    )

  ggsave(out_file, p, width = 6, height = 5, dpi = 300)
}

run_one_celltype <- function(cell_type, adult_label, top_n_var = 2000) {
  out_dir <- file.path(output_root, cell_type)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- load_celltype_data(cell_type, adult_label)
  vst_mat <- compute_vst(dat$counts, dat$metadata)

  selected <- select_variable_genes(vst_mat, top_n = top_n_var)
  vst_sel <- vst_mat[selected, , drop = FALSE]

  c_adult <- cor_matrix(vst_sel, dat$metadata, "adult")
  c_aged <- cor_matrix(vst_sel, dat$metadata, "aged")
  c_sham <- cor_matrix(vst_sel, dat$metadata, "sham")
  c_tac <- cor_matrix(vst_sel, dat$metadata, "TAC")

  delta_age <- c_aged - c_adult
  delta_tac <- c_tac - c_sham

  write.table(c_adult, file = file.path(out_dir, "cor_adult.tsv"), sep = "\t", quote = FALSE)
  write.table(c_aged, file = file.path(out_dir, "cor_aged.tsv"), sep = "\t", quote = FALSE)
  write.table(c_sham, file = file.path(out_dir, "cor_sham.tsv"), sep = "\t", quote = FALSE)
  write.table(c_tac, file = file.path(out_dir, "cor_TAC.tsv"), sep = "\t", quote = FALSE)
  write.table(delta_age, file = file.path(out_dir, "delta_age.tsv"), sep = "\t", quote = FALSE)
  write.table(delta_tac, file = file.path(out_dir, "delta_TAC.tsv"), sep = "\t", quote = FALSE)

  edges <- matrix_to_edges(delta_age, delta_tac)
  write.table(edges, file = file.path(out_dir, "edge_rewiring.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  top_shared <- head(edges, 1000)
  write.table(top_shared, file = file.path(out_dir, "edge_rewiring_top1000.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  plot_delta_scatter(edges, file.path(out_dir, "delta_age_vs_delta_TAC.png"))

  concordance <- suppressWarnings(cor.test(edges$delta_age, edges$delta_TAC, method = "spearman", exact = FALSE))

  data.frame(
    cell_type = cell_type,
    n_genes = nrow(vst_sel),
    n_edges = nrow(edges),
    rho_delta_age_vs_TAC = as.numeric(concordance$estimate),
    p_value = concordance$p.value,
    same_direction_frac = mean(edges$direction == "same"),
    stringsAsFactors = FALSE
  )
}

summary_list <- lapply(names(celltype_adult_map), function(ct) {
  message("Running rewiring analysis for ", ct)
  run_one_celltype(ct, celltype_adult_map[[ct]])
})

summary_df <- do.call(rbind, summary_list)
write.table(summary_df, file = file.path(output_root, "rewiring_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Done. Results written to: ", output_root)
