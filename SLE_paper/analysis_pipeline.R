suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggsignif)
  library(forcats)
  library(reshape2)
  library(UpSetR)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
  library(gridExtra)
  library(speckle)
  library(beeswarm)
  library(gplots)
  library(plyr)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  n_splits = 10,
  methods = c("boruta", "enet", "intersection", "combined"),
  models = c("logit", "RF", "SVM", "GBM", "MLP"),
  root_train = "perez_update",
  root_validation = "nehar_belaid_update",
  outputs_dir = "analysis_outputs"
)

dir.create(cfg$outputs_dir, recursive = TRUE, showWarnings = FALSE)

palettes <- list(
  gene_set = c(chrX = "#8A0798", SLE = "#D90750", autosome = "#D99750", HVG.autosome = "#D99799", HVG = "#399750"),
  method = c(boruta = "#BFFB00", enet = "#B875B1", intersection = "#D2EDF6", combined = "#4DB748"),
  model = c(logit = "#F4CE03", RF = "#BCEA9D", SVM = "#99B2F5", GBM = "#F5B29E", MLP = "#26779E", ensemble = "#F5A2F5")
)

read_required_inputs <- function() {
  load("metadata_pbmc_female.control_managed.Rdata")
  load("escapees.Rdata")
  load("celltype.colours.RData")

  x_immune <- read.delim("X.immune.txt", header = FALSE)$V1
  sle_genes <- read.delim("SLE_DisGeNet.tsv")$Gene
  chrx_tbl <- read.delim("chrX_biomaRt.txt")
  chrx_genes <- chrx_tbl %>% filter(Gene.name != "") %>% pull(Gene.name)

  list(
    metadata = metadata,
    escape = escape,
    celltype_colours = colours,
    X_immune = x_immune,
    SLE = sle_genes,
    chrX = chrx_genes
  )
}

clean_celltype_label <- function(x) {
  x %>%
    gsub("^(boruta_|enet_|intersection_|combined_)?(logit_|RF_|SVM_|GBM_|MLP_)", "", .) %>%
    gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", .) %>%
    gsub("_", " ", .)
}

extract_metric_annotations <- function(celltype_col) {
  tibble(
    model = str_extract(celltype_col, "logit|RF|SVM|GBM|MLP"),
    gene.set = str_extract(celltype_col, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
    method = str_extract(celltype_col, "boruta|enet|intersection|combined")
  )
}

read_model_metrics_split <- function(split_i, cfg) {
  metric_files <- unlist(lapply(cfg$methods, function(method_i) {
    list.files(
      file.path(cfg$root_train, paste0("split_", split_i), method_i, "metrics"),
      pattern = "metrics_",
      full.names = TRUE
    )
  }))

  metric_files <- metric_files[file.exists(metric_files) & file.size(metric_files) > 0]
  if (length(metric_files) == 0) return(NULL)

  metrics <- lapply(metric_files, read.csv)
  names(metrics) <- metric_files %>%
    gsub("split_\\d+/|metrics/|.csv", "", .) %>%
    gsub("_metrics", "", .) %>%
    gsub("/", "_", .)

  bind_rows(metrics, .id = "celltype") %>%
    bind_cols(extract_metric_annotations(.$celltype)) %>%
    mutate(
      celltype = gsub("^perez_update_", "", celltype),
      celltype = clean_celltype_label(celltype)
    )
}

read_model_metrics_all <- function(cfg) {
  out <- lapply(seq_len(cfg$n_splits), function(i) read_model_metrics_split(i, cfg))
  out <- out[!sapply(out, is.null)]
  if (length(out) == 0) stop("No model metric files found.")
  bind_rows(out, .id = "split")
}

summarise_metrics <- function(df, by_model = TRUE) {
  group_vars <- c("celltype", "gene.set")
  if (by_model) group_vars <- c(group_vars, "model")

  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      Accuracy = mean(Accuracy, na.rm = TRUE),
      Precision = mean(Precision, na.rm = TRUE),
      Recall = mean(Recall, na.rm = TRUE),
      F1 = mean(F1, na.rm = TRUE),
      F1_lower = mean(F1_lower, na.rm = TRUE),
      F1_upper = mean(F1_upper, na.rm = TRUE),
      AUC = mean(AUC, na.rm = TRUE),
      AUC_lower = mean(AUC_lower, na.rm = TRUE),
      AUC_upper = mean(AUC_upper, na.rm = TRUE),
      AUPRC = mean(AUPRC, na.rm = TRUE),
      AUPRC_lower = mean(AUPRC_lower, na.rm = TRUE),
      AUPRC_upper = mean(AUPRC_upper, na.rm = TRUE),
      Kappa = mean(Kappa, na.rm = TRUE),
      MCC = mean(MCC, na.rm = TRUE),
      MCC_lower = mean(MCC_lower, na.rm = TRUE),
      MCC_upper = mean(MCC_upper, na.rm = TRUE),
      n_features = round(mean(n_features, na.rm = TRUE), 0),
      .groups = "drop"
    )
}

read_ensemble_metrics_split <- function(split_i, cfg) {
  metric_files <- unlist(lapply(cfg$methods, function(method_i) {
    list.files(
      file.path(cfg$root_train, paste0("split_", split_i), method_i, "ensemble"),
      pattern = "metrics_",
      full.names = TRUE
    )
  }))

  metric_files <- metric_files[file.exists(metric_files) & file.size(metric_files) > 0]
  if (length(metric_files) == 0) return(NULL)

  metrics <- lapply(metric_files, read.csv)
  names(metrics) <- metric_files %>%
    gsub("split_\\d+/|ensemble/metrics_|.csv", "", .) %>%
    gsub("_combined", "", .) %>%
    gsub("/", "_", .)

  bind_rows(metrics, .id = "celltype") %>%
    mutate(
      gene.set = str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
      method = str_extract(celltype, "boruta|enet|intersection|combined"),
      celltype = gsub("perez_update_", "", celltype),
      celltype = gsub("combined_|.HVG.autosome", "", celltype),
      celltype = gsub("^(boruta_|enet_|intersection_|combined_)", "", celltype),
      celltype = gsub(".SLE|.chrX|.autosome|.HVG", "", celltype),
      celltype = gsub("_", " ", celltype),
      model = "ensemble"
    )
}

read_ensemble_metrics_all <- function(cfg) {
  out <- lapply(seq_len(cfg$n_splits), function(i) read_ensemble_metrics_split(i, cfg))
  out <- out[!sapply(out, is.null)]
  if (length(out) == 0) stop("No ensemble metric files found.")
  bind_rows(out, .id = "split")
}

plot_mcc_model_comparison <- function(df, out_file, palette_model) {
  comparisons <- list(c("logit", "MLP"), c("RF", "MLP"), c("SVM", "MLP"), c("GBM", "MLP"))
  p <- ggplot(df, aes(x = model, y = MCC, colour = model)) +
    geom_jitter(width = 0.2) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
    geom_signif(comparisons = comparisons, map_signif_level = TRUE, test = "wilcox.test", color = "black", step_increase = 0.1) +
    theme_minimal() +
    labs(x = "Model", y = "MCC") +
    theme(axis.text.x = element_blank()) +
    scale_colour_manual(values = palette_model, name = "Model")

  ggsave(out_file, p, width = 8, height = 5)
}

run_chrX_sle_tests <- function(ensemble_metrics_df) {
  x <- ensemble_metrics_df %>%
    filter(gene.set %in% c("chrX", "SLE"), method == "combined")

  pvals <- lapply(split(x, x$celltype), function(d) {
    if (length(unique(d$gene.set)) < 2) return(NA_real_)
    wilcox.test(MCC ~ gene.set, data = d)$p.value
  })

  tibble(celltype = names(pvals), p_value = unlist(pvals)) %>%
    mutate(
      FDR = p.adjust(p_value, method = "fdr"),
      p_signif = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

plot_chrX_sle_panel <- function(ensemble_metrics_df, pvals_df, out_file, palette_gene_set) {
  x <- ensemble_metrics_df %>% filter(gene.set %in% c("chrX", "SLE"), method == "combined")
  ord <- pvals_df %>% arrange(p_value) %>% pull(celltype)
  x$celltype <- fct_relevel(x$celltype, rev(ord))

  p <- ggplot(x, aes(y = MCC, x = celltype, fill = gene.set)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 1) +
    geom_text(data = pvals_df, aes(y = Inf, x = celltype, label = p_signif), hjust = 1.2, vjust = 0.5, size = 3, inherit.aes = FALSE) +
    geom_hline(yintercept = 0.7, linetype = "dotted", color = "red") +
    theme_minimal() +
    labs(y = "MCC", x = "") +
    theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12, angle = 90), legend.position = "top") +
    scale_fill_manual(values = palette_gene_set, name = "Gene Set")

  ggsave(out_file, p, width = 13, height = 6)
}

read_feature_files_split <- function(split_i, cfg) {
  files <- unlist(lapply(cfg$methods, function(method_i) {
    list.files(
      file.path(cfg$root_train, paste0("split_", split_i), "features"),
      pattern = method_i,
      full.names = TRUE
    )
  }))

  files <- files[file.exists(files) & file.size(files) > 0]
  if (length(files) == 0) return(NULL)

  features <- lapply(files, read.csv)
  names(features) <- files %>%
    gsub("split_\\d+/|metrics/|.csv", "", .) %>%
    gsub("_metrics", "", .) %>%
    gsub("/", "_", .)

  bind_rows(features, .id = "celltype") %>%
    mutate(
      gene.set = str_extract(celltype, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
      method = str_extract(celltype, "boruta|enet|intersection|combined"),
      celltype = gsub("^perez_update_features_", "", celltype),
      celltype = gsub("features.", "", celltype),
      celltype = gsub("(boruta_|enet_|intersection_|combined_)", "", celltype),
      celltype = gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", celltype),
      celltype = gsub("_", " ", celltype)
    )
}

build_feature_frequency_matrix <- function(features_df, celltypes) {
  features_temp <- lapply(celltypes, function(ct) {
    features_df %>% filter(method == "combined", celltype == ct)
  })
  names(features_temp) <- celltypes

  features_by_gene_set <- lapply(celltypes, function(ct) {
    tapply(features_temp[[ct]]$Feature, features_temp[[ct]]$gene.set, plyr::count)
  })
  names(features_by_gene_set) <- celltypes

  all_features <- unique(unlist(lapply(celltypes, function(ct) features_by_gene_set[[ct]]$chrX[, 1])))
  temp <- do.call(rbind, lapply(celltypes, function(ct) cbind(ct, features_by_gene_set[[ct]]$chrX)))

  m <- spread(temp, key = 1, value = 3, fill = 0)
  rownames(m) <- m[, 1]
  mx <- apply(m[, -1], 2, as.numeric)
  rownames(mx) <- rownames(m)

  list(all_features = all_features, feature_freq_mtx = mx, features_by_gene_set = features_by_gene_set)
}

plot_feature_heatmaps <- function(freq_mtx, escape, celltype_colours, out_binary, out_counts) {
  x_escape <- rownames(escape)[escape$Xie_Tukiainen2017 == 1]
  x_escape2 <- rownames(escape)[escape$Xive_Tukiainen2017 == 1]

  xi_cols <- rep("grey", nrow(freq_mtx))
  xi_cols[!is.na(match(rownames(freq_mtx), x_escape))] <- "purple"
  xi_cols[!is.na(match(rownames(freq_mtx), x_escape2))] <- "darkorchid4"

  top_binary <- 1 * (freq_mtx == 10)
  keep_rows <- rowSums(top_binary) > 0

  pdf(out_binary, width = 12, height = 8)
  heatmap.2(top_binary[keep_rows, ], density.info = "none", trace = "none", col = c(0, 1), RowSideColors = xi_cols[keep_rows], ColSideColors = celltype_colours[colnames(top_binary)])
  dev.off()

  cols_grey <- colorRampPalette(c("white", "black"))(10)
  pdf(out_counts, width = 12, height = 8)
  heatmap.2(freq_mtx[keep_rows, ], density.info = "none", trace = "none", col = cols_grey, RowSideColors = xi_cols[keep_rows], ColSideColors = celltype_colours[colnames(freq_mtx)])
  dev.off()
}

read_validation_metrics <- function(cfg) {
  metric_files <- list.files(cfg$root_validation, pattern = "metrics_", full.names = TRUE)
  metric_files <- metric_files[file.exists(metric_files) & file.size(metric_files) > 0]
  if (length(metric_files) == 0) stop("No validation metric files found.")

  metrics <- lapply(metric_files, read.csv)
  names(metrics) <- metric_files %>%
    gsub("split_\\d+/|metrics/|.csv", "", .) %>%
    gsub("_metrics", "", .) %>%
    gsub("/", "_", .)

  bind_rows(metrics, .id = "celltype") %>%
    mutate(
      gene.set = str_extract(celltype, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
      split = str_extract(celltype, "split_[0-9]*"),
      age_group = str_extract(celltype, "adult|child"),
      celltype = gsub("^nehar_belaid_update_metrics_", "", celltype),
      celltype = gsub("_child_split_[0-9]*", "", celltype),
      celltype = gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", celltype),
      celltype = gsub("adult|child|split_[0-9]*", "", celltype),
      celltype = gsub("_", " ", celltype),
      celltype = gsub("  ", " ", celltype),
      celltype = gsub(" $", "", celltype)
    )
}

plot_validation_mcc <- function(metrics_df, out_file, palette_gene_set) {
  p <- ggplot(metrics_df, aes(y = MCC, x = celltype, fill = gene.set)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 1) +
    geom_hline(yintercept = 0.7, linetype = "dotted", color = "red") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey") +
    facet_wrap(~age_group, ncol = 1) +
    theme_minimal() +
    labs(y = "MCC", x = "") +
    theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12, angle = 90), legend.position = "top") +
    scale_fill_manual(values = palette_gene_set, name = "Gene Set")

  ggsave(out_file, p, width = 13, height = 10)
}

main <- function(cfg, palettes) {
  inp <- read_required_inputs()

  model_metrics <- read_model_metrics_all(cfg)
  model_metrics$model <- factor(model_metrics$model, levels = cfg$models)
  model_metrics$method <- factor(model_metrics$method, levels = cfg$methods)

  celltypes <- sort(unique(model_metrics$celltype))
  celltype_colours <- inp$celltype_colours
  if (!all(celltypes %in% names(celltype_colours))) {
    missing_ct <- setdiff(celltypes, names(celltype_colours))
    add_cols <- viridisLite::turbo(length(missing_ct))
    names(add_cols) <- missing_ct
    celltype_colours <- c(celltype_colours, add_cols)
  }

  model_metrics$celltype <- factor(model_metrics$celltype, levels = celltypes)
  model_metrics$split <- factor(model_metrics$split, levels = as.character(seq_len(cfg$n_splits)))

  save(model_metrics, file = file.path(cfg$outputs_dir, "models_metrics_df.RData"))

  avg_by_method <- lapply(cfg$methods, function(m) {
    model_metrics %>% filter(method == m) %>% summarise_metrics(by_model = TRUE)
  })
  names(avg_by_method) <- cfg$methods

  plot_mcc_model_comparison(
    avg_by_method$combined,
    out_file = file.path(cfg$outputs_dir, "compare_models_mcc_combined.pdf"),
    palette_model = palettes$model
  )

  ensemble_metrics <- read_ensemble_metrics_all(cfg)
  ensemble_metrics$celltype <- factor(ensemble_metrics$celltype, levels = celltypes)

  avg_ensemble_by_method <- lapply(cfg$methods, function(m) {
    ensemble_metrics %>% filter(method == m) %>% summarise_metrics(by_model = TRUE)
  })
  names(avg_ensemble_by_method) <- cfg$methods

  pvals_chrX_sle <- run_chrX_sle_tests(ensemble_metrics)
  write.csv(pvals_chrX_sle, file.path(cfg$outputs_dir, "chrX_vs_SLE_wilcox_results.csv"), row.names = FALSE)

  plot_chrX_sle_panel(
    ensemble_metrics_df = ensemble_metrics,
    pvals_df = pvals_chrX_sle,
    out_file = file.path(cfg$outputs_dir, "compare_avg_mcc_combined_X_SLE.pdf"),
    palette_gene_set = palettes$gene_set
  )

  feature_splits <- lapply(seq_len(cfg$n_splits), function(i) read_feature_files_split(i, cfg))
  feature_splits <- feature_splits[!sapply(feature_splits, is.null)]
  if (length(feature_splits) > 0) {
    features_df <- bind_rows(feature_splits, .id = "split")
    features_df$celltype <- factor(features_df$celltype, levels = celltypes)

    freq_obj <- build_feature_frequency_matrix(features_df, celltypes)
    plot_feature_heatmaps(
      freq_mtx = freq_obj$feature_freq_mtx,
      escape = inp$escape,
      celltype_colours = celltype_colours,
      out_binary = file.path(cfg$outputs_dir, "all_celltypes_top_features_heatmap.pdf"),
      out_counts = file.path(cfg$outputs_dir, "all_celltypes_all_features_heatmap.pdf")
    )
  }

  validation_metrics <- read_validation_metrics(cfg)
  plot_validation_mcc(
    validation_metrics,
    out_file = file.path(cfg$outputs_dir, "validation_mcc.pdf"),
    palette_gene_set = palettes$gene_set
  )

  message("Analysis completed. Outputs in: ", cfg$outputs_dir)
}

main(cfg, palettes)
