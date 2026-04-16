library(ggplot2)

source('../../Postdoc/Allelome_PRO2_Bayes.R')

metadata <- read.csv('../SraRunTable.csv')
metadata <- subset(metadata, sex == 'female')
mean_posterior_AR <- lapply(split(metadata, metadata$source_name), function(ct) {
  per_sample <- lapply(ct$Run, function(sample) {
    locus_table <- read.delim(
      paste0("../", sample, "/2025_11_20_", sample, "_Aligned.sortedByCoord.out_annotation_us.bed_1/locus_table.txt")
    )
    locus_table <- subset(locus_table, total_reads >= 20 & chr == "chrX")

    results <- bayes_function(
      locus_table$A1_reads, locus_table$A2_reads,
      alpha = 19, beta = 1, null = 0.9, tail = "less"
    )

    data.frame(
      name = locus_table$name,
      sample = sample,
      prior_AR = locus_table$allelic_ratio,
      posterior_AR = results$mean,
      p_post = results$p_post,
      depth = locus_table$total_reads
    )
  })

  per_sample <- do.call(rbind, per_sample)

  # Unweighted mean across 3 samples
  out <- aggregate(posterior_AR ~ name, data = per_sample, FUN = function(x) mean(x, na.rm = TRUE))
  names(out)[2] <- "mean_posterior_AR"

  # Optional: weighted mean by read depth (often better)
  w <- aggregate(cbind(num = posterior_AR * depth, den = depth) ~ name, data = per_sample, FUN = sum, na.rm = TRUE)
  out$weighted_mean_posterior_AR <- w$num / w$den

  # Mean prior AR across samples (for reference)
  prior_agg <- aggregate(prior_AR ~ name, data = per_sample, FUN = function(x) mean(x, na.rm = TRUE))
  names(prior_agg)[2] <- "mean_prior_AR"
  out <- merge(out, prior_agg, by = "name")

  # Fisher's method to combine p_post across samples
  fisher_agg <- aggregate(p_post ~ name, data = per_sample, FUN = function(p) {
    k <- sum(!is.na(p))
    stat <- -2 * sum(log(p[!is.na(p)]))
    pchisq(stat, df = 2 * k, lower.tail = FALSE)
  })
  names(fisher_agg)[2] <- "fisher_p"
  out <- merge(out, fisher_agg, by = "name")

  out
})
names(mean_posterior_AR) <- names(split(metadata, metadata$source_name))

# Male samples analysis (same as female)
metadata_male <- read.csv('../SraRunTable.csv')
metadata_male <- subset(metadata_male, sex == 'male')

male_raw_AR <- lapply(split(metadata_male, metadata_male$source_name), function(ct) {
  per_sample <- lapply(ct$Run, function(sample) {
    locus_table <- read.delim(
      paste0("../", sample, "/2025_11_20_", sample, "_Aligned.sortedByCoord.out_annotation_us.bed_1/locus_table.txt")
    )
    locus_table <- subset(locus_table, total_reads >= 20 & chr == "chrX")
    data.frame(name = locus_table$name, raw_AR = locus_table$allelic_ratio)
  })
  per_sample <- do.call(rbind, per_sample)
  out <- aggregate(raw_AR ~ name, data = per_sample, FUN = function(x) mean(x, na.rm = TRUE))
  names(out)[2] <- "mean_raw_AR"
  out
})

# Get male false positives using raw AR
male_ar_named <- mapply(function(df, nm) {
  names(df)[names(df) == "mean_raw_AR"] <- nm
  df
}, male_raw_AR, names(male_raw_AR), SIMPLIFY = FALSE)

male_genes_combined <- Reduce(function(x, y) merge(x, y, by = "name", all = TRUE), male_ar_named)
male_genes_combined$mean <- rowMeans(male_genes_combined[, -1, drop = FALSE], na.rm = TRUE)
false_positives <- male_genes_combined$name[male_genes_combined$mean <= 0.9]

# Filter female escape genes
female_escape_genes <- lapply(mean_posterior_AR, function(df) {
  subset(df, weighted_mean_posterior_AR < 0.9 & fisher_p > 0.95 & !(name %in% false_positives))
})

lapply(female_escape_genes, function(x) cor.test(x$mean_prior_AR, x$weighted_mean_posterior_AR))

plot_df <- do.call(rbind, Map(function(df, ct) {
  if (nrow(df) == 0) return(NULL)
  df$cell_type <- ct
  df
}, female_escape_genes, names(female_escape_genes)))

pdf('prior_vs_posterior_AR_escape_genes_all_celltypes_beta_19_1.pdf', width = 10, height = 8)
ggplot(plot_df, aes(x = mean_prior_AR, y = weighted_mean_posterior_AR)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.9, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dashed") +
  facet_wrap(~ cell_type, ncol = 2) +
  labs(
    x = "Mean raw AR",
    y = "Weighted posterior AR",
    title = "Prior vs Posterior AR (escape genes) across cell types beta(19,1)"
  ) +
  theme_bw()
dev.off()

# Visualize prior, likelihood, posterior for one example
plot_prior_likelihood_posterior <- function(ref, alt, alpha = 19, beta = 1,
                                            gene_label = "example_gene",
                                            sample_label = "example_sample",
                                            outfile = "prior_likelihood_posterior_example.pdf") {
  k <- ref
  n <- ref + alt

  post_alpha <- alpha + k
  post_beta <- beta + (n - k)

  theta <- seq(0, 1, length.out = 1000)
  prior <- dbeta(theta, alpha, beta)
  likelihood <- dbeta(theta, k + 1, n - k + 1)
  posterior <- dbeta(theta, post_alpha, post_beta)

  d <- rbind(
    data.frame(theta = theta, density = prior, distribution = "Prior"),
    data.frame(theta = theta, density = likelihood, distribution = "Likelihood"),
    data.frame(theta = theta, density = posterior, distribution = "Posterior")
  )

  pdf(outfile, width = 8, height = 5)
  print(
    ggplot(d, aes(x = theta, y = density, color = distribution)) +
      geom_line(linewidth = 1.1) +
      geom_vline(xintercept = 0.9, linetype = "dashed", color = "black") +
      scale_color_manual(values = c(Prior = "#1f77b4", Likelihood = "#ff7f0e", Posterior = "#d62728")) +
      labs(
        title = paste0("Prior, Likelihood, Posterior: ", gene_label),
        subtitle = paste0(sample_label, " (A1=", k, ", A2=", n - k, ")"),
        x = "Allelic ratio (theta)",
        y = "Density"
      ) +
      theme_bw()
  )
  dev.off()
}

# Real-data example: choose one chrX gene from first female sample
example_sample <- metadata$Run[1]
example_locus <- read.delim(
  paste0("../", example_sample, "/2025_11_20_", example_sample, "_Aligned.sortedByCoord.out_annotation_us.bed_1/locus_table.txt")
)
example_locus <- subset(example_locus, total_reads >= 20 & chr == "chrX")

if (nrow(example_locus) > 0) {
  example_results <- bayes_function(
    example_locus$A1_reads, example_locus$A2_reads,
    alpha = 19, beta = 1, null = 0.9, tail = "less"
  )
  raw_ar <- example_locus$allelic_ratio
  post_ar <- example_results$mean
  shift <- abs(raw_ar - post_ar)

  # Prefer threshold-crossing examples for presentation: raw > 0.9 and posterior < 0.9
  idx_cross <- which(raw_ar > 0.9 & post_ar < 0.9)
  if (length(idx_cross) > 0) {
    i <- idx_cross[which.max(shift[idx_cross])]
  } else if (any(shift >= 0.10)) {
    idx_shift <- which(shift >= 0.10)
    i <- idx_shift[which.max(shift[idx_shift])]
  } else {
    i <- which.max(shift)
  }

  plot_prior_likelihood_posterior(
    ref = example_locus$A1_reads[i],
    alt = example_locus$A2_reads[i],
    alpha = 19,
    beta = 1,
    gene_label = paste0(
      as.character(example_locus$name[i]),
      " (raw=", round(raw_ar[i], 3),
      ", post=", round(post_ar[i], 3),
      ", |delta|=", round(shift[i], 3), ")"
    ),
    sample_label = as.character(example_sample),
    outfile = "prior_likelihood_posterior_real_example.pdf"
  )
}

# Synthetic example for presentation slides
plot_prior_likelihood_posterior(
  ref = 27,
  alt = 9,
  alpha = 19,
  beta = 1,
  gene_label = "Synthetic example",
  sample_label = '',
  outfile = "prior_likelihood_posterior_synthetic_example.pdf"
)

