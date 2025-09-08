# # Reading in all SNP files 

# SNP_files_path <- list.files('.', pattern='read_count_per_SNP.txt', recursive=TRUE)
# SNP_files <- lapply(SNP_files_path, function(x){
#     read.delim(x)
# })
# names(SNP_files) <- lapply(SNP_files_path, function(x){
#     strsplit(dirname(x), '_')[[1]][4]
# })


# metadata <- read.csv('../SraRunTable.csv')
# metadata <- subset(metadata, Run %in% names(SNP_files))
# metadata <- metadata[match(metadata$Run, names(SNP_files)),]
# metadata <- metadata[, c('Run', 'AGE', 'genotype', 'sex', 'tissue')]
# colnames(metadata)[1] <- 'sample'



bayes_function <- function(ref, alt, alpha, beta, null = 0.5, tail = c("two.sided","less","greater")){
  tail <- match.arg(tail)
  n <- ref + alt
  k <- ref
  a <- alpha + k
  b <- beta  + (n - k)

  post_mean <- a / (a + b)
  ci <- qbeta(c(0.025, 0.975), a, b)

  # posterior tail(s) at 'null'
  if (tail == "greater") {
    p <- 1 - pbeta(null, a, b)              # P(p > null)
  } else if (tail == "less") {
    p <- pbeta(null, a, b)                   # P(p < null)
  } else { # two.sided
    q_gt <- 1 - pbeta(null, a, b)
    p <- 2 * pmin(q_gt, 1 - q_gt)
  }
  p <- pmax(p, .Machine$double.xmin)  # avoid 0

  mag <- -log10(p)                         # magnitude
  score <- sign(post_mean - null) * mag         # signed score

  list(
    posterior = paste0("Beta(", a, ", ", b, ")"),
    mean = post_mean,
    ci_95 = ci,
    p_post = p,
    score_postP = score
  )
}

weighted_score <- function(results_list, ref_reads, alt_reads) {
  # Extract the posterior scores from the results list
  post_scores <- results_list$score_postP
  
  # Calculate the total read counts for each SNP to use as weights
  total_reads <- ref_reads + alt_reads
  
  # Calculate the weighted average of the scores
  weighted_score <- sum(post_scores * total_reads) / sum(total_reads)
  
  # Return the final weighted score
  return(weighted_score)
}

weighted_mean <- function(results_list, ref_reads, alt_reads) {
  # Extract the posterior means from the results list
  post_means <- results_list$mean
  
  # Calculate the total read counts for each SNP to use as weights
  total_reads <- ref_reads + alt_reads
  
  # Calculate the weighted mean
  weighted_mean <- sum(post_means * total_reads) / sum(total_reads)
  
  # Return the final weighted mean
  return(weighted_mean)
}

# final_list <- list()
# for(file in names(SNP_files)){
#     input <- SNP_files[[file]]

#     input_split_by_gene <- split(input, input$name)

#     result_list <- lapply(input_split_by_gene, function(x){
#     if(nrow(x) < 10){
#         return(data.frame())
#     }
#     if (x$chr[1] != 'chrX') {
#     result <- bayes_function(x$A1_reads, x$A2_reads, 10, 10, null = 0.5, tail = "two.sided")
#     data.frame(weighted_score = weighted_score(result, x$A1_reads, x$A2_reads), 
#               weighted_mean = weighted_mean(result, x$A1_reads, x$A2_reads), 
#               chr = x$chr[1])
#     } else {
#     result <- bayes_function(x$A1_reads, x$A2_reads, 0.5, 0.5, null = 0.5, tail = "two.sided")
#     data.frame(weighted_score = weighted_score(result, x$A1_reads, x$A2_reads), 
#               weighted_mean = weighted_mean(result, x$A1_reads, x$A2_reads),
#               chr = x$chr[1])
#         }
#     })
#     names(result_list) <- names(input_split_by_gene)
#     result_list <- result_list[sapply(result_list, function(x) nrow(x) > 0)]
#     result_df <- dplyr::bind_rows(result_list, .id='name')

#     final_list[[file]] <- result_df
#     print(paste('done', file))
# }

# save(final_list, file='analysis/bayesian_update.RData')


# ## Read in locus files
# locus_files_path <- list.files('.', pattern='locus_table.txt', recursive=TRUE)
# locus_files <- lapply(locus_files_path, function(x){
#     read.delim(x)
# })
# names(locus_files) <- lapply(locus_files_path, function(x){
#     strsplit(dirname(x), '_')[[1]][4]
# })

# subset(metadata, sex == 'female')

# merged <- merge(final_list[['SRR30223464']], locus_files[['SRR30223464']], by='name')

# library(ggplot2)
# library(dplyr)
# library(tidyr)

# # Put into long format for easy plotting
# df_long <- merged %>%
#   select(name, weighted_score, allelic_score, total_reads) %>%
#   pivot_longer(cols = c(weighted_score, allelic_score),
#                names_to = "method", values_to = "score")

# library(ggpubr)

# pdf('analysis/allelic_score_vs_weighted_score.pdf')
# # add spearman correlation
# ggplot(merged, aes(x = allelic_score, y = weighted_score)) +
#   geom_point(alpha = 0.5) +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   stat_cor(method = "spearman")
# dev.off()

# pdf('analysis/allelic_ratio_vs_weighted_mean.pdf')
# ggplot(merged, aes(x = allelic_ratio, y = weighted_mean)) +
#   geom_point(alpha = 0.5) +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   stat_cor(method = "spearman")
# dev.off()


# metadata

# all_combined <- bind_rows(final_list, .id='sample')

# ## Match all_combined$sample to metadata$organ
# all_combined <- left_join(all_combined, metadata, by = 'sample')

# all_combined_female_chrX <- subset(all_combined, sex == 'female' & chr == 'chrX')

# # Plot violin plot of weighted mean for each organ
# library(ggplot2)

# pdf('analysis/weighted_mean_by_tissue_female_chrX.pdf')
# ggplot(all_combined_female_chrX, aes(x = tissue, y = weighted_mean, fill = tissue)) +
#   geom_violin() +
#   theme_minimal() +
#   labs(title = "",
#        x = "",
#        y = "Weighted Mean") +
#   theme(legend.position = "none")
# dev.off()



# input <- read.delim('read_count_per_SNP.txt')