### Function to perform a bayesian update of allele specifc expression data ###
### Here we make some assumptions about the data. 
### For autosomal genes, we assume that the allelic ratio is centered around 0.5 with some flexibility
### For X chromosome genes in females, we assume that the allelic ratio is biased towards 1 or 0
### Therefore, we can use a beta distribution to model the allelic ratio. For autosomal genes a prior of Beta(10, 10) fits our expectation.
### Whereas, for X chromosome genes in females, we can use a prior of Beta(0.5, 0.5) to reflect the expected bias.

### NOTE: I restricted the use for genes with > 10 SNPs but this can be relaxed.

# The bayes function requires the following parameters:
# - ref: number of reference reads
# - alt: number of alternative reads
# - alpha: prior alpha parameter for the beta distribution
# - beta: prior beta parameter for the beta distribution
# - null: null hypothesis value for the allelic ratio
# - tail: direction of the test (two.sided, less, greater)

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
    p <- 1 - pbeta(null, a, b)
  } else if (tail == "less") {
    p <- pbeta(null, a, b)
  } else { # two.sided
    q_gt <- 1 - pbeta(null, a, b)
    p <- 2 * pmin(q_gt, 1 - q_gt)
  }
  p <- pmax(p, .Machine$double.xmin)

  mag <- -log10(p)
  score <- sign(post_mean - null) * mag

  list(
    posterior = paste0("Beta(", a, ", ", b, ")"),
    mean = post_mean,
    ci_95 = ci,
    p_post = p,
    score_postP = score
  )
}

# The weighted score aggregates the posterior scores using the total read counts as weights
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

# The weighted mean aggregates the posterior means using the total read counts as weights
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


# # Example use
# library(ggplot2)
# library(tidyr)
# library(dplyr)

# ref_reads <- c(18, 25, 30, 22, 27, 19, 35, 40, 28, 32)
# alt_reads <- c(12, 20, 15, 18, 23, 21, 30, 35, 25, 28)

# results <- bayes_function(ref_reads, alt_reads, alpha = 10, beta = 10, null = 0.5, tail = "two.sided")

# weighted_score(results, ref_reads, alt_reads)
# weighted_mean(results, ref_reads, alt_reads)


# k <- sum(ref_reads)
# n <- sum(ref_reads + alt_reads)

# # Prior parameters
# alpha <- 10
# beta <- 10

# # Posterior parameters
# post_alpha <- alpha + k
# post_beta <- beta + n - k

# # Likelihood parameters (for visualization as Beta)
# like_alpha <- k + 1
# like_beta <- n - k + 1

# # Create theta grid
# theta <- seq(0, 1, length.out = 1000)

# # Calculate densities
# prior <- dbeta(theta, alpha, beta)
# likelihood <- dbeta(theta, like_alpha, like_beta)
# posterior <- dbeta(theta, post_alpha, post_beta)

# # Combine into a data frame for ggplot
# df <- tibble(
#   theta = theta,
#   Prior = prior,
#   Likelihood = likelihood,
#   Posterior = posterior
# ) %>%
#   pivot_longer(-theta, names_to = "Distribution", values_to = "Density")

# # Plot
# pdf('Bayes_example.pdf')
# ggplot(df, aes(x = theta, y = Density, color = Distribution)) +
#   geom_line(size = 1.2) +
#   labs(
#     title = "Prior, Likelihood, and Posterior Distributions",
#     x = "Allelic Ratio",
#     y = "Density"
#   ) +
#   theme_minimal() +
#   scale_color_manual(values = c("Prior" = "blue", "Likelihood" = "orange", "Posterior" = "red"))
# dev.off()