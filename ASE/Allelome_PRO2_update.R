
reads <- read.delim("read_count_per_SNP.txt")

autosomal <- subset(reads, chr != "chrX" & chr != "chrY")
chrX <- subset(reads, chr == "chrX")

gene_split <- split(autosomal, autosomal$name)

bayes_function <- function(ref, alt, alpha, beta){
    n <- ref + alt
    k <- ref
    alpha_post <- alpha + k
    beta_post  <- beta + n - k
    post_mean <- alpha_post / (alpha_post + beta_post)
    ci <- qbeta(c(0.025, 0.975), alpha_post, beta_post) 
    result_list <- list(
        posterior = paste0("Beta(", alpha_post, ", ", beta_post, ")"),
        mean = post_mean,
        ci_95 = ci)
    return(result_list)
}

gene_mean_posteriors <- list()
for(i in 1:length(gene_split)){
    tmp <- gene_split[[i]]
    mean_list <- lapply(1:nrow(tmp), function(i){
        bayes_function(tmp$A1_reads[i], tmp$A2_reads[i], 10, 10)$mean
    })
    gene_mean_posteriors[[i]] <- mean(unlist(mean_list))
}
names(gene_mean_posteriors) <- names(gene_split)

gene_mean_posteriors[unlist(lapply(gene_mean_posteriors, function(x) x < 0.3 | x > 0.7))]

gene_split_chrX <- split(chrX, chrX$name)

gene_mean_posteriors_chrX <- list()
for(i in 1:length(gene_split_chrX)){
    tmp <- gene_split_chrX[[i]]
    mean_list <- lapply(1:nrow(tmp), function(i){
        bayes_function(tmp$A1_reads[i], tmp$A2_reads[i], 0.5, 0.5)$mean
    })
    gene_mean_posteriors_chrX[[i]] <- mean(unlist(mean_list))
}
names(gene_mean_posteriors_chrX) <- names(gene_split_chrX)

gene_mean_posteriors_chrX[unlist(lapply(gene_mean_posteriors_chrX, function(x) x < 0.9 | x > 0.1))]

lapply(1:nrow(gene_split_chrX[[586]]), function(x){
    bayes_function(gene_split_chrX[[586]]$A1_reads[x], gene_split_chrX[[586]]$A2_reads[x], 0.5, 0.5)
})


# prior parameters
alpha <- 2
beta <- 2

# observed data
n <- test$A1_reads[1] + test$A2_reads[1]
k <- test$A1_reads[1]

# posterior parameters
alpha_post <- alpha + k
beta_post  <- beta + n - k

# Posterior mean
post_mean <- alpha_post / (alpha_post + beta_post)

# 95% equal-tailed credible interval
ci <- qbeta(c(0.025, 0.975), alpha_post, beta_post)

list(
  posterior = paste0("Beta(", alpha_post, ", ", beta_post, ")"),
  mean = post_mean,
  ci_95 = ci
)

pdf('../bayes.pdf')
curve(dbeta(x, alpha, beta), 0, 1, lwd=2, ylab="Density", xlab="p")
curve(dbeta(x, alpha_post, beta_post), add=TRUE, lwd=2, lty=2)
legend("topright", c("Prior", "Posterior"), lwd=2, lty=c(1,2))
dev.off()