library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(ggalluvial)
library(forcats)
library(broom)
library(lme4)
library(emmeans)
library(ggpubr)

setwd("/Users/graylachlan/LRZ Sync+Share/LGray/GTEx")

gtex_link <- data.frame(read_excel('07_SUPPLEMENTARY_TABLE.xlsx', sheet = 17))
gtex_link$Link <- paste(gtex_link$ncRNA, gtex_link$Name_pcGene, sep = "_")

metadata <- read.delim('GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt')

# Extract the relevant part of the Sample
gtex_link$SUBJID <- unlist(lapply(gtex_link$Sample, function(x) {
    paste(strsplit(x, '-')[[1]][1], strsplit(x, '-')[[1]][2], sep = '-')
}))
# Merge gtex_link with metadata to add AGE
gtex_link <- merge(gtex_link, metadata[, c("SUBJID", "AGE", "SEX")], by = "SUBJID", all.x = TRUE)

# Set factors for age
gtex_link$AGE <- factor(gtex_link$AGE, levels = c('20-29', '30-39', '40-49', '50-59', '60-69', '70-79'))
gtex_link$SEX <- ifelse(gtex_link$SEX == 1, 'Male', 'Female')

# Collapse ages to increase sample size
gtex_link$AGE_collapsed <- ifelse(gtex_link$AGE %in% c('20-29', '30-39', '40-49'), '20-49',
                              ifelse(gtex_link$AGE %in% c('60-69', '70-79'), '60-79', '50-59'))
gtex_link$AGE_collapsed <- factor(gtex_link$AGE_collapsed, levels = c('20-49', '50-59', '60-79'))


# Count the number of enhancing and repressing links per individual and return long format
summarise_links <- function(df) {
    summary_df <- df %>%
        group_by(SUBJID) %>%
        summarise(enhancing = sum(Mechanism == 'enhancing'),
                  repressive = sum(Mechanism == 'repressing'),
                  total = n(),
                  AGE_collapsed = first(AGE_collapsed),
                  SEX = first(SEX)) %>%
        ungroup()
    summary_df$AGE_collapsed <- factor(summary_df$AGE_collapsed, levels = c('20-49', '50-59', '60-79'))
    summary_df$SEX <- factor(summary_df$SEX, levels = c('Male', 'Female'))
    summary_df <- pivot_longer(summary_df, cols = c(enhancing, repressive), names_to = "LinkType", values_to = "Count")
    return(data.frame(summary_df))
}

# Function to create link information
link_key <- function(df) paste(df$ncRNA, df$Name_pcGene, df$Mechanism, sep=" | ")

# Function to calculate Jaccard index between two sets
jaccard_index <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    if (union == 0) return(0)
    return(intersection / union)
}

# Function to create matrix of Jaccard indices
create_jaccard_matrix <- function(df) {
    link_vector <- lapply(unique(df$AGE_collapsed), function(age_group) {
        subset_df <- subset(df, AGE_collapsed == age_group)
        unique(link_key(subset_df))
    })
    names(link_vector) <- unique(df$AGE_collapsed)
    jaccard_matrix <- matrix(0, nrow = length(link_vector), ncol = length(link_vector))
    rownames(jaccard_matrix) <- c('20-49', '50-59', '60-79')
    colnames(jaccard_matrix) <- c('20-49', '50-59', '60-79')
    for (i in names(link_vector)) {
        for (j in names(link_vector)) {
            jaccard_matrix[i, j] <- jaccard_index(link_vector[[i]], link_vector[[j]])
        }
    }
    return(jaccard_matrix)
}

################################################
# Show the sample size across tissues and ages #
################################################
sample_sizes <- gtex_link %>%
    select(SUBJID, Tissue, AGE_collapsed) %>%
    distinct() %>%
    group_by(Tissue, AGE_collapsed) %>%
    summarise(SampleSize = n(), .groups = 'drop')

pdf('Sample_Sizes_by_Tissue_and_Age.pdf')
ggplot(sample_sizes, aes(x = AGE_collapsed, y = Tissue, size = SampleSize)) +
    geom_point() +
    scale_size_area(max_size = 4) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
dev.off()


###########################################################
# Plot bar plots of enhancing and repressing links by age #
###########################################################

# Subset for heart
left_ventricle <- subset(gtex_link, Tissue == 'Heart - Left Ventricle')
atrial_appendage <- subset(gtex_link, Tissue == 'Heart - Atrial Appendage')


pdf('LV_links_by_age.pdf')
left_ventricle_summary <- summarise_links(left_ventricle)
ggplot(left_ventricle_summary, aes(x = AGE_collapsed, y = Count, fill = LinkType)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~SEX) +
  labs(x = "Age Group", y = "Number of Links", fill = "Link Type") +
  theme_minimal()
dev.off()


pdf('AA_links_by_age.pdf')
atrial_appendage_summary <- summarise_links(atrial_appendage)
ggplot(atrial_appendage_summary, aes(x = AGE_collapsed, y = Count, fill = LinkType)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~SEX) +
  labs(x = "Age Group", y = "Number of Links", fill = "Link Type") +
  theme_minimal()
dev.off()

pdf('All_links_by_age.pdf')
all_organs_summary <- summarise_links(gtex_link)
ggplot(all_organs_summary, aes(x = AGE_collapsed, y = Count, fill = LinkType)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~SEX) +
  labs(x = "Age Group", y = "Number of Links", fill = "Link Type") +
  theme_minimal()
dev.off()


##################
# Jaccard matrix #
##################

# Calculate Jaccard indices for LV
lv_jaccard_matrix <- create_jaccard_matrix(left_ventricle)

# Plot heatmap for the upper triangle of the Jaccard matrix
lv_jaccard_matrix[lower.tri(lv_jaccard_matrix)] <- NA
col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "blue", "red"))
pdf('LV_Jaccard_Heatmap.pdf')
Heatmap(lv_jaccard_matrix, name = "Jaccard Index", col = col_fun, na_col = "white",
        cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(lv_jaccard_matrix[i, j])) {
                grid.text(sprintf("%.2f", lv_jaccard_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
            }
        })
dev.off()


# Calculate Jaccard indices for AA
aa_jaccard_matrix <- create_jaccard_matrix(atrial_appendage)
# Plot heatmap for the upper triangle of the Jaccard matrix
aa_jaccard_matrix[lower.tri(aa_jaccard_matrix)] <- NA
pdf('AA_Jaccard_Heatmap.pdf')
Heatmap(aa_jaccard_matrix, name = "Jaccard Index", col = col_fun, na_col = "white",
        cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(aa_jaccard_matrix[i, j])) {
                grid.text(sprintf("%.2f", aa_jaccard_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
            }
        })
dev.off()

# Calculate Jaccard indices for all organs
all_jaccard_matrix <- create_jaccard_matrix(gtex_link)
# Plot heatmap for the upper triangle of the Jaccard matrix
all_jaccard_matrix[lower.tri(all_jaccard_matrix)] <- NA
pdf('All_Jaccard_Heatmap.pdf')
Heatmap(all_jaccard_matrix, name = "Jaccard Index", col = col_fun, na_col = "white",
        cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(all_jaccard_matrix[i, j])) {
                grid.text(sprintf("%.2f", all_jaccard_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
            }
        })
dev.off()


#################################
# How many switch during aging? #
#################################

pair_key <- function(df) paste(df$ncRNA, df$Name_pcGene, sep="|")

lv_vector <- lapply(split(left_ventricle, left_ventricle$AGE_collapsed), function(x){
    x %>% 
    mutate(pair = pair_key(.), mech = Mechanism) %>%
    select(pair, mech) %>%
    unique()
})
age_names <- names(lv_vector)
lv_merged <- lv_vector[[1]]
for (i in 2:length(lv_vector)) {
  lv_merged <- full_join(
    lv_merged,
    lv_vector[[i]],
    by = "pair",
    suffix = c(paste0("_", age_names[i-1]), paste0("_", age_names[i])),
    relationship = "many-to-many"
  )
}

colnames(lv_merged) <- gsub('mech_', '', colnames(lv_merged))
lv_merged[,-1] <- apply(lv_merged[,-1], 2, function(x) ifelse(x == 'enhancing', 1, -1))
lv_merged[,-1][is.na(lv_merged[,-1])] <- 0

lv_merged_long <- pivot_longer(lv_merged, cols = -pair, names_to = "Age", values_to = "LinkType")
lv_merged_long$Age <- factor(lv_merged_long$Age, levels = age_names)

# Line plot of link type changes
ggplot(lv_merged_long, aes(x = Age, y = LinkType)) +
  geom_line(alpha = 0.1) +
  geom_point(alpha = 0.1) +
  scale_y_continuous(breaks = c(-1, 0, 1), labels = c("Repressing", "Absent", "Enhancing")) +
  labs(x = "Age Group", y = "Link Type") +
  theme_minimal()
ggsave('LV_link_switching_lineplot.pdf')

col_fun <- c("enhancing" = "blue", "repressing" = "red", 'NA' = "grey")
pdf('LV_link_switching.pdf')
Heatmap(
  as.matrix(lv_merged[,-1]),
  name = "Link Type",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE
)
dev.off()

#############################################
# Alluvial plot to show changes across ages #
##############################################
age_levels <- c("20-49", "50-59", "60-79")
status_by_age <- left_ventricle %>%
  group_by(AGE_collapsed, Link) %>%
  summarise(
    n_enh = sum(Mechanism == "enhancing"),
    n_rep = sum(Mechanism == "repressing"),
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      n_enh == 0 & n_rep == 0 ~ "No link",
      n_enh >  n_rep          ~ "Enhancing",
      n_rep >  n_enh          ~ "Repressing",
      TRUE                    ~ "Mixed"     # tie within the bin
    ),
    status = factor(status, levels = c("No link","Repressing","Enhancing","Mixed"))
  ) %>%
  complete(AGE_collapsed = age_levels, Link,
           fill = list(n_enh = 0, n_rep = 0,
                       status = factor("No link",
                         levels = c("No link","Repressing","Enhancing","Mixed"))))

## 3) Lodes for ggalluvial: one row per (Link, AGE)
lodes <- status_by_age %>%
  arrange(Link, AGE_collapsed) %>%
  mutate(alluvium = Link) %>%
  transmute(x = AGE_collapsed, alluvium, stratum = status, y = 1)

## 4) Plot across ages (no sex facets)
pdf('LV_link_switching_alluvial.pdf')
ggplot(lodes,
       aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = stratum)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(colour = "grey40", fill = "grey92") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  labs(x = "", y = "Number of links", fill = "Status") +
  theme_minimal() +
  theme(panel.grid = element_blank())
dev.off()



# Repeat alluvial for AA
status_by_age_aa <- atrial_appendage %>%
  group_by(AGE_collapsed, Link) %>%
  summarise(
    n_enh = sum(Mechanism == "enhancing"),
    n_rep = sum(Mechanism == "repressing"),
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      n_enh == 0 & n_rep == 0 ~ "No link",
      n_enh >  n_rep          ~ "Enhancing",
      n_rep >  n_enh          ~ "Repressing",
      TRUE                    ~ "Mixed"     # tie within the bin
    ),
    status = factor(status, levels = c("No link","Repressing","Enhancing","Mixed"))
  ) %>%
  complete(AGE_collapsed = age_levels, Link,
           fill = list(n_enh = 0, n_rep = 0,
                       status = factor("No link",
                         levels = c("No link","Repressing","Enhancing","Mixed"))))
## 3) Lodes for ggalluvial: one row per (Link, AGE)
lodes_aa <- status_by_age_aa %>%
  arrange(Link, AGE_collapsed) %>%
  mutate(alluvium = Link) %>%
  transmute(x = AGE_collapsed, alluvium, stratum = status, y = 1) 
## 4) Plot
pdf('AA_link_switching_alluvial.pdf', width = 8, height = 5)
ggplot(lodes_aa,
       aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = stratum)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(colour = "grey40", fill = "grey92") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  labs(x = "", y = "Number of links", fill = "Status") +
  theme_minimal() +
  theme(panel.grid = element_blank())
dev.off()



# Repeat alluvial for all organs
status_by_age_all <- gtex_link %>%
  group_by(AGE, Link) %>%
  summarise(
    n_enh = sum(Mechanism == "enhancing"),
    n_rep = sum(Mechanism == "repressing"),
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      n_enh == 0 & n_rep == 0 ~ "No link",
      n_enh >  n_rep          ~ "Enhancing",
      n_rep >  n_enh          ~ "Repressing",
      TRUE                    ~ "Mixed"     # tie within the bin
    ),
    status = factor(status, levels = c("No link","Repressing","Enhancing","Mixed"))
  ) %>%
  complete(AGE = age_levels, Link,
           fill = list(n_enh = 0, n_rep = 0,
                       status = factor("No link",
                         levels = c("No link","Repressing","Enhancing","Mixed"))))
## 3) Lodes for ggalluvial: one row per (Link, AGE)
lodes_all <- status_by_age_all %>%
  arrange(Link, AGE) %>%
  mutate(alluvium = Link) %>%
  transmute(x = AGE, alluvium, stratum = status, y = 1) 
## 4) Plot
ggplot(lodes_all,
       aes(x = x, stratum = stratum, alluvium = alluvium, y = y, fill = stratum)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(colour = "grey40", fill = "grey92") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  labs(x = "", y = "Number of links", fill = "Status") +
  theme_minimal() +
  theme(panel.grid = element_blank())
ggsave('All_link_switching_alluvial.pdf', width = 8, height = 5)

tmp <- unique(gtex_link[,c('SUBJID', 'AGE')])
table(tmp$AGE)

###################################
# Link accumulation curves by age #
###################################

library(dplyr)
library(purrr)
library(ggplot2)
set.seed(42)

for (tissue in unique(gtex_link$Tissue)) {

    lv <- gtex_link %>%
    filter(Tissue == tissue) %>%
    mutate(Key = link_key(.)) %>%
    distinct(SUBJID, AGE_collapsed, Key)

    # list of links per subject within each age
    sbj_links <- lv %>%
    group_by(AGE_collapsed, SUBJID) %>%
    summarise(links = list(Key), .groups = "drop")

    #--- 2) Function: accumulation curve by random subject order -------------------
    acc_curve <- function(link_lists, n_perm = 200) {
    # link_lists: list of character vectors, one per subject
    S <- length(link_lists)
    if (S == 0) return(NULL)

    # pre-allocate (rows: k=1..S, cols: permutations)
    out <- matrix(NA_integer_, nrow = S, ncol = n_perm)

    for (p in seq_len(n_perm)) {
        ord <- sample.int(S)
        seen <- character(0)
        for (k in seq_len(S)) {
        seen <- union(seen, link_lists[[ord[k]]])
        out[k, p] <- length(seen)
        }
    }

    # summarise across permutations for each k (no. of individuals)
    tibble(
        n_individuals = seq_len(S),
        mean_links    = rowMeans(out),
        lwr_links     = apply(out, 1, quantile, probs = 0.05),
        upr_links     = apply(out, 1, quantile, probs = 0.95)
    )
    }

    #--- 3) Run curves per age group ----------------------------------------------
    curves <- sbj_links %>%
    group_by(AGE_collapsed) %>%
    summarise(
        curve = list(acc_curve(links)),
        n_subjects = n_distinct(SUBJID),
        .groups = "drop"
    ) %>%
    tidyr::unnest(curve)

    total_links_lv <- n_distinct(lv$Key)

    #--- 4) Plot -------------------------------------------------------------------
    p <- ggplot(curves, aes(x = n_individuals, y = mean_links, colour = AGE_collapsed, fill = AGE_collapsed)) +
    geom_line(linewidth = 0.9) +
    labs(
        title = paste0(tissue, " (Total unique links n = ", total_links_lv, ")"),
        x = "No. of individuals",
        y = "No. of linkages"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.title = element_blank())
    print(p)

    # Optionally save
    ggsave(paste0('accumulation_plot/', tissue, "_link_accumulation_by_age.pdf"), width = 6.5, height = 5.5)
}

# # Comparing AUC of accumulation lines between age groups - subset to tissues with smallest sample sizes

# library(dplyr)
# library(purrr)
# library(ggplot2)
# library(tidyr)
# library(DescTools)
# set.seed(42)

# #--- Modified accumulation function: return all permutations -------------------
# acc_curve_perm <- function(link_lists, n_perm = 200) {
#   S <- length(link_lists)
#   if (S == 0) return(NULL)

#   out <- vector("list", n_perm)

#   for (p in seq_len(n_perm)) {
#     ord <- sample.int(S)
#     seen <- character(0)
#     counts <- integer(S)
#     for (k in seq_len(S)) {
#       seen <- union(seen, link_lists[[ord[k]]])
#       counts[k] <- length(seen)
#     }
#     out[[p]] <- tibble(
#       perm = p,
#       n_individuals = seq_len(S),
#       n_links = counts
#     )
#   }
#   bind_rows(out)
# }

# #--- Run per tissue ------------------------------------------------------------
# results <- list()

# for (tissue in unique(gtex_link$Tissue)) {

#   lv <- gtex_link %>%
#     filter(Tissue == tissue) %>%
#     mutate(Key = link_key(.)) %>%
#     distinct(SUBJID, AGE_collapsed, Key)

#   sbj_links <- lv %>%
#     group_by(AGE_collapsed, SUBJID) %>%
#     summarise(links = list(Key), .groups = "drop")

#   curves <- sbj_links %>%
#     group_by(AGE_collapsed) %>%
#     summarise(curve = list(acc_curve_perm(links)), .groups = "drop") %>%
#     unnest(curve)

#   # truncate at min sample size across groups
#   max_depth <- curves %>%
#     group_by(AGE_collapsed, perm) %>%
#     summarise(max_n = max(n_individuals), .groups = "drop") %>%
#     summarise(min_n = min(max_n)) %>%
#     pull(min_n)

#   curves_trunc <- curves %>%
#     filter(n_individuals <= max_depth)

#   # compute AUC per permutation
#   auc_df <- curves_trunc %>%
#     group_by(AGE_collapsed, perm) %>%
#     summarise(
#       auc = AUC(n_individuals, n_links, method = "trapezoid"),
#       .groups = "drop"
#     ) %>%
#     mutate(Tissue = tissue)

#   results[[tissue]] <- auc_df
# }

# auc_all <- bind_rows(results)

# #--- Compare groups ------------------------------------------------------------
# pdf('AUC_by_age_group.pdf', width = 8, height = 6)
# ggplot(auc_all, aes(x = Tissue, y = auc, fill = AGE_collapsed)) +
#   geom_boxplot(alpha = 0.6) +
#   theme_minimal(base_size = 12) +
#   theme(legend.position = "none") +
#   labs(y = "AUC (truncated)", x = "Tissue")
# dev.off()

# # Example statistical test (within one tissue)
# library(rstatix)
# test_results <- lapply(unique(auc_all$Tissue), function(tissue) {
#   tissue_df <- subset(auc_all, Tissue == tissue)
#   pairwise_wilcox_test(auc ~ AGE_collapsed, data = tissue_df, p.adjust.method = "fdr")
# })

##############################################################
# Bootstrap to test for age-related changes in link presence #
##############################################################
# link_presence <- function(df) {
#   df %>%
#     mutate(key = link_key(.), present = 1) %>%
#     distinct(SUBJID, AGE, key, present) %>%
#     tidyr::complete(tidyr::nesting(SUBJID, AGE), key,
#                     fill = list(present = 0))
# }

# tissue <- "Heart - Atrial Appendage"
# tmp <- subset(gtex_link, Tissue == tissue)

# library(future.apply)
# plan(multicore)
# future.seed = NULL
# perm_result <- future_lapply(1:1000, function(i) {
#     # Find age group with minimum number of subjects
#     min_age_group <- tmp %>%
#     distinct(SUBJID, AGE) %>%
#     group_by(AGE) %>%
#     summarise(n = n()) %>%
#     arrange(n) %>%
#     slice(1)
#     min_n <- min_age_group$n

#     sampled_ids <- tmp %>%
#     group_by(AGE) %>%
#     distinct(SUBJID) %>%
#     slice_sample(n = min_n) %>%
#     ungroup()
#     sampled_ids <- sampled_ids$SUBJID

#     tmp_subset <- tmp %>%
#     filter(SUBJID %in% sampled_ids)

#     lnc_binary <- link_presence(tmp_subset)

#     link_age_tests <- lnc_binary %>%
#         group_by(key) %>%
#         do({
#             fit <- glm(present ~ AGE, data = ., family = binomial)
#             em <- emmeans(fit, ~ AGE, type = "response") # gives predicted probs
#             tidy(pairs(em, adjust = "fdr"))
#         }) %>% 
#         ungroup()
#     return(link_age_tests)
# })

# unlist(lapply(perm_result, function(df) {
#     nrow(subset(df, adj.p.value < 0.05))
# }))


##################################################
# Modelling age-related changes in link presence #
#########################################
link_presence_mechanism <- function(df) {
  df %>%
    mutate(key = link_key(.)) %>%
    mutate(present = 1) %>%
    select(SUBJID, AGE_collapsed, SEX, key, present) %>%
    distinct() %>%
    complete(SUBJID, AGE_collapsed, SEX, key, fill = list(present = 0))
}

link_presence <- function(df) {
    df %>%
    mutate(present = 1) %>%
    select(SUBJID, AGE_collapsed, SEX, Link, present) %>%
    distinct() %>%
    complete(SUBJID, AGE_collapsed, SEX, Link, fill = list(present = 0))
}

age_change_links <- lapply(unique(gtex_link$Tissue), function(tissue){
    tissue_df <- subset(gtex_link, Tissue == tissue)
    lnc_binary <- link_presence_mechanism(tissue_df)

    link_age_tests <- lnc_binary %>%
        group_by(key) %>%
        do({
            fit <- glm(present ~ AGE_collapsed, data = ., family = binomial)
            em <- emmeans(fit, ~ AGE_collapsed, type = "response")
            tidy(pairs(em, adjust = "fdr"))
        }) %>%
        ungroup() %>%
        data.frame()
    return(link_age_tests)
})
names(age_change_links) <- unique(gtex_link$Tissue)

# Matrix of significant changes
sig_links <- lapply(age_change_links, function(df){
    unique(subset(df, adj.p.value < 0.05)$key)
})

sig_links[sapply(sig_links, function(x) length(x) > 0)]

links_mtx <- fromList(sig_links)
colnames(links_mtx) <- names(sig_links)
rownames(links_mtx) <- unique(unlist(sig_links))

pdf('Link_Age_Changes_collapsed_heatmap.pdf')
Heatmap(t(links_mtx), col = c("white", "red"),
        show_heatmap_legend = FALSE,
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8))
dev.off()


# Plot effects for significant links
lapply(names(sig_links), function(tissue) {
    if(sig_links[[tissue]] %>% length() > 0) {
        for(i in 1:length(sig_links[[tissue]])) {
            link <- sig_links[[tissue]][i]
            tissue_df <- subset(gtex_link, Tissue == tissue)
            lnc_binary <- link_presence_collapsed(tissue_df)
            sig_df <- lnc_binary %>% filter(key %in% sig_links[[tissue]])
            fit <- glm(present ~ AGE_collapsed, data = sig_df, family = binomial)
            em <- emmeans(fit, ~ AGE_collapsed, type = "response")
            em_df <- as.data.frame(em)
            pdf(paste0('link_age_effects/', tissue, '_', i, '_age_effects.pdf'), width = 8, height = 6)
            p <- ggplot(em_df, aes(x = AGE_collapsed, y = prob)) +
                geom_point(size = 3) +
                geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
                labs(x = "Age Group", y = "Predicted Probability", title = paste(tissue, ': ', link)) +
                theme_minimal()
            print(p)
            dev.off()
        }
    }
})


# Clustered bar plot of significant links to show mechanism shifting over age
sig_results <- lapply(age_change_links, function(df){
    subset(df, adj.p.value < 0.05)[, c('key', 'term')]
})
sig_results <- unique(bind_rows(sig_results, .id = 'Tissue'))

lapply(1:nrow(sig_results), function(i){
    tissue <- sig_results$Tissue[i]
    key <- sig_results$key[i]
    lnc <- strsplit(key, ' \\| ')[[1]][1]
    pcgene <- strsplit(key, ' \\| ')[[1]][2]
    link <- paste(lnc, pcgene, sep = '_')
    tmp <- subset(gtex_link, Tissue == tissue & Link %in% link)
    tmp_summary <- tmp %>%
        group_by(AGE_collapsed, Mechanism) %>%
        summarise(Count = n(), .groups = 'drop')
    pdf(paste0(tissue, '_', link, '_link_counts_by_age.pdf'))
    p <- ggplot(tmp_summary, aes(x = AGE_collapsed, y = Count, fill = Mechanism)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "Age Group", y = "Number of Links", fill = "Link Type",
             title = paste(tissue, ': ', key)) +
        theme_minimal()
    print(p)
    dev.off()
})




###############################################################
# Compare expression of lncRNAs and pcGenes across age groups #
################################################################

gtex_TPM <- read.delim('GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz', skip = 2)
colnames(gtex_TPM)[3:ncol(gtex_TPM)] <- gsub('\\.', '-', colnames(gtex_TPM)[3:ncol(gtex_TPM)])
gtex_TPM_metadata <- read.delim('GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt')
gtex_TPM_metadata$SUBJID <- unlist(lapply(gtex_TPM_metadata$SAMPID, function(x) {
    paste(strsplit(x, '-')[[1]][1], strsplit(x, '-')[[1]][2], sep = '-')
}))
gtex_TPM_metadata <- merge(gtex_TPM_metadata, metadata[, c("SUBJID", "AGE", "SEX")], by = "SUBJID")
gtex_TPM_metadata$AGE_collapsed <- ifelse(gtex_TPM_metadata$AGE %in% c('20-29', '30-39', '40-49'), '20-49',
                                              ifelse(gtex_TPM_metadata$AGE %in% c('60-69', '70-79'), '60-79', '50-59'))
gtex_TPM_metadata$AGE_collapsed <- factor(gtex_TPM_metadata$AGE_collapsed, levels = c('20-49', '50-59', '60-79'))
gtex_TPM_metadata <- subset(gtex_TPM_metadata, SAMPID %in% colnames(gtex_TPM)[-c(1,2)])

sig_links_flt <- sig_links[sapply(sig_links, function(x) length(x) > 0)]
lapply(names(sig_links_flt), function(tissue) {
    x <- sig_links_flt[[tissue]]
    if(length(x) > 0) {
        for(i in 1:length(x)) {
            link <- x[i]
            lnc <- strsplit(link, ' \\| ')[[1]][1]
            pcgene <- strsplit(link, ' \\| ')[[1]][2]

            expr <- gtex_TPM[,subset(gtex_TPM_metadata, SMTSD == tissue)$SAMPID]
            expr <- cbind(gtex_TPM[,1:2], expr)

            lnc_test <- subset(expr, Description == lnc) %>% pivot_longer(cols = -c(Name, Description), names_to = "SAMPID", values_to = "TPM")
            lnc_test <- merge(lnc_test, gtex_TPM_metadata[, c("SAMPID", "AGE_collapsed")], by = "SAMPID")
            tmp <- pairwise_wilcox_test(TPM ~ AGE_collapsed, data = lnc_test, p.adjust.method = "fdr")
            tmp$gene <- lnc
            tmp$tissue <- tissue
            print(data.frame(tmp))

            pcgene_test <- subset(expr, Description == pcgene) %>% pivot_longer(cols = -c(Name, Description), names_to = "SAMPID", values_to = "TPM")
            pcgene_test <- merge(pcgene_test, gtex_TPM_metadata[, c("SAMPID", "AGE_collapsed")], by = "SAMPID")
            tmp <- pairwise_wilcox_test(TPM ~ AGE_collapsed, data = pcgene_test, p.adjust.method = "fdr")
            tmp$gene <- pcgene
            tmp$tissue <- tissue
            print(data.frame(tmp))

            plt_data <- subset(expr, Description %in% c(lnc, pcgene)) %>%
                pivot_longer(cols = -c(Name, Description), names_to = "SAMPID", values_to = "TPM")
            plt_data <- merge(plt_data, gtex_TPM_metadata[, c("SAMPID", "AGE_collapsed")], by = "SAMPID")

            
            pdf(paste0('link_expression_by_age/', tissue, '_', i, '_expression_by_age.pdf'))
                  p <- ggplot(plt_data, aes(x = AGE_collapsed, y = TPM, fill = AGE_collapsed)) +
                    geom_boxplot() +
                    facet_wrap(~Description, scales = 'free_y') +
                    theme_minimal() +
                    labs(x = "Age Group", y = "TPM", title = paste(tissue, ': ', link)) +
                    theme(legend.position = "none") +
                    stat_compare_means(method = "wilcox.test", label = "p.signif", 
                              comparisons = combn(levels(plt_data$AGE_collapsed), 2, simplify = FALSE))
                  print(p)
                  dev.off()
        }
    }
})