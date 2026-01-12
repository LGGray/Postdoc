library(NMF)
library(mulea)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Extract top genes function
top_genes <- function(W, comp, n=50){
    ord <- order(W[, comp], decreasing = TRUE)
    data.frame(
        gene = rownames(W)[ord][1:n],
        loading = W[ord, comp][1:n]
    )
}

# Set up pathway data
msigdb_BP <- read_gmt('pathways/m5.go.bp.v2025.1.Mm.symbols.gmt')
# Filtering ontologies with too few or too many genes
msigdb_BP_filtered <- filter_ontology(gmt = msigdb_BP,
                                        min_nr_of_elements = 3,
                                        max_nr_of_elements = 400)
msigsb_hallmark <- read_gmt('pathways/mh.all.v2025.1.Mm.symbols.gmt')
msigsb_hallmark_filtered <- filter_ontology(gmt = msigsb_hallmark,
                                        min_nr_of_elements = 3,
                                        max_nr_of_elements = 400)


# Set Seed
set.seed(42)

# Adult vs Aged
adult_vs_age <- read.delim("adult_aged_bodymap/DEG/He_GEM.txt")
adult_vs_age <- adult_vs_age[rowSums(adult_vs_age) > 0, ]

res <- nmf(adult_vs_age, 2:6, nrun = 10)

pdf('adult_aged_bodymap/NMF/NMF_rank.pdf')
plot(res)
dev.off()

adult_age_fit3 <- nmf(adult_vs_age, 3, nrun=50)

pdf('adult_aged_bodymap/NMF/NMF_basis_heatmap_rank3.pdf', onefile=FALSE)
basismap(adult_age_fit3, tracks=':basis', annColor=list(basis=1:3))
dev.off()

adult_age_W3 <- basis(adult_age_fit3)

# Extract top genes for each component
comp1_adult_age <- top_genes(adult_age_W3, 1, 50)
comp2_adult_age <- top_genes(adult_age_W3, 2, 50)
comp3_adult_age <- top_genes(adult_age_W3, 3, 50)


# Over rerepresentation analysis
ora_comp1_adult_age <- ora(
  gmt = msigdb_BP_filtered,
  element_names = comp1_adult_age$gene,                 
  background_element_names = rownames(adult_vs_age),
  p_value_adjustment_method = "eFDR",      
  number_of_permutations = 1000,         
  random_seed = 42
)
ora_results_comp1_adult_age <- run_test(ora_comp1_adult_age)

ora_comp2_adult_age <- ora(
  gmt = msigdb_BP_filtered,
  element_names = comp2_adult_age$gene,                 
  background_element_names = rownames(adult_vs_age),
  p_value_adjustment_method = "eFDR",      
  number_of_permutations = 1000,        
  random_seed = 42
)
ora_results_comp2_adult_age <- run_test(ora_comp2_adult_age)

ora_comp3_adult_age <- ora(
  gmt = msigdb_BP_filtered,
  element_names = comp3_adult_age$gene,                
  background_element_names = rownames(adult_vs_age),
  p_value_adjustment_method = "eFDR",      
  number_of_permutations = 1000,        
  random_seed = 42
)
ora_results_comp3_adult_age <- run_test(ora_comp3_adult_age)

component_mtx <- as.matrix(fromList(list(Comp1=subset(ora_results_comp1_adult_age, eFDR < 0.05)$ontology_id,
              Comp2=subset(ora_results_comp2_adult_age, eFDR < 0.05)$ontology_id,
              Comp3=subset(ora_results_comp3_adult_age, eFDR < 0.05)$ontology_id)))
pdf('adult_aged_bodymap/NMF/heatmap_NMF_pathway_overlap.pdf')
Heatmap(component_mtx, name='Pathway Overlap', 
        col = colorRamp2(c(0, 1), c("white", "blue")),
        cluster_rows = TRUE,
        cluster_columns = TRUE
)
dev.off()

comp1_unique <- setdiff(subset(ora_results_comp1_adult_age, eFDR < 0.01)$ontology_id, c(subset(ora_results_comp2_adult_age, eFDR < 0.01)$ontology_id, subset(ora_results_comp3_adult_age, eFDR < 0.01)$ontology_id))
comp2_unique <- setdiff(subset(ora_results_comp2_adult_age, eFDR < 0.01)$ontology_id, c(subset(ora_results_comp1_adult_age, eFDR < 0.01)$ontology_id, subset(ora_results_comp3_adult_age, eFDR < 0.01)$ontology_id))
comp3_unique <- setdiff(subset(ora_results_comp3_adult_age, eFDR < 0.01)$ontology_id, c(subset(ora_results_comp1_adult_age, eFDR < 0.01)$ontology_id, subset(ora_results_comp2_adult_age, eFDR < 0.01)$ontology_id))

unique_pathways <- bind_rows(
  subset(ora_results_comp1_adult_age, ontology_id %in% comp1_unique),
  subset(ora_results_comp2_adult_age, ontology_id %in% comp2_unique),
  subset(ora_results_comp3_adult_age, ontology_id %in% comp3_unique)
, .id = 'component'
)
unique_pathways$gene_ratio <- unique_pathways$nr_common_with_tested_elements / unique_pathways$nr_common_with_background_elements
unique_pathways$ontology_id <- gsub('GOBP_', '', unique_pathways$ontology_id)

ht <- Heatmap(
  matrix = as.matrix(
    unique_pathways %>%
      select(component, ontology_id, gene_ratio) %>%
      pivot_wider(names_from = component, values_from = gene_ratio, values_fill = 0) %>%
      column_to_rownames('ontology_id')
  ),
  name = 'Gene Ratio',
  col = colorRamp2(c(0, max(unique_pathways$gene_ratio)), c("white", "blue")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 12, rot = 0),
  column_title = NULL,
  width = unit(4, "cm")
)

pdf('adult_aged_bodymap/NMF/heatmap_NMF_unique_pathways.pdf', width = 15, height = 8)
draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()
# Sham vs TAC #
###############
sham_vs_TAC <- read.delim('F1_TAC_Sarah/DEG/TAC_GEM.txt')
sham_vs_TAC <- sham_vs_TAC[rowSums(sham_vs_TAC) > 0, ]
res <- nmf(sham_vs_TAC, 2:6, nrun = 10)

pdf('F1_TAC_Sarah/NMF/NMF_rank_sham_TAC.pdf')
plot(res)
dev.off()

sham_TAC_fit3 <- nmf(sham_vs_TAC, 3, nrun=50)

pdf('F1_TAC_Sarah/NMF/NMF_basis_heatmap_rank3_sham_TAC.pdf', onefile=FALSE)
basismap(sham_TAC_fit3, tracks=':basis', annColor=list(basis=1:3))
dev.off()

sham_TAC_W3 <- basis(sham_TAC_fit3)
comp1_sham_TAC <- top_genes(sham_TAC_W3, 1, 50)
comp2_sham_TAC <- top_genes(sham_TAC_W3, 2, 50)
comp3_sham_TAC <- top_genes(sham_TAC_W3, 3, 50)

ora_comp1_sham_TAC <- ora(
  gmt = msigdb_BP_filtered,
  element_names = comp1_sham_TAC$gene,                 
  background_element_names = rownames(sham_vs_TAC),
  p_value_adjustment_method = "eFDR",      
  number_of_permutations = 1000,         
  random_seed = 42
)
ora_results_comp1_sham_TAC <- run_test(ora_comp1_sham_TAC)

ora_comp2_sham_TAC <- ora(
  gmt = msigdb_BP_filtered,
  element_names = comp2_sham_TAC$gene,                 
  background_element_names = rownames(sham_vs_TAC),
  p_value_adjustment_method = "eFDR",      
  number_of_permutations = 1000,        
  random_seed = 42
)
ora_results_comp2_sham_TAC <- run_test(ora_comp2_sham_TAC)

ora_comp3_sham_TAC <- ora(
  gmt = msigdb_BP_filtered,
  element_names = comp3_sham_TAC$gene,                 
  background_element_names = rownames(sham_vs_TAC),
  p_value_adjustment_method = "eFDR",      
  number_of_permutations = 1000,        
  random_seed = 42
)
ora_results_comp3_sham_TAC <- run_test(ora_comp3_sham_TAC)

comp1_unique <- setdiff(subset(ora_results_comp1_sham_TAC, eFDR < 0.01)$ontology_id, c(subset(ora_results_comp2_sham_TAC, eFDR < 0.01)$ontology_id, subset(ora_results_comp3_sham_TAC, eFDR < 0.01)$ontology_id))
comp2_unique <- setdiff(subset(ora_results_comp2_sham_TAC, eFDR < 0.01)$ontology_id, c(subset(ora_results_comp1_sham_TAC, eFDR < 0.01)$ontology_id, subset(ora_results_comp3_sham_TAC, eFDR < 0.01)$ontology_id))
comp3_unique <- setdiff(subset(ora_results_comp3_sham_TAC, eFDR < 0.01)$ontology_id, c(subset(ora_results_comp1_sham_TAC, eFDR < 0.01)$ontology_id, subset(ora_results_comp2_sham_TAC, eFDR < 0.01)$ontology_id))

unique_pathways <- bind_rows(
  subset(ora_results_comp1_sham_TAC, ontology_id %in% comp1_unique),
  subset(ora_results_comp2_sham_TAC, ontology_id %in% comp2_unique),
  subset(ora_results_comp3_sham_TAC, ontology_id %in% comp3_unique)
, .id = 'component'
)
unique_pathways$gene_ratio <- unique_pathways$nr_common_with_tested_elements / unique_pathways$nr_common_with_background_elements
unique_pathways$ontology_id <- gsub('GOBP_', '', unique_pathways$ontology_id)

ht <- Heatmap(
  matrix = as.matrix(
    unique_pathways %>%
      select(component, ontology_id, gene_ratio) %>%
      pivot_wider(names_from = component, values_from = gene_ratio, values_fill = 0) %>%
      column_to_rownames('ontology_id')
  ),
  name = 'Gene Ratio',
  col = colorRamp2(c(0, max(unique_pathways$gene_ratio)), c("white", "blue")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 12, rot = 0),
  column_title = NULL,
  width = unit(4, "cm")
)

pdf('F1_TAC_Sarah/NMF/heatmap_NMF_unique_pathways.pdf', width = 15, height = 8)
draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()






# Compare genes in each component between adult vs aged and sham vs TAC
gene_jaccard <- sapply(list(comp1_adult_age=comp1_adult_age$gene, comp2_adult_age=comp2_adult_age$gene, comp3_adult_age=comp3_adult_age$gene), function(g1)
sapply(list(comp1_sham_TAC=comp1_sham_TAC$gene, comp2_sham_TAC=comp2_sham_TAC$gene, comp3_sham_TAC=comp3_sham_TAC$gene), function(g2)
    length(intersect(g1, g2)) / length(union(g1, g2))
)
)

pdf('LINKS_study/NMF_gene_jaccard.pdf')
Heatmap(gene_jaccard, name='Jaccard Index', 
        col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE
)
dev.off()


# Compare pathways in each component between adult vs aged and sham vs TAC
pathway_jaccard <- sapply(list(ora_results_comp1_adult_age=subset(ora_results_comp1_adult_age, eFDR < 0.05)$ontology_id,
            ora_results_comp2_adult_age=subset(ora_results_comp2_adult_age, eFDR < 0.05)$ontology_id,
            ora_results_comp3_adult_age=subset(ora_results_comp3_adult_age, eFDR < 0.05)$ontology_id), function(p1)
sapply(list(ora_results_comp1_sham_TAC=subset(ora_results_comp1_sham_TAC, eFDR < 0.05)$ontology_id,
            ora_results_comp2_sham_TAC=subset(ora_results_comp2_sham_TAC, eFDR < 0.05)$ontology_id,
            ora_results_comp3_sham_TAC=subset(ora_results_comp3_sham_TAC, eFDR < 0.05)$ontology_id), function(p2)
    length(intersect(p1, p2)) / length(union(p1, p2))
)
)

pdf('LINKS_study/NMF_pathway_jaccard.pdf')
Heatmap(pathway_jaccard, name='Jaccard Index', 
        col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE
)
dev.off()


# Overlap links identifed 
load('LINKS_study/adult_aged_TAC_unique_link_overlaps.RData')
base <- unlist(lapply(unlist(overlap_lists),function(x) strsplit(x, '\\|')[[1]][1]))
target <- unlist(lapply(unlist(overlap_lists),function(x) strsplit(x, '\\|')[[1]][2]))

rank_adult <- apply(adult_age_W3, 2, function(v) rank(-v))   # high loading → rank 1
rank_TAC   <- apply(sham_TAC_W3, 2, function(v) rank(-v))
zs_adult <- scale(adult_age_W3)
zs_TAC   <- scale(sham_TAC_W3)

rank_adult[rownames(rank_adult) %in% c(base, target), ]
rank_TAC[rownames(rank_TAC) %in% c(base, target), ]
zs_adult[rownames(zs_adult) %in% c(base, target), ]
zs_TAC[rownames(zs_TAC) %in% c(base, target), ]



adult_age_W3[rownames(adult_age_W3) %in% base, ]
sham_TAC_W3[rownames(sham_TAC_W3) %in% target, ]