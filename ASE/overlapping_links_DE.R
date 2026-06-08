library(dplyr)
library(ggplot2)
library(ggrepel)

link_key <- function(df) paste(df$name_base, df$name_target, df$mechanism, sep="|")
sample_colours <- c("adult" = "#E69F00", "aged" = "#56B4E9", "TAC" = "#009E73", "Sham" = "#CC79A7")

extract_link_genes <- function(link_ids) {
	if (length(link_ids) == 0) {
		return(list(lncRNA = character(0), pcGene = character(0)))
	}

	split_links <- strsplit(link_ids, "\\|")
	list(
		lncRNA = unique(vapply(split_links, `[`, character(1), 1)),
		pcGene = unique(vapply(split_links, `[`, character(1), 2))
	)
}

compute_link_sets <- function(group_a_df, group_b_df) {
	group_a_links <- link_key(group_a_df)
	group_b_links <- link_key(group_b_df)

	group_a_only <- setdiff(group_a_links, group_b_links)
	group_b_only <- setdiff(group_b_links, group_a_links)
	common_links <- intersect(group_a_links, group_b_links)

	list(
		group_a_only = group_a_only,
		group_b_only = group_b_only,
		common = common_links,
		group_a_genes = extract_link_genes(group_a_only),
		group_b_genes = extract_link_genes(group_b_only),
		common_genes = extract_link_genes(common_links)
	)
}

compute_de_sets <- function(de_file, lfc_cutoff = 1, padj_cutoff = 0.05) {
	de_table <- read.delim(de_file, header = TRUE)

	list(
		de_table = de_table,
		deg = subset(de_table, abs(log2FoldChange) > lfc_cutoff & padj < padj_cutoff)$gene,
		up_group_a = subset(de_table, log2FoldChange < -lfc_cutoff & padj < padj_cutoff)$gene,
		up_group_b = subset(de_table, log2FoldChange > lfc_cutoff & padj < padj_cutoff)$gene
	)
}

hits_to_df <- function(cell_type, hit_group, genes) {
	if (length(genes) == 0) {
		return(data.frame(cell_type = character(0), hit_group = character(0), gene = character(0)))
	}

	data.frame(
		cell_type = rep(cell_type, length(genes)),
		hit_group = rep(hit_group, length(genes)),
		gene = genes,
		stringsAsFactors = FALSE
	)
}

analyse_celltype_overlap <- function(group_a_df, group_b_df, de_file, cell_type,
																		 group_a_label = "adult", group_b_label = "aged") {
	link_sets <- compute_link_sets(group_a_df, group_b_df)
	de_sets <- compute_de_sets(de_file)

	hits <- list(
		group_a_only_lncRNA_up_group_a = link_sets$group_a_genes$lncRNA[link_sets$group_a_genes$lncRNA %in% de_sets$up_group_a],
		group_a_only_pcGene_up_group_a = link_sets$group_a_genes$pcGene[link_sets$group_a_genes$pcGene %in% de_sets$up_group_a],
		group_b_only_lncRNA_up_group_b = link_sets$group_b_genes$lncRNA[link_sets$group_b_genes$lncRNA %in% de_sets$up_group_b],
		group_b_only_pcGene_up_group_b = link_sets$group_b_genes$pcGene[link_sets$group_b_genes$pcGene %in% de_sets$up_group_b],
		common_lncRNA_in_deg = link_sets$common_genes$lncRNA[link_sets$common_genes$lncRNA %in% de_sets$deg],
		common_pcGene_in_deg = link_sets$common_genes$pcGene[link_sets$common_genes$pcGene %in% de_sets$deg]
	)

	summary <- data.frame(
		cell_type = rep(cell_type, 6),
		comparison = rep(paste(group_a_label, "vs", group_b_label), 6),
		metric = c(
			"group_a_only_lncRNA_up_group_a",
			"group_a_only_pcGene_up_group_a",
			"group_b_only_lncRNA_up_group_b",
			"group_b_only_pcGene_up_group_b",
			"common_lncRNA_in_deg",
			"common_pcGene_in_deg"
		),
		n = c(
			length(hits$group_a_only_lncRNA_up_group_a),
			length(hits$group_a_only_pcGene_up_group_a),
			length(hits$group_b_only_lncRNA_up_group_b),
			length(hits$group_b_only_pcGene_up_group_b),
			length(hits$common_lncRNA_in_deg),
			length(hits$common_pcGene_in_deg)
		),
		stringsAsFactors = FALSE
	)

	detailed_hits <- bind_rows(
		hits_to_df(cell_type, "group_a_only_lncRNA_up_group_a", hits$group_a_only_lncRNA_up_group_a),
		hits_to_df(cell_type, "group_a_only_pcGene_up_group_a", hits$group_a_only_pcGene_up_group_a),
		hits_to_df(cell_type, "group_b_only_lncRNA_up_group_b", hits$group_b_only_lncRNA_up_group_b),
		hits_to_df(cell_type, "group_b_only_pcGene_up_group_b", hits$group_b_only_pcGene_up_group_b),
		hits_to_df(cell_type, "common_lncRNA_in_deg", hits$common_lncRNA_in_deg),
		hits_to_df(cell_type, "common_pcGene_in_deg", hits$common_pcGene_in_deg)
	)

	list(
		link_sets = link_sets,
		de_sets = de_sets,
		hits = hits,
		summary = summary,
		detailed_hits = detailed_hits
	)
}

################################
# Load in adult and aged links #
################################
load('adult_aged_bodymap/bodymap_links.RData')
adult_heart <- link_key(bodymap[['He_9w']])
aged_heart <- link_key(bodymap[['He_78w']])

adult_links <- setdiff(adult_heart, aged_heart)
adult_lncRNA <- strsplit(adult_links, "\\|") %>% sapply("[", 1) %>% unique()
adult_pcGene <- strsplit(adult_links, "\\|") %>% sapply("[", 2) %>% unique()
aged_links <- setdiff(aged_heart, adult_heart)
aged_lncRNA <- strsplit(aged_links, "\\|") %>% sapply("[", 1) %>% unique()
aged_pcGene <- strsplit(aged_links, "\\|") %>% sapply("[", 2) %>% unique()
adult_aged_common_links <- intersect(adult_heart, aged_heart)
adult_aged_common_lncRNA <- strsplit(adult_aged_common_links, "\\|") %>% sapply("[", 1) %>% unique()
adult_aged_common_pcGene <- strsplit(adult_aged_common_links, "\\|") %>% sapply("[", 2) %>% unique()

# Read in adult and aged DE results
adult_aged_deseq2 <- read.delim('adult_aged_bodymap/DESeq2/He_9w_vs_78w.txt', header = TRUE)
adult_aged_deg <- subset(adult_aged_deseq2, abs(log2FoldChange) > 1 & padj < 0.05)$gene
up_adult <- subset(adult_aged_deseq2, log2FoldChange < -1 & padj < 0.05)$gene
up_aged <- subset(adult_aged_deseq2, log2FoldChange > 1 & padj < 0.05)$gene


adult_lncRNA[adult_lncRNA %in% up_adult]
adult_pcGene[adult_pcGene %in% up_adult]
aged_lncRNA[aged_lncRNA %in% up_aged]
aged_pcGene[aged_pcGene %in% up_aged]
adult_aged_common_lncRNA[adult_aged_common_lncRNA %in% adult_aged_deg]
adult_aged_common_pcGene[adult_aged_common_pcGene %in% adult_aged_deg]

#########################################
# Load in TAC and Sham bulk heart links #
#########################################
load('F1_TAC_Sarah/TAC_SHAM_links.RData')
sham_heart <- link_key(TAC_SHAM[["He_TACSham28d_RNA_XistBxC_-+_XX"]])
tac_heart <- link_key(TAC_SHAM[["He_TAC28d_RNA_XistBxC_-+_XX"]])

tac_links <- setdiff(tac_heart, sham_heart)
tac_lncRNA <- strsplit(tac_links, "\\|") %>% sapply("[", 1) %>% unique()
tac_pcGene <- strsplit(tac_links, "\\|") %>% sapply("[", 2) %>% unique()
sham_links <- setdiff(sham_heart, tac_heart)
sham_lncRNA <- strsplit(sham_links, "\\|") %>% sapply("[", 1) %>% unique()
sham_pcGene <- strsplit(sham_links, "\\|") %>% sapply("[", 2) %>% unique()
tac_sham_common_links <- intersect(tac_heart, sham_heart)
tac_sham_common_lncRNA <- strsplit(tac_sham_common_links, "\\|") %>% sapply("[", 1) %>% unique()
tac_sham_common_pcGene <- strsplit(tac_sham_common_links, "\\|") %>% sapply("[", 2) %>% unique()

# Read in TAC and Sham DE results
tac_sham_deseq2 <- read.delim('F1_TAC_Sarah/DEG/TAC_vs_sham.txt', header = TRUE)
tac_sham_deg <- subset(tac_sham_deseq2, abs(log2FoldChange) > 1 & padj < 0.05)$gene
up_tac <- subset(tac_sham_deseq2, log2FoldChange > 1 & padj < 0.05)$gene
up_sham <- subset(tac_sham_deseq2, log2FoldChange < -1 & padj < 0.05)$gene


tac_lncRNA[tac_lncRNA %in% up_tac]
tac_pcGene[tac_pcGene %in% up_tac]
sham_lncRNA[sham_lncRNA %in% up_sham]
sham_pcGene[sham_pcGene %in% up_sham]
tac_sham_common_lncRNA[tac_sham_common_lncRNA %in% tac_sham_deg]
tac_sham_common_pcGene[tac_sham_common_pcGene %in% tac_sham_deg]

# Load in adult FACS links
load('cardiac_RNAseq/adult_facs_links.RData')
load('aged_cardiac_RNAseq/aged_facs_links.RData')

# Configure once, then add more rows here for additional cell types.
facs_celltype_config <- data.frame(
	cell_type = c("CF"),
	adult_key = c("CF_adult"),
	aged_key = c("CF_aged"),
	de_file = c("LINKS_study/DEG/CF_aged_vs_adult.txt"),
	stringsAsFactors = FALSE
)

facs_link_results <- lapply(seq_len(nrow(facs_celltype_config)), function(i) {
	cfg <- facs_celltype_config[i, ]
	analyse_celltype_overlap(
		group_a_df = adult_facs_links[[cfg$adult_key]],
		group_b_df = aged_facs_links[[cfg$aged_key]],
		de_file = cfg$de_file,
		cell_type = cfg$cell_type,
		group_a_label = "adult",
		group_b_label = "aged"
	)
})
names(facs_link_results) <- facs_celltype_config$cell_type

facs_output_dir <- "LINKS_study/DEG/facs_link_overlap"
if (!dir.exists(facs_output_dir)) {
	dir.create(facs_output_dir, recursive = TRUE)
}

facs_summary <- bind_rows(lapply(facs_link_results, `[[`, "summary"))
facs_detailed_hits <- bind_rows(lapply(facs_link_results, `[[`, "detailed_hits"))

write.table(
	facs_summary,
	file = file.path(facs_output_dir, "all_celltypes_summary.tsv"),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE
)
write.table(
	facs_detailed_hits,
	file = file.path(facs_output_dir, "all_celltypes_hits.tsv"),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE
)

for (cell_type in names(facs_link_results)) {
	result <- facs_link_results[[cell_type]]
	write.table(
		result$summary,
		file = file.path(facs_output_dir, paste0(cell_type, "_summary.tsv")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)
	write.table(
		result$detailed_hits,
		file = file.path(facs_output_dir, paste0(cell_type, "_hits.tsv")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)
}

# Keep CF objects available for downstream code compatibility.
CF_results <- facs_link_results[["CF"]]

adult_CF_links <- CF_results$link_sets$group_a_only
adult_CF_lncRNA <- CF_results$link_sets$group_a_genes$lncRNA
adult_CF_pcGene <- CF_results$link_sets$group_a_genes$pcGene
aged_CF_links <- CF_results$link_sets$group_b_only
aged_CF_lncRNA <- CF_results$link_sets$group_b_genes$lncRNA
aged_CF_pcGene <- CF_results$link_sets$group_b_genes$pcGene
adult_aged_CF_common_links <- CF_results$link_sets$common
adult_aged_CF_common_lncRNA <- CF_results$link_sets$common_genes$lncRNA
adult_aged_CF_common_pcGene <- CF_results$link_sets$common_genes$pcGene

adult_aged_CF_deseq2 <- CF_results$de_sets$de_table
adult_aged_CF_deg <- CF_results$de_sets$deg
up_adult_CF <- CF_results$de_sets$up_group_a
up_aged_CF <- CF_results$de_sets$up_group_b

adult_CF_lncRNA[adult_CF_lncRNA %in% up_adult_CF]
adult_CF_pcGene[adult_CF_pcGene %in% up_adult_CF]
aged_CF_lncRNA[aged_CF_lncRNA %in% up_aged_CF]
aged_CF_pcGene[aged_CF_pcGene %in% up_aged_CF]
adult_aged_CF_common_lncRNA[adult_aged_CF_common_lncRNA %in% adult_aged_CF_deg]
adult_aged_CF_common_pcGene[adult_aged_CF_common_pcGene %in% adult_aged_CF_deg]

# Load in aged FACS links
load('aged_cardiac_RNAseq/aged_facs_links.RData')




# Load in TAC and Sham FACS links
load('TAC_cardiac_RNAseq/sham_facs_links.RData')