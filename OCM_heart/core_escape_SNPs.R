library(tidyverse)
library(data.table)

annotation <- fread("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_us.bed")
SNPfile <- fread("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/SNPs/SNPfile_C57BL6_NJ_CAST_EiJ.bed")

# create annotation summed up for all chromosome
annotation_summed <- annotation %>% group_by(V1) %>% summarise(V2 = min(V2), V3 = max(V3), .groups = "drop")


# Subset SNP file for all genes excluding Xist
Xist_annotation <- annotation[V4 == "Xist",]
Xist_span <- seq(Xist_annotation$V2, Xist_annotation$V3, by = 1)
SNPfile_no_Xist <- SNPfile %>%  dplyr::filter(!(V1 == "chrX" & V2 %in% Xist_span))
# Save file
write.table(SNPfile_no_Xist, "SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed",row.names=F,col.names=F,quote=F,sep="\t")

# Generate annotation file for just the X chromosome genes
annotation_chrX <- annotation[V1 == "chrX",]
write.table(annotation_chrX, "annotation_us_mm39_chrX.bed", row.names=F, col.names=F, quote=F, sep="\t")

# Subset SNP file for core escape genes
core_escape_genes <- c("Kdm5c", "Kdm6a", "Ddx3x", "Eif2s3x")

#load SNP file
core_escape <- annotation[V4 %in% core_escape_genes,]

# create a vector with all positions
core_escape_span <- unlist(Map(seq, core_escape$V2, core_escape$V3))

#filter SNPfile
SNPfile_core_escape <- SNPfile %>%  dplyr::filter(V1 == "chrX" & V2 %in% core_escape_span)

#Save file
write.table(SNPfile_core_escape, "SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_core_escape.bed",row.names=F,col.names=F,quote=F,sep="\t")



# Collapse transcript rows to one gene interval per gene
core_escape_gene_level <- core_escape %>%
	dplyr::group_by(V1, V4, V6) %>%
	dplyr::summarise(
		V2 = min(V2),
		V3 = max(V3),
		V5 = 0L,
		.groups = "drop"
	) %>%
	dplyr::select(V1, V2, V3, V4, V5, V6)

# Save gene-level annotation file
write.table(core_escape_gene_level, "annotation_us_mm39_core_escape.bed", row.names=F, col.names=F, quote=F, sep="\t")
