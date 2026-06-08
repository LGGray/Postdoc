# Export Seurat objects to h5ad for scVI / Google Colab
# Run this BEFORE uploading data to Google Drive
#
# Required package for robust h5ad export:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("zellkonverter")

library(Seurat)
library(zellkonverter)

get_joined_layer_matrix <- function(obj, assay = 'RNA', base_layer = 'counts') {
	layers <- Layers(obj[[assay]])

	# Prefer canonical layer if present.
	if (base_layer %in% layers) {
		return(LayerData(obj, assay = assay, layer = base_layer))
	}

	# Otherwise merge split layers (e.g. counts.1, counts.2, ...).
	pat <- paste0('^', base_layer, '(\\.|$)')
	matched <- grep(pat, layers, value = TRUE)
	if (length(matched) == 0) {
		stop('No layer found for base layer: ', base_layer)
	}

	mats <- lapply(matched, function(ln) LayerData(obj, assay = assay, layer = ln))
	merged <- do.call(cbind, mats)

	# Reorder columns to match Seurat object cell order.
	common_cells <- intersect(colnames(obj), colnames(merged))
	if (length(common_cells) != ncol(obj)) {
		stop('Merged layer matrix does not contain all Seurat cells for layer ', base_layer)
	}
	merged <- merged[, colnames(obj), drop = FALSE]

	merged
}

ensure_rna_data <- function(obj, object_name) {
	DefaultAssay(obj) <- 'RNA'

	# Keep only RNA for clean export.
	obj <- DietSeurat(obj, assays = 'RNA', dimreducs = NULL, graphs = NULL)
	obj <- JoinLayers(obj)

	# Keep a data layer for compatibility with downstream tools.
	has_data_layer <- tryCatch({
		'data' %in% Layers(obj[['RNA']])
	}, error = function(e) {
		FALSE
	})

	if (!has_data_layer) {
		cat(object_name, ': creating RNA data layer with NormalizeData()\n', sep = '')
		obj <- NormalizeData(obj, assay = 'RNA', verbose = FALSE)
	}

	obj
}

read_first_existing_rds <- function(paths) {
	for (p in paths) {
		if (file.exists(p)) {
			cat('Loading: ', p, '\n', sep = '')
			return(readRDS(p))
		}
	}
	stop('Could not find any of these files:\n', paste('  -', paths, collapse = '\n'))
}

write_h5ad_counts_x <- function(obj, outfile, object_name) {
	# Build SingleCellExperiment manually to avoid GetAssayData() failures on Seurat v5 multilayer assays.
	counts_mat <- get_joined_layer_matrix(obj, assay = 'RNA', base_layer = 'counts')
	data_layers <- Layers(obj[['RNA']])
	has_data <- any(grepl('^data(\\.|$)', data_layers))

	sce <- SingleCellExperiment::SingleCellExperiment(
		assays = list(counts = counts_mat),
		colData = S4Vectors::DataFrame(obj@meta.data)
	)

	if (has_data) {
		logcounts_mat <- get_joined_layer_matrix(obj, assay = 'RNA', base_layer = 'data')
		SummarizedExperiment::assay(sce, 'logcounts') <- logcounts_mat
	}

	zellkonverter::writeH5AD(sce, file = outfile, X_name = 'counts')
	cat(object_name, ': wrote ', outfile, '\n', sep = '')
}

# -------------------------------------------------------------------
# Reference: adult_aged_merged.RDS (ASE/)
#   - 9w adult + 78w aged mouse cardiac snRNA-seq
#   - Annotated cell types in the 'celltype' column:
#     FB, ventCM, EC, SMC, Macro, atrialCM, LEC, Epicard, Adipo, Glial
#   - Batch variable: 'age' (values: "adult", "aged")
# -------------------------------------------------------------------
reference <- read_first_existing_rds(c(
	'adult_aged_merged.RDS',
	'adult_aged_merged.RDS'
))
DefaultAssay(reference) <- 'RNA'
reference <- ensure_rna_data(reference, 'Reference')

cat('Reference cells:', ncol(reference), '\n')
cat('Cell types:\n'); print(table(reference$celltype))
cat('Ages:\n');       print(table(reference$age))

# -------------------------------------------------------------------
# Query: heart_seurat_object_SCT.rds (OCM_heart/)
#   - 9w, 78w, TAC (28d), Sham (28d) mouse cardiac snRNA-seq
#   - No cell-type labels — these will be predicted by scANVI
#   - Batch variable: 'sample' (values: "9w", "78w", "TAC", "Sham")
# -------------------------------------------------------------------
query <- read_first_existing_rds(c(
	'heart_seurat_object_SCT.rds',
	'../OCM_heart/heart_seurat_object_SCT.rds',
	'../OCM/heart_seurat_object_SCT.rds'
))
DefaultAssay(query) <- 'RNA'
query <- ensure_rna_data(query, 'Query')

cat('Query cells:', ncol(query), '\n')
cat('Samples:\n'); print(table(query$sample))

# -------------------------------------------------------------------
# Export directly to h5ad using zellkonverter.
# AnnData X is set to raw counts (required by scVI).
# -------------------------------------------------------------------

# Reference
write_h5ad_counts_x(reference, outfile = 'reference_adult_aged.h5ad', object_name = 'Reference')
# -> produces reference_adult_aged.h5ad

# Query
write_h5ad_counts_x(query, outfile = 'query_OCM.h5ad', object_name = 'Query')
# -> produces query_OCM.h5ad

cat('\nDone. Upload the following files to Google Drive:\n')
cat('  reference_adult_aged.h5ad\n')
cat('  query_OCM.h5ad\n')
