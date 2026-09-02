# Shared input handling for the spASE analyses. Sourced by spase_scase.R
# (per-gene escape) and spase_spatial.R (2D spatial pattern per gene).
#
# It exists so the two scripts cannot drift. Both need the same gene x pixel
# matrices built the same way, and this project already has one case of the
# same-named quantity meaning three different things in three files
# (NEXT_ANALYSIS task 10, `allelic_ratio`). The matrix orientation in
# particular must never differ between them: matrix1 is ALT.
#
# Not meant to be run directly.

.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE     <- Sys.getenv("BASE",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
SAMPLE   <- Sys.getenv("SAMPLE", "9w")
X_CHROM  <- Sys.getenv("X_CHROM", "chrX")

# Which SNP mask the counting pass used, in the naming convention the rest of
# the spatial scripts follow: "no_Xist" keeps every path as it was, anything
# else reads and writes alongside it so two masks can be compared rather than
# one silently replacing the other.
SNP_LABEL <- Sys.getenv("SNP_LABEL", "no_Xist")
SUF       <- if (SNP_LABEL == "no_Xist") "" else paste0("_", SNP_LABEL)

GENE_BIN_UM <- as.integer(Sys.getenv("GENE_BIN_UM", "16"))
ASE_DIR   <- Sys.getenv("ASE_DIR", file.path(BASE, "ase", SAMPLE))
GENE_BINS <- Sys.getenv("GENE_BINS",
  file.path(ASE_DIR, sprintf("gene_bins_%dum%s.tsv.gz", GENE_BIN_UM, SUF)))
OUT_DIR   <- Sys.getenv("OUT_DIR", ASE_DIR)

# Gene must be seen on at least this many pixels, and carry at least this many
# informative UMIs in total, to be fitted. A beta-binomial with a free
# overdispersion needs pixels more than it needs depth per pixel, which is why
# the pixel grid is 16um and not 64um.
MIN_PIXELS <- as.integer(Sys.getenv("MIN_PIXELS", "100"))
MIN_UMI    <- as.integer(Sys.getenv("MIN_UMI", "50"))
CORES      <- as.integer(Sys.getenv("CORES", "1"))
# spASE's progress bars write several thousand carriage returns into a slurm
# log, where nothing consumes them. Off by default; VERBOSE=1 interactively.
VERBOSE    <- Sys.getenv("VERBOSE", "0") == "1"

GENE_TAG_OK <- Sys.getenv("GENE_TAG_OK", "0") == "1"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))
# ggplot does not wrap a subtitle, it clips it at the device edge, and every
# subtitle in these scripts carries numbers that have to stay readable.
wrap <- function(x, w = 96) paste(strwrap(x, w), collapse = "\n")

# Validated categorical slots 1 and 2. chrX is ALWAYS slot 1 and the autosomal
# control ALWAYS slot 2, in every panel of both scripts, so the same entity
# keeps its colour when the gene set changes.
C_X <- "#2a78d6"; C_A <- "#eb6834"; INK <- "#0b0b0b"; INK2 <- "#52514e"

spase_theme <- function() {
  theme_set(theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
          plot.background  = element_rect(fill = "#fcfcfb", colour = NA),
          plot.title    = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = INK2, size = 8),
          axis.text = element_text(colour = INK2)))
}

##### ------------------- DEPENDENCIES ------------------- #####
# spASE declares no Imports in its DESCRIPTION, so nothing is pulled in
# automatically and a missing dependency shows up as "lazy loading failed"
# during install rather than as a missing-package error here. The full set,
# read off the package source, is in slurm/spatial_spase.slurm.
load_spase <- function() {
  if (!requireNamespace("spASE", quietly = TRUE))
    stop("spASE is not installed. See the install block at the top of\n",
         "  slurm/spatial_spase.slurm - it lists the undeclared dependencies,\n",
         "  which must be installed BEFORE install_github or the build fails.",
         call. = FALSE)
  suppressPackageStartupMessages({
    library(spacexr)  # order matters: spASE overrides several spacexr functions
    library(spASE)
  })
}

##### ---------------------- INPUT ----------------------- #####
# Returns everything both scripts need, so neither builds a matrix of its own:
#   m_alt, m_ref  gene x pixel, matrix1 = ALT = CAST = the inactive X
#   coords        pixel, x, y - the covariates argument spase() expects
#   per_gene      one row per gene in the table, before filtering
#   keep          the genes that passed
#   gene_assignment / prov_path  for the provenance sidecar
load_spase_input <- function() {
  if (!file.exists(GENE_BINS))
    stop("No gene x pixel table at ", GENE_BINS, "\n",
         "  Produce it with ase_bin_allele_counts.py --gene-bin-out ",
         "--gene-bed ...\n  (slurm/spatial_spase.slurm runs the counting pass ",
         "for you).", call. = FALSE)

  # How the counting pass assigned molecules to genes. A GX-tag table is exonic
  # only and is not the thing this project measures elsewhere, so it is refused
  # by default rather than quietly analysed.
  prov_path <- Sys.getenv("PROVENANCE", "")
  if (!nzchar(prov_path)) {
    cand <- c(sub("\\.tsv\\.gz$", ".provenance.tsv", GENE_BINS),
              file.path(ASE_DIR, "bin_allele_counts.tsv.provenance.tsv"))
    prov_path <- cand[file.exists(cand)][1]
  }
  gene_assignment <- NA_character_
  if (!is.na(prov_path) && nzchar(prov_path) && file.exists(prov_path)) {
    prov <- fread(prov_path, sep = "\t", header = TRUE,
                  colClasses = "character")
    gv <- prov$value[prov$key == "gene_assignment"]
    if (length(gv)) gene_assignment <- gv[1]
  }
  if (!is.na(gene_assignment) && gene_assignment != "gene_bed" && !GENE_TAG_OK)
    stop("The gene table at\n  ", GENE_BINS, "\nwas built by gene assignment '",
         gene_assignment, "', not a gene-body bed. GX marks exonic assignment\n",
         "only and this data is mostly intronic, so most molecules are missing.\n",
         "Recount with --gene-bed, or set GENE_TAG_OK=1 to analyse it anyway.",
         call. = FALSE)

  cat(sprintf("Sample %s | pixels %dum | SNP mask %s | gene assignment %s\n",
              SAMPLE, GENE_BIN_UM, SNP_LABEL,
              if (is.na(gene_assignment)) "unrecorded" else gene_assignment))

  gb <- fread(GENE_BINS)
  setnames(gb, tolower(names(gb)))
  need <- c("chrom", "gene", "row", "col", "ref", "alt")
  if (!all(need %in% names(gb)))
    stop("Expected columns ", paste(need, collapse = ", "), " in ", GENE_BINS,
         "; got ", paste(names(gb), collapse = ", "), call. = FALSE)

  # A gene name appearing on two chromosomes would be silently pooled by the
  # matrix build, and the autosomal calibration set is defined by chromosome,
  # so this has to be an error rather than a warning.
  dup <- gb[, .(nchrom = uniqueN(chrom)), by = gene][nchrom > 1]
  if (nrow(dup))
    stop(nrow(dup), " gene name(s) occur on more than one chromosome in the ",
         "annotation, e.g. ", paste(head(dup$gene, 5), collapse = ", "),
         ". De-duplicate the gene bed.", call. = FALSE)

  gb[, pixel := paste0("r", row, "_c", col)]
  gb[, n := ref + alt]

  per_gene <- gb[, .(chrom = chrom[1], n_pixels = .N, umi = sum(n),
                     ref = sum(ref), alt = sum(alt)), by = gene]
  keep <- per_gene[n_pixels >= MIN_PIXELS & umi >= MIN_UMI]

  cat(sprintf(paste0("%d genes in the table (%d chrX); %d pass ",
                     ">=%d pixels and >=%d UMIs (%d chrX)\n"),
              nrow(per_gene), sum(per_gene$chrom == X_CHROM),
              nrow(keep), MIN_PIXELS, MIN_UMI, sum(keep$chrom == X_CHROM)))
  if (sum(keep$chrom == X_CHROM) == 0)
    stop("No chrX gene passes the filters - nothing to test.", call. = FALSE)
  n_auto <- sum(keep$chrom != X_CHROM)
  if (n_auto < 20)
    warning("Only ", n_auto, " autosomal genes pass the filters. The autosomal ",
            "control below is derived from them and will be unstable.",
            call. = FALSE)

  gbk    <- gb[gene %in% keep$gene]
  genes  <- sort(unique(gbk$gene))
  pixels <- sort(unique(gbk$pixel))
  gi <- match(gbk$gene, genes); pj <- match(gbk$pixel, pixels)
  dims <- c(length(genes), length(pixels))
  dn   <- list(genes, pixels)

  # matrix1 = ALT = CAST = the inactive X, so scase()'s p is the escape
  # fraction directly and spase()'s surface is an escape surface.
  m_alt <- sparseMatrix(i = gi, j = pj, x = as.numeric(gbk$alt),
                        dims = dims, dimnames = dn)
  m_ref <- sparseMatrix(i = gi, j = pj, x = as.numeric(gbk$ref),
                        dims = dims, dimnames = dn)
  cat(sprintf("Matrices: %d genes x %d pixels, %.1f%% non-zero\n",
              dims[1], dims[2], 100 * nrow(gbk) / prod(dims)))

  # spase() wants pixel ID first, then the two coordinates. Coordinates are the
  # PIXEL indices (row//k, col//k), so one unit is GENE_BIN_UM microns; spase()
  # standardises them internally anyway, but anything that reads a distance off
  # a fitted surface has to know the unit.
  px <- unique(gbk[, .(pixel, row, col)])
  coords <- data.frame(pixel = px$pixel, x = px$col, y = px$row,
                       stringsAsFactors = FALSE)
  coords <- coords[match(pixels, coords$pixel), ]
  rownames(coords) <- NULL

  list(m_alt = m_alt, m_ref = m_ref, coords = coords, per_gene = per_gene,
       keep = keep, genes = genes, pixels = pixels,
       gene_assignment = gene_assignment, prov_path = prov_path)
}

# Provenance keys both scripts write. Task 9's rule: every input that shaped a
# result is named next to it.
spase_provenance_common <- function(dat) {
  md5 <- tryCatch(tools::md5sum(GENE_BINS)[[1]],
                  error = function(e) NA_character_)
  c(paste0("sample\t", SAMPLE),
    paste0("gene_bins\t", GENE_BINS),
    paste0("gene_bins_md5\t", md5),
    paste0("gene_assignment\t", dat$gene_assignment),
    paste0("counter_provenance\t",
           if (is.na(dat$prov_path)) "" else dat$prov_path),
    paste0("snp_label\t", SNP_LABEL),
    paste0("gene_bin_um\t", GENE_BIN_UM),
    paste0("min_pixels\t", MIN_PIXELS),
    paste0("min_umi\t", MIN_UMI),
    paste0("n_pixels\t", length(dat$pixels)),
    paste0("spASE_version\t", as.character(packageVersion("spASE"))),
    paste0("R_version\t", paste(R.version$major, R.version$minor, sep = ".")),
    paste0("run_utc\t", format(Sys.time(), tz = "UTC")))
}
