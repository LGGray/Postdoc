# Patch for Banksy:::computeHarmonics  (Banksy issue prabhakarlab/Banksy#77)
#
# PROBLEM
#   gcm_chunk is a sparse S4 Matrix, and the harmonic multiplier
#   exp(j * M * phi) is complex, so
#       gcm_chunk[, to] %*% (weight * exp(j * M * phi))
#   returns a complex S4 Matrix. abs() takes the modulus but the result is
#   still S4, and data.table rejects S4 values in j:
#       Error in `[.data.table`(knn_df, , abs(gcm_chunk[, to, drop = FALSE] %*% :
#         j evaluates to type 'S4'. Must evaluate to atomic vector or list.
#
#   Note this requires Matrix >= 1.8-0 to even get this far; on Matrix 1.7-5
#   the %*% fails earlier with 'zgeMatrix is not a defined class' (issue #76).
#
# FIX
#   Wrap each abs(...) in as.numeric(). This is the ONLY change - the function
#   below is otherwise a verbatim copy of Banksy:::computeHarmonics as
#   installed. Search for "PATCHED" to find the two edited lines.
#
#   Why it is safe: per group (by = from, i.e. per cell) the product is an
#   n_genes x 1 column vector. as.numeric() flattens column-major, which is a
#   no-op on a single column, so the value order is unchanged. The caller
#   immediately does matrix(chunk_results[[i]], nrow = ..., ncol = ncol(gcm)),
#   refilling column-major, so column i is still cell i. dim is discarded
#   there regardless. abs() already returned real moduli, so no numerical
#   change occurs - only the S4 wrapper is removed.
#
# USE
#   library(Banksy); source("spatial/banksy_patch.R")   # before RunBanksy()

.computeHarmonics_patched <- function(gcm, knn_df, M, center, verbose,
                                      chunk_size = NULL, parallel = FALSE,
                                      num_cores = NULL, row_limit_factor = 0.75)
{
    from <- to <- weight <- phi <- .N <- count <- . <- NULL
    j <- sqrt(as.complex(-1))
    mean_k <- round(mean(knn_df[, .(count = .N), by = from]$count), 1)
    total_rows <- as.double(nrow(gcm)) * ncol(gcm)
    max_rows <- (2^31 - 1) * row_limit_factor
    if (total_rows > max_rows || !is.null(chunk_size)) {
        if (verbose)
            message("Computing neighborhood matrices in chunks...")
        max_chunk_size <- floor(max_rows/ncol(gcm))
        if (!is.null(chunk_size)) {
            if (chunk_size > max_chunk_size) {
                stop("Specified chunk_size too large. Must be smaller than ",
                     floor(max_rows/ncol(gcm)))
            }
            max_chunk_size <- chunk_size
        }
    }
    else {
        max_chunk_size <- nrow(gcm)
    }
    num_chunks <- ceiling(nrow(gcm)/max_chunk_size)
    if (verbose && num_chunks > 1) {
        message("Processing in ", num_chunks, " chunks of max ",
                max_chunk_size, " genes each")
        if (parallel) {
            if (.Platform$OS.type == "windows") {
                message("Parallel processing not supported on Windows - using sequential processing")
            }
            else {
                if (!is.null(num_cores)) {
                    message("Using parallel processing with ", num_cores, " cores")
                }
                else {
                    message("Using parallel processing with default backend")
                }
            }
        }
    }
    if (parallel) {
        if (.Platform$OS.type == "windows") {
            if (verbose)
                message("Parallel processing not supported on Windows. Using sequential processing.")
            apply_fun <- lapply
        }
        else {
            if (!requireNamespace("BiocParallel", quietly = TRUE)) {
                stop("BiocParallel package is required for parallel processing. Please install it with: BiocManager::install('BiocParallel')")
            }
            if (!is.null(num_cores)) {
                bp_param <- BiocParallel::MulticoreParam(workers = num_cores)
                apply_fun <- function(x, fun) BiocParallel::bplapply(x, fun,
                                                                    BPPARAM = bp_param)
            }
            else {
                apply_fun <- BiocParallel::bplapply
            }
        }
    }
    else {
        apply_fun <- lapply
    }
    process_chunk <- function(chunk) {
        start_idx <- (chunk - 1) * max_chunk_size + 1
        end_idx <- min(chunk * max_chunk_size, nrow(gcm))
        if (verbose && num_chunks > 1 && (!parallel || .Platform$OS.type ==
            "windows")) {
            message("Processing chunk ", chunk, "/", num_chunks,
                    " (rows ", start_idx, " to ", end_idx, ")")
        }
        gcm_chunk <- gcm[start_idx:end_idx, , drop = FALSE]
        if (center) {
            if (verbose && chunk == 1 && (!parallel || .Platform$OS.type ==
                "windows"))
                message("Centering")
            ## PATCHED: as.numeric() wrapped around abs()
            chunk_aggr <- knn_df[, as.numeric(abs(fscale(gcm_chunk[, to,
                drop = FALSE]) %*% (weight * exp(j * M * phi)))),
                by = from]
        }
        else {
            ## PATCHED: as.numeric() wrapped around abs()
            chunk_aggr <- knn_df[, as.numeric(abs(gcm_chunk[, to, drop = FALSE] %*%
                (weight * exp(j * M * phi)))), by = from]
        }
        result_vector <- chunk_aggr$V1
        rm(gcm_chunk, chunk_aggr)
        return(result_vector)
    }
    chunk_results <- apply_fun(1:num_chunks, process_chunk)
    ncm <- matrix(0, nrow = nrow(gcm), ncol = ncol(gcm))
    for (i in seq_along(chunk_results)) {
        start_idx <- (i - 1) * max_chunk_size + 1
        end_idx <- min(i * max_chunk_size, nrow(gcm))
        ncm[start_idx:end_idx, ] <- matrix(chunk_results[[i]],
            nrow = end_idx - start_idx + 1, ncol = ncol(gcm))
    }
    rm(chunk_results)
    rownames(ncm) <- rownames(gcm)
    colnames(ncm) <- colnames(gcm)
    if (verbose)
        message("Done")
    return(ncm)
}

# Bind into the Banksy namespace so internals (fscale) and data.table's
# imported [ method resolve correctly.
environment(.computeHarmonics_patched) <- asNamespace("Banksy")
assignInNamespace("computeHarmonics", .computeHarmonics_patched, ns = "Banksy")

message("Banksy:::computeHarmonics patched (as.numeric around abs; issue #77)")
