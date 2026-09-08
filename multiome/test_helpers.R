# ---------------------------------------------------------------------------
# Tests for multiome/00_helpers.R. No Seurat required - it stubs an S4 object
# with a meta.data slot, so this runs anywhere R does:
#
#   Rscript multiome/test_helpers.R
#
# Worth having because ONE bug shape cost three job submissions from three
# call sites: a named-vector lookup carries the LOOKUP KEYS as names, Seurat's
# `obj$col <-` matches metadata on names, and the result is always
# "No cell overlap between new meta data and Seurat object". Run this after
# touching the helpers.
# ---------------------------------------------------------------------------
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_helpers.R"))
fail <- 0
chk <- function(ok, what) { cat(if (ok) "PASS  " else "FAIL  ", what, "\n"); if (!ok) fail <<- fail + 1 }

# ---- short_labels ----
full <- c("Ventricular Cardiomyocytes", "Fibroblasts",
          "Pericytes - Smooth muscle cells", "Some Unmapped Type")
s <- short_labels(full)
cat("short_labels ->", paste(s, collapse=" | "), "\n")
chk(is.null(names(s)),                      "short_labels returns UNNAMED (the bug)")
chk(length(s) == length(full),              "length preserved")
chk(s[1] == "Ventricular CM",               "known label shortened")
chk(s[3] == "Pericyte/SMC",                 "long label shortened")
chk(s[4] == "Some Unmapped Type",           "unmapped label passes through")

# ---- short_labels on a factor, which is what Seurat metadata often holds ----
sf <- short_labels(factor(full))
chk(is.null(names(sf)) && identical(unname(sf), unname(s)), "factor input handled")

# ---- set_meta, against a stub S4 object with a meta.data slot ----
setClass("FakeSeurat", representation(meta.data = "data.frame"))
o <- new("FakeSeurat",
         meta.data = data.frame(cluster = c("0","0","1"),
                                row.names = paste0("9w_BC", 1:3)))
named <- c("Ventricular Cardiomyocytes" = "Ventricular CM",
           "Ventricular Cardiomyocytes" = "Ventricular CM",
           "Fibroblasts"                = "Fibroblasts")
o2 <- set_meta(o, "celltype_short", named)
chk(is.null(names(o2@meta.data$celltype_short)), "set_meta strips names")
chk(identical(unname(o2@meta.data$celltype_short), unname(named)), "values preserved in order")
chk(rownames(o2@meta.data)[1] == "9w_BC1",       "rownames untouched")

len_err <- tryCatch({ set_meta(o, "x", c("a","b")); "no error" },
                    error = function(e) conditionMessage(e))
chk(grepl("2 values for 3 cells", len_err), "set_meta rejects a length mismatch")
cat("  ->", len_err, "\n")

# ---- celltype_scale fixed order / no cycling ----
chk(length(OKABE_ITO) == 8, "8 fixed colours")
chk(is.null(celltype_scale(paste0("t", 1:9))), "9 categories -> NULL, no invented hue")

cat("\n", if (fail == 0) "ALL PASS" else paste(fail, "FAILURES"), "\n")
quit(status = if (fail == 0) 0 else 1)
