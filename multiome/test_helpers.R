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

# ---- sample order ----
# The default is WRONG rather than arbitrary here: sorted as text "78w" comes
# before "9w", so anything left to ggplot presents aged before adult.
chk(identical(SAMPLE_LEVELS, c("9w", "78w")),        "adult before aged")
chk(identical(sort(c("9w", "78w")), c("78w", "9w")), "and plain sorting reverses it (the bug)")

sx <- as_sample(c("78w", "9w", "78w"))
chk(is.factor(sx),                                "as_sample returns a factor")
chk(identical(levels(sx), c("9w", "78w")),        "levels in display order, not data order")
chk(identical(levels(as_sample(factor(c("78w", "9w")))), c("9w", "78w")),
    "a factor arriving in the wrong order is re-levelled")
chk(is.na(as_sample("Sham")),                     "an unknown sample becomes NA, not a silent level")

# limits, not just values: that is what pins the order when the column reached
# ggplot as plain character.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  for (a in c("colour", "fill")) {
    sc <- sample_scale(a)
    chk(identical(sc$limits, SAMPLE_LEVELS), paste0("sample_scale('", a, "') pins limits"))
  }
  chk(identical(sample_x()$limits, SAMPLE_LEVELS), "sample_x pins limits")
} else {
  cat("SKIP   sample_scale tests (ggplot2 not installed)\n")
}

cat("\n", if (fail == 0) "ALL PASS" else paste(fail, "FAILURES"), "\n")
quit(status = if (fail == 0) 0 else 1)
