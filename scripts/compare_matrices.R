#!/usr/bin/env Rscript
#
# compare_matrices.R -- cell-by-cell comparison of two feature matrices.
#
# Use it to check a pipeline run against a reference matrix (for example the
# published one), or two runs against each other.
#
# Usage:
#   Rscript scripts/compare_matrices.R <reference.txt> <ours.txt> [--csv out.csv]
#
# Both files need a `sample` column; everything else is treated as a feature.
# Comparison is on the intersection of samples and of features, and the summary
# separates three different things that "differ" can mean:
#
#   * both sides have a value and the values disagree  -- a real disagreement
#   * one side is NA                                   -- an availability
#                                                         difference, usually an
#                                                         assay that was not done
#   * both NA                                          -- agreement
#
# hg19-era gene symbols in the reference are translated to current ones, so a
# matrix written with FAM46C / HIST1H1E lines up with one written with
# TENT5C / H1.4.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: Rscript scripts/compare_matrices.R <reference.txt> <ours.txt> [--csv out.csv]",
       call. = FALSE)
}
ref.file <- args[[1]]
our.file <- args[[2]]
csv.out <- if ("--csv" %in% args) args[[which(args == "--csv") + 1L]] else NA_character_

rd <- function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample" %in% colnames(d)) {
    #### Some outputs carry the sample as row names instead.
    if (!is.null(rownames(d)) && !identical(rownames(d), as.character(seq_len(nrow(d))))) {
      d$sample <- rownames(d)
    } else {
      stop(f, ": no 'sample' column and no usable row names", call. = FALSE)
    }
  }
  d
}

#### Gene symbols renamed between hg19 and hg38, so either spelling lines up.
ALIASES <- c(FAM46C = "TENT5C", HIST1H1E = "H1.4", HIST1H1C = "H1.2",
             HIST1H1D = "H1.3", HIST1H1B = "H1.5", HIST1H2BK = "H2BC12",
             TBC1D29 = "TBC1D29P")
canonical <- function(x) {
  for (old in names(ALIASES)) x <- sub(paste0("(^|_)", old, "$"),
                                      paste0("\\1", ALIASES[[old]]), x)
  x
}

A0 <- rd(ref.file); B0 <- rd(our.file)
colnames(A0) <- canonical(colnames(A0))
colnames(B0) <- canonical(colnames(B0))

samples  <- intersect(A0$sample, B0$sample)
features <- setdiff(intersect(colnames(A0), colnames(B0)), "sample")
if (!length(samples))  stop("the two matrices share no samples", call. = FALSE)
if (!length(features)) stop("the two matrices share no features", call. = FALSE)

num <- function(d) {
  m <- as.matrix(d[match(samples, d$sample), features, drop = FALSE])
  matrix(suppressWarnings(as.numeric(m)), nrow = length(samples),
         dimnames = list(samples, features))
}
A <- num(A0); B <- num(B0)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("COVERAGE\n"); cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("reference : %s  %d samples x %d features\n", basename(ref.file),
            nrow(A0), ncol(A0) - 1L))
cat(sprintf("ours      : %s  %d samples x %d features\n", basename(our.file),
            nrow(B0), ncol(B0) - 1L))
cat(sprintf("compared  : %d samples x %d features\n", length(samples), length(features)))
only.a <- setdiff(A0$sample, B0$sample); only.b <- setdiff(B0$sample, A0$sample)
if (length(only.a)) cat("  reference-only samples:", length(only.a), "\n")
if (length(only.b)) cat("  ours-only samples     :", length(only.b), "\n")
fa <- setdiff(setdiff(colnames(A0), "sample"), features)
fb <- setdiff(setdiff(colnames(B0), "sample"), features)
if (length(fa)) cat("  reference-only features:", paste(fa, collapse = ", "), "\n")
if (length(fb)) cat("  ours-only features     :", paste(fb, collapse = ", "), "\n")

na.a <- is.na(A); na.b <- is.na(B)
both <- !na.a & !na.b
agree.val <- both & A == B

cat("\n", strrep("=", 74), "\n", sep = "")
cat("WHAT \"DIFFERENT\" MEANS HERE\n"); cat(strrep("=", 74), "\n", sep = "")
n <- length(A)
row <- function(lab, k) cat(sprintf("  %-38s %9s  %6.2f%%\n", lab,
                                    format(k, big.mark = ","), 100 * k / n))
row("both have a value, values agree",   sum(agree.val))
row("both NA",                            sum(na.a & na.b))
row("both have a value, values DIFFER",   sum(both & A != B))
row("reference NA, ours has a value",     sum(na.a & !na.b))
row("ours NA, reference has a value",     sum(!na.a & na.b))
cat(sprintf("  %-38s %9s\n", "total cells", format(n, big.mark = ",")))

cat("\n", strrep("=", 74), "\n", sep = "")
cat("AGREEMENT ON COMPARABLE CELLS (both sides carry a value)\n")
cat(strrep("=", 74), "\n", sep = "")
nc <- sum(both); ac <- sum(agree.val)
cat(sprintf("  %s of %s identical  (%.2f%%)\n", format(ac, big.mark = ","),
            format(nc, big.mark = ","), 100 * ac / nc))
d <- both & A != B
cat(sprintf("  reference higher: %s   ours higher: %s\n",
            format(sum(d & A > B), big.mark = ","),
            format(sum(d & B > A), big.mark = ",")))

per.f <- data.frame(
  feature    = features,
  comparable = colSums(both),
  identical  = colSums(agree.val),
  differing  = colSums(d),
  na_mismatch = colSums(xor(na.a, na.b)),
  stringsAsFactors = FALSE)
per.f$pct <- ifelse(per.f$comparable > 0,
                    round(100 * per.f$identical / per.f$comparable, 2), NA_real_)

cat(sprintf("\n  features fully identical: %d of %d\n",
            sum(per.f$differing == 0), nrow(per.f)))

worst <- per.f[per.f$differing > 0, ]
if (nrow(worst)) {
  worst <- worst[order(-worst$differing), ]
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("FEATURES THAT DISAGREE (worst first)\n"); cat(strrep("=", 74), "\n", sep = "")
  print(head(worst, 30), row.names = FALSE)
}

per.s.n <- rowSums(both); per.s.d <- rowSums(d)
ok <- per.s.n > 0
cat("\n", strrep("=", 74), "\n", sep = "")
cat("PER SAMPLE\n"); cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  fully identical : %d of %d (%.1f%%)\n", sum(per.s.d[ok] == 0), sum(ok),
            100 * mean(per.s.d[ok] == 0)))
cat(sprintf("  median agreement: %.2f%%\n",
            100 * median((per.s.n[ok] - per.s.d[ok]) / per.s.n[ok])))
cat(sprintf("  below 90%%       : %d\n",
            sum((per.s.n[ok] - per.s.d[ok]) / per.s.n[ok] < 0.90)))
o <- sort(per.s.d, decreasing = TRUE)
if (sum(o) > 0) {
  cat(sprintf("  top 50 samples hold %.1f%% of all disagreement\n",
              100 * sum(head(o, 50)) / sum(o)))
}

if (!is.na(csv.out)) {
  write.csv(per.f, csv.out, row.names = FALSE)
  cat("\nper-feature detail written to", csv.out, "\n")
}
cat("\n")
