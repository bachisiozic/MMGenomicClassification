#### . . . - Input QC for copy-number calls - . . . ####
#
# The pipeline is built for whole-genome sequencing: allele-specific, purity- and
# ploidy-corrected segments from ASCAT, Battenberg, FACETS or equivalent. Given
# that, the published gain definition (`total >= 2.5 & minor >= 0.5`, i.e.
# REQUIRE.MINOR.ALLELE = TRUE) is the right one, because the minor allele is
# available for every segment.
#
# What goes wrong is silent. Two failure modes cost real samples without raising
# an error anywhere:
#
#   * No minor allele. A segment with `minor = NA` can never satisfy
#     `minor >= 0.5`, so it never counts towards a gain. A sample whose calls are
#     total-copy-number only therefore comes back as non-hyperdiploid, not as
#     unassessable. This is what array or panel data looks like when it is mixed
#     into a WGS cohort.
#
#   * Overlapping segments. HYidentification() sums gained bases per arm; when
#     segments overlap, the same bases are counted more than once and the gained
#     fraction is inflated, which can push a chromosome over the 60% threshold on
#     far less than 60% of real coverage. Allele-specific WGS callers emit disjoint
#     segments, so overlap means the input was merged from more than one source or
#     not collapsed.
#
# checkCNV() looks for both, plus the structural problems that make thresholds
# meaningless, and says which samples are safe to analyse.


#### . . . - .mmgc_weighted_median() - . . . ####
#### . . . - Length-weighted median, ignoring segments with a missing value. Used
####         for the ploidy estimate, where a plain median over segments would be
####         dominated by the many short ones.

.mmgc_weighted_median <- function(X, W) {
  keep <- !is.na(X) & !is.na(W) & W > 0
  if (!any(keep)) return(NA_real_)
  X <- X[keep]; W <- W[keep]
  o <- order(X)
  cw <- cumsum(W[o]) / sum(W[o])
  X[o][which(cw >= 0.5)[1L]]
}


#### . . . - checkCNV() - . . . ####

#### . . . - Check a copy-number table for WGS analysis, before running anything.
####
####         CNV     = output of inputCNV()
####         GENOME  = build to size the genome against; defaults to MMGC_setup()'s
####         MAX.COVERAGE = flag samples whose segments cover more than this
####                   fraction of the autosomes (default 1.05; above 1 means the
####                   same base is in more than one segment)
####         MIN.COVERAGE = flag samples covering less than this (default 0.8)
####         PLOIDY.RANGE = plausible length-weighted median total copy number
####                   (default 1.2 to 5). Outside it, the calls were almost
####                   certainly not normalised for purity and ploidy, and every
####                   absolute threshold downstream is meaningless.
####         REPORT  = print the summary
####
####         Returns, invisibly, a list:
####           $per_sample  one row per sample, with every metric and a verdict
####           $problems    the samples that need attention
####           $advice      the recommended REQUIRE.MINOR.ALLELE for this cohort

checkCNV <- function(CNV,
                     GENOME = NULL,
                     MAX.COVERAGE = 1.05,
                     MIN.COVERAGE = 0.80,
                     PLOIDY.RANGE = c(1.2, 5),
                     REPORT = TRUE) {

  genome <- if (is.null(GENOME)) MMGC_genome() else GENOME
  coding <- MMGC_CHROM_CODING[[genome]]

  cnv <- as.data.frame(CNV, stringsAsFactors = FALSE)
  need <- c("sample", "chrom", "start", "end", "total", "minor")
  absent <- setdiff(need, colnames(cnv))
  if (length(absent)) {
    stop("checkCNV(): missing column(s): ", paste(absent, collapse = ", "),
         "\n  Run the table through inputCNV() first.", call. = FALSE)
  }

  sizes <- .mmgc_chrom_sizes()
  autosomes <- coding$autosomes
  autosome.bp <- sum(sizes$Size[sizes$Offset %in% autosomes])

  cnv <- cnv[order(cnv$sample, cnv$chrom, cnv$start), , drop = FALSE]
  cnv$width <- cnv$end - cnv$start

  #### . . . - Overlap: within a sample and chromosome, does a segment start before
  ####         the previous one ends?
  key <- paste(cnv$sample, cnv$chrom, sep = "\r")
  prev.end <- c(NA_real_, cnv$end[-nrow(cnv)])
  same.grp <- c(FALSE, key[-1L] == key[-length(key)])
  cnv$overlaps <- same.grp & !is.na(prev.end) & cnv$start < prev.end

  by.sample <- split(seq_len(nrow(cnv)), cnv$sample)

  per <- do.call(rbind, lapply(names(by.sample), function(s) {
    d <- cnv[by.sample[[s]], , drop = FALSE]
    a <- d[d$chrom %in% autosomes, , drop = FALSE]
    data.frame(
      sample             = s,
      n_segments         = nrow(d),
      frac_minor_missing = mean(is.na(d$minor)),
      n_overlaps         = sum(d$overlaps),
      autosome_coverage  = if (nrow(a)) sum(pmax(0, a$width, na.rm = TRUE)) / autosome.bp else 0,
      #### . . . - Length-weighted median total copy number: a ploidy sanity check.
      ploidy             = .mmgc_weighted_median(a$total, a$width),
      n_zero_width       = sum(d$width <= 0, na.rm = TRUE),
      n_negative_cn      = sum(d$total < 0, na.rm = TRUE) +
                           sum(d$minor < 0, na.rm = TRUE),
      n_duplicated       = sum(duplicated(d[, c("chrom", "start", "end")])),
      #### . . . - Reported, not treated as an error: with fractional (subclonal or
      ####         purity-scaled) copy number, minor > total - minor is routine --
      ####         total = 1.5 with minor = 1 is a subclonal loss, not a malformed
      ####         call. Nothing in the pipeline uses the major allele, so this
      ####         cannot break a threshold; it is here only as a caller-convention
      ####         signal.
      n_minor_gt_major   = sum(d$minor > d$total - d$minor + 1e-9, na.rm = TRUE),
      stringsAsFactors   = FALSE)
  }))
  rownames(per) <- NULL

  #### . . . - Independent flags, then one verdict per sample: the most serious
  ####         thing found, worst first. A sample can trip several.
  per$allele_specific <- per$frac_minor_missing == 0
  per$malformed  <- per$n_zero_width + per$n_negative_cn + per$n_duplicated > 0
  per$ploidy_ok  <- !is.na(per$ploidy) &
                      per$ploidy >= PLOIDY.RANGE[1] & per$ploidy <= PLOIDY.RANGE[2]
  per$overlapping <- per$n_overlaps > 0 | per$autosome_coverage > MAX.COVERAGE
  per$sparse      <- per$autosome_coverage < MIN.COVERAGE

  per$verdict <- ifelse(
    per$malformed,                   "malformed segments",
    ifelse(!per$ploidy_ok,           "implausible ploidy",
    ifelse(per$frac_minor_missing == 1, "no allele-specific data",
    ifelse(per$overlapping,          "overlapping segments",
    ifelse(per$frac_minor_missing > 0, "partial allele-specific data",
    ifelse(per$sparse,               "sparse coverage", "ok"))))))

  problems <- per[per$verdict != "ok", , drop = FALSE]
  problems <- problems[order(problems$verdict, problems$sample), , drop = FALSE]

  #### . . . - Cohort-level recommendation.
  n.no.allelic <- sum(per$frac_minor_missing == 1)
  n.partial    <- sum(per$frac_minor_missing > 0 & per$frac_minor_missing < 1)
  #### . . . - TRUE is right for allele-specific WGS, which is what this pipeline
  ####         expects. Only recommend relaxing it when the data cannot support it.
  advice <- if (n.no.allelic == 0) TRUE else "non-HRD"

  if (isTRUE(REPORT)) {
    line <- function() cat(strrep("=", 74), "\n", sep = "")
    cat("\n"); line()
    cat("checkCNV(): ", nrow(per), " samples, ", format(nrow(cnv), big.mark = ","),
        " segments, ", genome, "\n", sep = "")
    line()

    #### . . . - One field width for every counted line, so the report scans.
    kv <- function(label, value) cat(sprintf("  %-44s: %s\n", label, value))

    cat("\nAllele-specific copy number\n")
    kv("complete (minor allele on every segment)", sum(per$frac_minor_missing == 0))
    kv("partial", n.partial)
    kv("absent (total copy number only)", n.no.allelic)

    cat("\nSegment structure\n")
    kv("samples with overlapping segments", sum(per$n_overlaps > 0))
    kv(sprintf("samples covering >%.0f%% of the autosomes", 100 * MAX.COVERAGE),
       sum(per$autosome_coverage > MAX.COVERAGE))
    kv(sprintf("samples covering <%.0f%% of the autosomes", 100 * MIN.COVERAGE),
       sum(per$autosome_coverage < MIN.COVERAGE))
    kv("samples with malformed segments", sum(per$malformed))
    if (sum(per$n_minor_gt_major) > 0) {
      kv("segments with minor > major",
         sprintf("%s in %d sample(s)", format(sum(per$n_minor_gt_major), big.mark = ","),
                 sum(per$n_minor_gt_major > 0)))
      cat("      not an error -- routine for subclonal / purity-scaled calls, and\n")
      cat("      the pipeline never uses the major allele\n")
    }

    cat("\nPloidy (length-weighted median total copy number)\n")
    pl <- per$ploidy[!is.na(per$ploidy)]
    if (length(pl)) {
      kv("median", sprintf("%.2f", stats::median(pl)))
      kv("range", sprintf("%.2f to %.2f", min(pl), max(pl)))
      kv(sprintf("outside the plausible range %.1f - %.1f",
                 PLOIDY.RANGE[1], PLOIDY.RANGE[2]), sum(!per$ploidy_ok))
    }

    cat("\nVerdict -- the most serious finding per sample. The checks above are\n")
    cat("independent and a sample can trip several; here it appears once.\n")
    v <- table(per$verdict)
    for (nm in names(sort(v, decreasing = TRUE))) kv(nm, v[[nm]])

    if (nrow(problems)) {
      cat("\n"); line()
      cat("WHAT TO DO\n"); line()

      n.ploidy <- sum(!per$ploidy_ok)
      if (n.ploidy) {
        bad <- per[!per$ploidy_ok, ]
        bad <- bad[order(-abs(bad$ploidy - 2)), ]
        cat("\n", n.ploidy, " sample(s) have an implausible ploidy (outside ",
            PLOIDY.RANGE[1], " - ", PLOIDY.RANGE[2], ").\n", sep = "")
        cat("  Every threshold in the pipeline is an ABSOLUTE copy number -- gain at\n")
        cat("  total >= 2.5, loss below 1.5 -- so calls that were not corrected for\n")
        cat("  purity and ploidy are scored genome-wide. A sample sitting at total 20\n")
        cat("  is gained on every chromosome and comes out hyperdiploid automatically.\n")
        cat("  Worst: ",
            paste(sprintf("%s (ploidy %.0f)", utils::head(bad$sample, 5),
                          utils::head(bad$ploidy, 5)), collapse = ", "), "\n", sep = "")
        cat("  Re-call or drop these before running.\n")
      }

      if (n.no.allelic) {
        cat("\n", n.no.allelic, " sample(s) have NO allele-specific data at all.\n", sep = "")
        cat("  With the default REQUIRE.MINOR.ALLELE = TRUE no segment can satisfy\n")
        cat("  `minor >= 0.5`, so every one of them is reported as non-hyperdiploid\n")
        cat("  rather than as unassessable. Either\n")
        cat("    - re-call them with an allele-specific WGS caller (ASCAT, Battenberg,\n")
        cat("      FACETS), which is what this pipeline is built for; or\n")
        cat("    - drop them; or\n")
        cat("    - pass HYidentification(REQUIRE.MINOR.ALLELE = \"non-HRD\"), which scores\n")
        cat("      trisomies on total copy number and keeps gains elsewhere\n")
        cat("      allele-aware. See docs/KNOWN_ISSUES.md.\n")
      }
      if (n.partial) {
        pf <- per$frac_minor_missing[per$frac_minor_missing > 0 &
                                       per$frac_minor_missing < 1]
        cat("\n", n.partial, " sample(s) are missing the minor allele on SOME segments ",
            "(median ", sprintf("%.1f%%", 100 * stats::median(pf)), ",\n", sep = "")
        cat("  worst ", sprintf("%.1f%%", 100 * max(pf)), ").\n", sep = "")
        cat("  Those segments cannot contribute to a gain, so a trisomy is missed\n")
        cat("  whenever the gained span is one of them. Usually a handful of\n")
        cat("  low-confidence regions; check the worst samples before trusting their\n")
        cat("  hyperdiploidy calls.\n")
      }

      n.ov <- sum(per$n_overlaps > 0 | per$autosome_coverage > MAX.COVERAGE)
      if (n.ov) {
        cat("\n", n.ov, " sample(s) have overlapping segments.\n", sep = "")
        cat("  Gained bases are summed per arm, so overlap double-counts them and\n")
        cat("  inflates the gained fraction -- a chromosome can pass the 60% threshold\n")
        cat("  on much less than 60% of real coverage. Allele-specific WGS callers emit\n")
        cat("  disjoint segments, so this usually means the table was merged from more\n")
        cat("  than one source. Resolve the overlaps before running.\n")
      }

      n.sparse <- sum(per$autosome_coverage < MIN.COVERAGE)
      if (n.sparse) {
        cat("\n", n.sparse, " sample(s) cover less than ", round(100 * MIN.COVERAGE),
            "% of the autosomes.\n", sep = "")
        cat("  Gained fractions are computed against full chromosome length, so\n")
        cat("  uncovered regions count as not gained and trisomies are under-called.\n")
      }

      n.mal <- sum(per$n_zero_width + per$n_negative_cn + per$n_duplicated > 0)
      if (n.mal) {
        cat("\n", n.mal, " sample(s) have malformed segments ",
            "(zero or negative width,\n", sep = "")
        cat("  negative copy number, or exact duplicate coordinates). ",
            sum(per$n_zero_width), " zero-width,\n  ", sep = "")
        cat(sum(per$n_negative_cn), " negative-CN, ", sum(per$n_duplicated),
            " duplicated segments in total.\n", sep = "")
      }
    }

    cat("\n"); line()
    cat("Recommended setting for this cohort: REQUIRE.MINOR.ALLELE = ",
        if (isTRUE(advice)) "TRUE (the default)" else paste0("\"", advice, "\""), "\n", sep = "")
    if (isTRUE(advice)) {
      cat("  Every sample has allele-specific calls, so the published gain definition\n")
      cat("  applies as intended.\n")
    }
    cat(sprintf("Samples with no findings: %d of %d\n", sum(per$verdict == "ok"), nrow(per)))
    line(); cat("\n")
  }

  invisible(list(per_sample = per, problems = problems, advice = advice))
}
