#### . . . - Dependencies - . . . ####
#
# The original function files opened with 25 library() calls. Ten of those
# packages are never actually used by any function in the pipeline (gtools, plyr,
# ComplexHeatmap, VariantAnnotation, Biostrings, annotate, org.Hs.eg.db,
# clonevol, survival, survminer), so they are not required here. pheatmap is not
# required either: MMGenoCla() returns the matrix, annotation and colour list for
# pheatmap::pheatmap() but never calls it, so drawing the heatmap is the caller's
# job -- see examples/.

# Packages every build needs, attached because the pipeline code calls them
# unqualified (%>%, arrange, setDT, GRanges, str_split_fixed, ...).
MMGC_PACKAGES_CORE <- c(
  "dplyr",         # arrange(), distinct(), %>%
  "data.table",    # setDT(), dcast()
  "GenomicRanges", # GRanges(), values()
  "IRanges",       # IRanges()
  "stringr",       # str_split_fixed()
  "stringi",       # stri_sub()
  "RColorBrewer"   # brewer.pal()
)

# Needed only by APOBECanalyzer(). Kept separate so the copy-number half of the
# pipeline runs on a machine without the signature-fitting stack installed.
MMGC_PACKAGES_SIGNATURES <- c(
  "mmsig",              # mm_fit_signatures() and its internals
  "deconstructSigs",    # mut.to.sigs.input()
  "MutationalPatterns", # mut_matrix_stranded(), strand_bias_test()
  "reshape2",           # melt(), dcast()
  "tibble",             # rownames_to_column()
  "tidyr"               # separate()
)

# Genome packages, per build. mmsig hardcodes hg19 internally, so the pipeline
# supplies these explicitly -- see R/09-APOBECanalyzer.R.
MMGC_PACKAGES_GENOME <- list(
  hg19 = c(bsgenome = "BSgenome.Hsapiens.UCSC.hg19",
           txdb     = "TxDb.Hsapiens.UCSC.hg19.knownGene"),
  hg38 = c(bsgenome = "BSgenome.Hsapiens.UCSC.hg38",
           txdb     = "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  mm10 = c(bsgenome = "BSgenome.Mmusculus.UCSC.mm10",
           txdb     = "TxDb.Mmusculus.UCSC.mm10.knownGene")
)


#### . . . - MMGC_check_packages() - . . . ####
#### . . . - Report which dependencies are present without loading them.
####         WHICH = "core", "signatures", "genome" or "all"
####         GENOME = build to check the genome packages for; defaults to the
####                  build set by MMGC_setup(), if any.

MMGC_check_packages <- function(WHICH = "all", GENOME = NULL) {
  WHICH <- match.arg(WHICH, c("all", "core", "signatures", "genome"))

  wanted <- character(0)
  if (WHICH %in% c("all", "core"))       wanted <- c(wanted, MMGC_PACKAGES_CORE)
  if (WHICH %in% c("all", "signatures")) wanted <- c(wanted, MMGC_PACKAGES_SIGNATURES)
  if (WHICH %in% c("all", "genome")) {
    build <- if (!is.null(GENOME)) GENOME else tryCatch(MMGC_genome(), error = function(e) NULL)
    if (!is.null(build)) wanted <- c(wanted, unname(MMGC_PACKAGES_GENOME[[build]]))
  }
  wanted <- unique(wanted)

  ok <- vapply(wanted, function(p) requireNamespace(p, quietly = TRUE), logical(1))
  out <- data.frame(package = wanted, installed = unname(ok), stringsAsFactors = FALSE)

  if (any(!ok)) {
    message("Missing package(s): ", paste(out$package[!ok], collapse = ", "))
    message("Install with: Rscript scripts/install_dependencies.R")
  }
  out
}


#### . . . - Attach the core packages now -- these are needed by every function.
local({
  missing <- character(0)
  for (p in MMGC_PACKAGES_CORE) {
    if (!suppressPackageStartupMessages(require(p, character.only = TRUE, quietly = TRUE))) {
      missing <- c(missing, p)
    }
  }
  if (length(missing)) {
    stop("MMGenomicClassification needs these packages, which are not installed: ",
         paste(missing, collapse = ", "),
         "\n  Install them with: Rscript scripts/install_dependencies.R",
         call. = FALSE)
  }
})
