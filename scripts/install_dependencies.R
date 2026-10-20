#!/usr/bin/env Rscript
#
# install_dependencies.R -- install everything the pipeline needs.
#
# Usage:
#   Rscript scripts/install_dependencies.R            # core + signatures + human genomes
#   Rscript scripts/install_dependencies.R hg38       # core + signatures + hg38 only
#   Rscript scripts/install_dependencies.R mm10       # ... and the mouse packages
#   Rscript scripts/install_dependencies.R --core     # copy-number half only, no genomes
#
# Genome packages are ~1 GB each, so install only the builds you will use.

args <- commandArgs(trailingOnly = TRUE)
core.only <- "--core" %in% args
builds <- setdiff(args, "--core")
if (!length(builds) && !core.only) builds <- c("hg19", "hg38")

cran <- c("dplyr", "data.table", "stringr", "stringi", "RColorBrewer",
          "reshape2", "tibble", "tidyr", "pheatmap", "ggplot2",
          "deconstructSigs", "devtools", "BiocManager")

bioc <- c("GenomicRanges", "IRanges", "GenomeInfoDb", "GenomicFeatures",
          "Biostrings", "BSgenome", "MutationalPatterns")

genome.pkgs <- list(
  hg19 = c("BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene"),
  hg38 = c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  mm10 = c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene"))

# Not on CRAN or Bioconductor.
github <- c(mmsig      = "UM-Myeloma-Genomics/mmsig",
            dndscv     = "im3sanger/dndscv")

need <- function(p) !requireNamespace(p, quietly = TRUE)

install_if_needed <- function(pkgs, installer) {
  todo <- pkgs[vapply(pkgs, need, logical(1))]
  if (!length(todo)) { message("  all present"); return(invisible()) }
  message("  installing: ", paste(todo, collapse = ", "))
  installer(todo)
}

message("\n== CRAN packages")
install_if_needed(cran, function(p) install.packages(p, repos = "https://cloud.r-project.org"))

message("\n== Bioconductor packages")
if (need("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
install_if_needed(bioc, function(p) BiocManager::install(p, ask = FALSE, update = FALSE))

if (!core.only) {
  for (b in builds) {
    if (is.null(genome.pkgs[[b]])) {
      message("\n!! unknown build: ", b, " (expected hg19, hg38 or mm10)")
      next
    }
    message("\n== Genome packages for ", b, " (large downloads)")
    install_if_needed(genome.pkgs[[b]],
                      function(p) BiocManager::install(p, ask = FALSE, update = FALSE))
  }
}

message("\n== GitHub packages")
for (nm in names(github)) {
  if (need(nm)) {
    message("  installing ", nm, " from ", github[[nm]])
    devtools::install_github(github[[nm]], upgrade = "never")
  } else {
    message("  ", nm, " present")
  }
}

message("\n== Check")
ok <- TRUE
report <- function(label, pkgs) {
  miss <- pkgs[vapply(pkgs, need, logical(1))]
  if (length(miss)) {
    ok <<- FALSE
    message("  ", label, ": MISSING ", paste(miss, collapse = ", "))
  } else message("  ", label, ": ok")
}
report("core",       c("dplyr", "data.table", "GenomicRanges", "IRanges",
                       "stringr", "stringi", "RColorBrewer"))
report("signatures", c("mmsig", "deconstructSigs", "MutationalPatterns",
                       "reshape2", "tibble", "tidyr"))
report("dndscv",     "dndscv")
if (!core.only) for (b in builds) if (!is.null(genome.pkgs[[b]])) report(b, genome.pkgs[[b]])

message("")
if (ok) {
  message("All set. Next:  Rscript scripts/selftest.R")
} else {
  message("Some packages are missing -- see above. Bioconductor packages sometimes ")
  message("need system libraries; check the install log for the underlying error.")
  quit(status = 1)
}
