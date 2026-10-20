#!/usr/bin/env Rscript
#
# build_slim_gene_references.R
#
# Produce the slim UCSC gene-model references that ship with this repository.
#
# The pipeline only ever reads four columns out of the UCSC gene tables
# (see MutationsMatrixGenerator(): chromosome, transcription start, transcription
# end, gene symbol). The full UCSC dumps are large -- the hg38 table is 78 MB,
# which is past GitHub's 50 MB advisory limit -- so the repository ships the
# four-column projection instead, under one unified set of column names:
#
#     chrom   txStart   txEnd   geneSymbol
#
# Row order is preserved exactly. This matters: MutationsMatrixGenerator() calls
# dplyr::distinct(geneSymbol, .keep_all = TRUE), which keeps the FIRST transcript
# listed for each gene symbol. Re-sorting the table would silently change which
# transcript's coordinates are used for every gene.
#
# Usage:
#   Rscript scripts/build_slim_gene_references.R <source-dir>
#
# where <source-dir> holds the original full UCSC dumps:
#   <source-dir>/hg19/UCSC_hg19_ref_genes.dms
#   <source-dir>/hg38/UCSC_hg38_ref_genes.dms
#   <source-dir>/mm10/UCSC_mm10_ref_genes.dms
#
# Outputs are written to reference/<build>/UCSC_<build>_gene_models.txt.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript scripts/build_slim_gene_references.R <source-dir>", call. = FALSE)
}
src.dir <- normalizePath(args[[1]], mustWork = TRUE)

# Resolve the repository root relative to this script, not to getwd().
this.file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
repo.root <- if (length(this.file)) {
  dirname(dirname(normalizePath(this.file)))
} else {
  getwd()
}

# The column that carries each field, per source table. UCSC has renamed these
# over the years, so each build is listed explicitly rather than guessed.
column.map <- list(
  hg19 = c(chrom      = "UCSC.hg19.knownGene.chrom",
           txStart    = "UCSC.hg19.knownGene.txStart",
           txEnd      = "UCSC.hg19.knownGene.txEnd",
           geneSymbol = "UCSC.hg19.kgXref.geneSymbol"),
  hg38 = c(chrom      = "chrom",
           txStart    = "txStart",
           txEnd      = "txEnd",
           geneSymbol = "geneSymbol"),
  mm10 = c(chrom      = "chrom",
           txStart    = "txStart",
           txEnd      = "txEnd",
           geneSymbol = "geneSymbol")
)

source.file <- c(hg19 = "hg19/UCSC_hg19_ref_genes.dms",
                 hg38 = "hg38/UCSC_hg38_ref_genes.dms",
                 mm10 = "mm10/UCSC_mm10_ref_genes.dms")

for (build in names(column.map)) {
  infile <- file.path(src.dir, source.file[[build]])
  if (!file.exists(infile)) {
    message(sprintf("-- %s: source table not found, skipping (%s)", build, infile))
    next
  }

  wanted <- column.map[[build]]
  full <- fread(infile, sep = "\t", header = TRUE, quote = "",
                select = unname(wanted), colClasses = "character",
                showProgress = FALSE)

  missing <- setdiff(unname(wanted), names(full))
  if (length(missing)) {
    stop(sprintf("%s: expected column(s) absent from %s: %s",
                 build, infile, paste(missing, collapse = ", ")), call. = FALSE)
  }

  # fread(select=) does not guarantee the requested order; impose it, then rename.
  full <- full[, unname(wanted), with = FALSE]
  setnames(full, unname(wanted), names(wanted))

  outdir <- file.path(repo.root, "reference", build)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(outdir, sprintf("UCSC_%s_gene_models.txt", build))

  fwrite(full, outfile, sep = "\t", quote = FALSE, na = "NA")

  message(sprintf("-- %s: %s rows, %s genes -> %s (%.1f MB, was %.1f MB)",
                  build,
                  format(nrow(full), big.mark = ","),
                  format(length(unique(full$geneSymbol)), big.mark = ","),
                  file.path("reference", build, basename(outfile)),
                  file.size(outfile) / 1024^2,
                  file.size(infile) / 1024^2))
}
