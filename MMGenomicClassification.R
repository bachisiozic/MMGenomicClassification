#####################################################################################
####                                                                             ####
####   MMGenomicClassification -- entry point                                    ####
####                                                                             ####
####   Genomic classification of multiple myeloma, as published in:              ####
####     Maura F, Rustad EH, Boyle EM, et al.                                    ####
####     J Clin Oncol. doi:10.1200/JCO.23.01277                                  ####
####                                                                             ####
####   Source this one file, then pick a genome build:                           ####
####                                                                             ####
####       source("MMGenomicClassification.R")                                   ####
####       MMGC_setup("hg38")        # or "hg19", or "mm10"                      ####
####                                                                             ####
####   Every pipeline function then resolves its own reference files from the     ####
####   build you chose, so paths no longer have to be typed out by hand. Passing  ####
####   a path explicitly still overrides the default, so scripts written against  ####
####   the original function files keep working.                                  ####
####                                                                             ####
#####################################################################################

MMGC_HOME <- local({
  # Where this file lives, however it was reached: source()d from the console,
  # source()d from another script, or run with Rscript. source() records the path
  # as `ofile` in its own evaluation frame, so look for the innermost frame that
  # has one; fall back to --file=, then to getwd().
  this <- NA_character_
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frames()[[i]]$ofile
    if (!is.null(of) && !is.na(of)) { this <- as.character(of); break }
  }
  if (is.na(this)) {
    argv <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", grep("^--file=", argv, value = TRUE))
    if (length(f)) this <- f[[1]]
  }
  if (is.na(this)) getwd() else dirname(normalizePath(this))
})

local({
  files <- c(
    "00-dependencies.R",
    "01-genome.R",
    "02-utils.R",
    "03-inputCNV.R",
    "04-SegmentCollapsing.R",
    "05-FocalLargeMatrix.R",
    "06-HYidentification.R",
    "07-CytogeneticMatrixGenerator.R",
    "08-MutationsMatrixGenerator.R",
    "09-APOBECanalyzer.R",
    "10-FinalMatrixGenerator.R",
    "11-dis3_hs.R",
    "12-MMGenoCla.R",
    "13-checkCNV.R"
  )
  for (f in files) {
    path <- file.path(MMGC_HOME, "R", f)
    if (!file.exists(path)) {
      stop("MMGenomicClassification: cannot find ", path,
           "\n  Source MMGenomicClassification.R from inside the cloned repository.",
           call. = FALSE)
    }
    source(path, local = FALSE)
  }
})

message("MMGenomicClassification loaded from ", MMGC_HOME)
message("Next: MMGC_setup(\"hg38\")   # or \"hg19\", or \"mm10\"")
