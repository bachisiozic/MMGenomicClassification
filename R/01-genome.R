#### . . . - Genome registry - . . . ####
#
# One place that knows, for each supported build:
#   * which reference files to use,
#   * which BSgenome / TxDb packages to load,
#   * how chromosomes are numbered,
#   * which gene symbols changed name between builds.
#
# Everything else in the pipeline asks this registry instead of hardcoding a
# build, which is what lets a single set of functions serve hg19 and hg38.

.MMGC <- new.env(parent = emptyenv())
.MMGC$genome <- NULL
.MMGC$refdir <- NULL
.MMGC$refs   <- NULL


#### . . . - Reference file names, per build.
#### . . . - Keys are stable across builds; only the file names differ, so callers
####         can ask for MMGC_ref("cytoband") and not care which build is active.
MMGC_REFERENCE_FILES <- list(
  hg19 = c(
    cytoband     = "cytoBand_hg19.txt",
    gistic_peaks = "Genomic_Coordinates_GISTICpeaks_hg19.txt",
    arms         = "arm_hg19.txt",
    chrom_sizes  = "chrom_sizes_hg19.txt",
    gistic_arms  = "arm_level_hg19.txt",
    gene_models  = "UCSC_hg19_gene_models.txt",
    drivers      = "Driver_summary_hg19.txt"
  ),
  hg38 = c(
    cytoband     = "cytoBand_hg38.txt",
    gistic_peaks = "Genomic_Coordinates_GISTICpeaks_hg38.txt",
    arms         = "arm_hg38.txt",
    chrom_sizes  = "chrom_sizes_hg38.txt",
    gistic_arms  = "arm_level_hg38.txt",
    gene_models  = "UCSC_hg38_gene_models.txt",
    drivers      = "Driver_summary_hg38.txt"
  ),
  mm10 = c(
    cytoband     = "cytoBand_mm10.txt",
    gistic_peaks = "Genomic_Coordinates_GISTICpeaks_mm10.txt",
    arms         = "arm_mm10.txt",
    chrom_sizes  = "chrom_sizes_mm10.txt",
    gene_models  = "UCSC_mm10_gene_models.txt",
    drivers      = "Driver_summary_mm10.txt",
    hrd_map      = "hrd_chromosome_map.txt"
  )
)


#### . . . - Chromosome numbering, per build.
#### . . . - inputCNV() turns chromosome names into integers so that GRanges
####         merges line up. Human: 1-22, X = 23, Y = 24. Mouse: 1-19, X = 20,
####         Y = 21.
MMGC_CHROM_CODING <- list(
  hg19 = list(autosomes = 1:22, X = 23L, Y = 24L),
  hg38 = list(autosomes = 1:22, X = 23L, Y = 24L),
  mm10 = list(autosomes = 1:19, X = 20L, Y = 21L)
)


#### . . . - Gene symbols that were renamed between hg19 and hg38.
#### . . . - The hg19 reference files and any hg19-era mutation calls use the left
####         name; the hg38 files use the right one. The pipeline works internally
####         in hg38 nomenclature and translates hg19 input on the way in, so the
####         feature vectors and output column names are the same for both builds.
####         (The original code carried two hand-edited copies of every one of
####         those vectors instead.)
MMGC_SYMBOL_ALIASES <- c(
  FAM46C    = "TENT5C",
  HIST1H1E  = "H1.4",
  HIST1H1C  = "H1.2",
  HIST1H1D  = "H1.3",
  HIST1H1B  = "H1.5",
  HIST1H2BK = "H2BC12",
  TBC1D29   = "TBC1D29P"
)


#### . . . - MMGC_setup() - . . . ####

#### . . . - Choose the genome build and make its references available.
####         GENOME        = "hg19", "hg38" or "mm10"
####         REFERENCE.DIR = folder holding the per-build reference subfolders.
####                         Defaults to reference/ inside the repository.
####         LOAD.GENOME.PACKAGES = attach the BSgenome/TxDb packages for this
####                         build. Only APOBECanalyzer() needs them, so set this
####                         to FALSE to skip a slow load when running the
####                         copy-number half of the pipeline on its own.
####         QUIET         = suppress the summary message.

MMGC_setup <- function(GENOME,
                       REFERENCE.DIR = NULL,
                       LOAD.GENOME.PACKAGES = TRUE,
                       QUIET = FALSE) {

  if (missing(GENOME) || length(GENOME) != 1L || is.na(GENOME)) {
    stop("MMGC_setup(): give one genome build -- ",
         paste(sQuote(names(MMGC_REFERENCE_FILES)), collapse = ", "), call. = FALSE)
  }
  GENOME <- as.character(GENOME)
  if (!GENOME %in% names(MMGC_REFERENCE_FILES)) {
    stop("MMGC_setup(): unsupported genome ", sQuote(GENOME), ". Supported: ",
         paste(sQuote(names(MMGC_REFERENCE_FILES)), collapse = ", "), call. = FALSE)
  }

  if (is.null(REFERENCE.DIR)) {
    REFERENCE.DIR <- file.path(
      if (exists("MMGC_HOME", envir = globalenv())) get("MMGC_HOME", envir = globalenv()) else getwd(),
      "reference")
  }
  if (!dir.exists(REFERENCE.DIR)) {
    stop("MMGC_setup(): reference directory not found: ", REFERENCE.DIR,
         "\n  Pass REFERENCE.DIR = \"/path/to/reference\" if it lives elsewhere.",
         call. = FALSE)
  }

  build.dir <- file.path(REFERENCE.DIR, GENOME)
  if (!dir.exists(build.dir)) {
    stop("MMGC_setup(): no reference folder for ", GENOME, " at ", build.dir, call. = FALSE)
  }

  wanted <- MMGC_REFERENCE_FILES[[GENOME]]
  paths  <- file.path(build.dir, wanted)
  names(paths) <- names(wanted)

  absent <- names(paths)[!file.exists(paths)]
  if (length(absent)) {
    msg <- paste(sprintf("  %-13s %s", absent, unname(paths[absent])), collapse = "\n")
    stop("MMGC_setup(): reference file(s) missing for ", GENOME, ":\n", msg,
         if (GENOME == "mm10" && any(c("arms", "hrd_map") %in% absent))
           "\n  Build the mouse-specific tables with:\n    Rscript scripts/build_mm10_arms_and_hrd_map.R"
         else
           "\n  Build the gene models with:\n    Rscript scripts/build_slim_gene_references.R <dir-with-full-UCSC-dumps>",
         call. = FALSE)
  }

  .MMGC$genome <- GENOME
  .MMGC$refdir <- normalizePath(REFERENCE.DIR)
  .MMGC$refs   <- paths

  if (isTRUE(LOAD.GENOME.PACKAGES)) {
    for (p in MMGC_PACKAGES_GENOME[[GENOME]]) {
      if (!suppressPackageStartupMessages(
            require(p, character.only = TRUE, quietly = TRUE))) {
        warning("MMGC_setup(): ", p, " is not installed. APOBECanalyzer() will ",
                "not run for ", GENOME, " until it is.", call. = FALSE)
      }
    }
  }

  if (!isTRUE(QUIET)) {
    message("MMGenomicClassification: genome set to ", GENOME)
    message("  references: ", build.dir)
    if (GENOME == "mm10") {
      message("  NOTE mouse build. Hyperdiploidy has no mouse equivalent; HRD ",
              "features default to 0.")
      message("       See reference/mm10/hrd_chromosome_map.txt and docs/MOUSE.md.")
    }
  }

  invisible(paths)
}


#### . . . - MMGC_genome() - . . . ####
#### . . . - The build currently in force. Errors if MMGC_setup() has not run.

MMGC_genome <- function() {
  if (is.null(.MMGC$genome)) {
    stop("No genome build selected yet. Call MMGC_setup(\"hg38\") ",
         "(or \"hg19\" / \"mm10\") first.", call. = FALSE)
  }
  .MMGC$genome
}


#### . . . - MMGC_ref() - . . . ####
#### . . . - Path to a reference file for the active build.
####         KEY = one of the names in MMGC_REFERENCE_FILES, e.g. "cytoband".
####         Called with no arguments, returns every path as a named vector.

MMGC_ref <- function(KEY) {
  genome <- MMGC_genome()
  if (missing(KEY)) return(.MMGC$refs)
  if (length(KEY) != 1L) stop("MMGC_ref(): ask for one key at a time.", call. = FALSE)
  if (!KEY %in% names(.MMGC$refs)) {
    stop("MMGC_ref(): ", sQuote(KEY), " is not a reference for build ", genome,
         ". Available: ", paste(sQuote(names(.MMGC$refs)), collapse = ", "),
         call. = FALSE)
  }
  unname(.MMGC$refs[[KEY]])
}


#### . . . - MMGC_is_mouse() - . . . ####

MMGC_is_mouse <- function() identical(MMGC_genome(), "mm10")


#### . . . - MMGC_canonical_symbols() - . . . ####
#### . . . - Rewrite hg19-era gene symbols to their current (hg38) names, so that
####         one set of feature vectors serves both builds.
####         X      = character vector of symbols, or of column names.
####         PREFIX = TRUE when X holds prefixed column names ("CNV.SNV_FAM46C",
####                  "SNV_HIST1H1E"); the prefix is preserved.

MMGC_canonical_symbols <- function(X, PREFIX = FALSE) {
  if (!length(X)) return(X)
  aliases <- MMGC_SYMBOL_ALIASES

  if (!isTRUE(PREFIX)) {
    hit <- X %in% names(aliases)
    X[hit] <- unname(aliases[X[hit]])
    return(X)
  }

  # Column names look like <prefix>_<symbol>; split on the LAST underscore that
  # separates a known prefix from the symbol.
  for (old in names(aliases)) {
    pattern <- paste0("(^|_)", old, "$")
    X <- sub(pattern, paste0("\\1", aliases[[old]]), X)
  }
  X
}
