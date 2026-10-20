#### . . . - Internal helpers - . . . ####
#
# Three patterns were written out by hand several times each in the original
# function files (interval overlap, collapsing multi-value cells to their
# maximum, dropping columns by pattern). They are factored out here so the
# pipeline functions read as the analysis they describe, and so a fix lands once.
#
# Names start with a dot: these are implementation details, not part of the
# public pipeline API.


#### . . . - .mmgc_overlap() - . . . ####
#### . . . - Every pair of rows from A and B that sit on the same chromosome and
####         overlap by at least one base.
####
####         Reproduces the original four-way range comparison exactly:
####           1  B starts at/before A and ends inside A
####           2  B starts inside A and ends at/after A
####           3  B entirely inside A
####           4  B entirely contains A
####         Their union is "any overlap of non-zero width"; the cases are not
####         disjoint, so duplicated rows are removed as before.
####
####         A, B          = data frames
####         A.CHROM, ...  = column names holding chromosome / start / end
####         A.KEEP, B.KEEP = extra columns to carry through
####
####         Returns a data frame with seqnames, startA/endA, startB/endB and the
####         requested columns -- the same shape the original inline code built.

.mmgc_overlap <- function(A, B,
                          A.CHROM = "chrom", A.START = "start", A.END = "end",
                          B.CHROM = "chr",   B.START = "start", B.END = "end",
                          A.KEEP, B.KEEP) {

  gr.a <- GenomicRanges::GRanges(as.character(A[[A.CHROM]]),
                                 IRanges::IRanges(start = as.numeric(A[[A.START]]),
                                                  end   = as.numeric(A[[A.END]])))
  GenomicRanges::values(gr.a) <- A[, A.KEEP, drop = FALSE]

  gr.b <- GenomicRanges::GRanges(as.character(B[[B.CHROM]]),
                                 IRanges::IRanges(start = as.numeric(B[[B.START]]),
                                                  end   = as.numeric(B[[B.END]])))
  GenomicRanges::values(gr.b) <- B[, B.KEEP, drop = FALSE]

  ranges <- merge(as.data.frame(gr.a), as.data.frame(gr.b),
                  by = "seqnames", suffixes = c("A", "B"))
  if (!nrow(ranges)) return(ranges)

  keep <- with(ranges,
    (startB <= startA & endB <= endA & endB > startA) |
    (startB >= startA & endB >= endA & startB < endA) |
    (startB >= startA & endB <= endA) |
    (startB <  startA & endB >  endA))

  unique(ranges[keep, , drop = FALSE])
}


#### . . . - .mmgc_max_by_cell() - . . . ####
#### . . . - Reduce a sample x feature table whose cells hold comma-separated
####         values ("1, 2, 0", from data.table::dcast(..., toString)) to the
####         maximum value per cell.
####
####         Matches the original nested-loop behaviour, including two edge cases
####         that the loop reached by accident and that real data does hit:
####           * an empty cell (no overlapping segment) is NA, not 0;
####           * a cell whose FIRST value does not parse as a number is NA. This
####             is how the "qqq" sentinel that FocalLargeMatrix() writes for
####             total >= 2 with an unknown minor allele ends up as NA.
####
####         Vectorised rather than looping over every cell: the original walked
####         the table one cell at a time, rebuilding a data frame each step,
####         which dominated the runtime on cohort-sized inputs.

.mmgc_max_by_cell <- function(X) {
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  if (!nrow(X) || !ncol(X)) return(X)

  for (j in seq_len(ncol(X))) {
    cells <- as.character(X[[j]])
    cells[!is.na(cells) & cells == ""] <- NA_character_

    parts <- strsplit(cells, ",", fixed = TRUE)
    X[[j]] <- vapply(parts, function(p) {
      if (!length(p)) return(NA_real_)
      v <- suppressWarnings(as.numeric(trimws(p)))
      # First element unparseable (or the cell was NA) -> NA, as in the original.
      if (is.na(v[1L])) return(NA_real_)
      max(v, na.rm = TRUE)
    }, numeric(1))
  }
  X
}


#### . . . - .mmgc_drop_cols() - . . . ####
#### . . . - Drop columns whose names match a pattern, safely.
####
####         The original code used df[, -grep(pattern, colnames(df))]. When the
####         pattern matches nothing, grep() returns integer(0) and the negative
####         index silently returns a table with ZERO columns instead of leaving
####         the input alone -- so a cohort that happens not to carry the gene
####         being dropped loses its whole matrix. This guards that.

.mmgc_drop_cols <- function(X, PATTERN, FIXED = TRUE) {
  hit <- grep(PATTERN, colnames(X), fixed = FIXED)
  if (!length(hit)) return(X)
  X[, -hit, drop = FALSE]
}


#### . . . - .mmgc_chrom_sizes() - . . . ####
#### . . . - Chromosome sizes for the active build, keyed by the same integer
####         codes inputCNV() assigns (human X = 23, Y = 24; mouse X = 20,
####         Y = 21).
####
####         The original code appended hardcoded rows  c(23, 155270560)  and
####         c(24, 59373566)  -- the hg19 X and Y lengths -- to whichever size
####         file was passed, including the hg38 one. The size files already
####         carry X and Y rows, so those are read instead of being reinvented.
####         Only per-chromosome gain fractions for X/Y were affected, and those
####         feed no HRD call, but there is no reason to keep the wrong numbers.

.mmgc_chrom_sizes <- function(FILE = NULL) {
  if (is.null(FILE)) FILE <- MMGC_ref("chrom_sizes")
  coding <- MMGC_CHROM_CODING[[MMGC_genome()]]

  sizes <- read.delim(FILE, stringsAsFactors = FALSE)
  if (!all(c("Offset", "Size") %in% colnames(sizes))) {
    stop(".mmgc_chrom_sizes(): ", FILE,
         " must have columns 'Offset' and 'Size'.", call. = FALSE)
  }

  off <- gsub("^chr", "", as.character(sizes$Offset))
  off[off == "X"] <- as.character(coding$X)
  off[off == "Y"] <- as.character(coding$Y)

  out <- data.frame(Offset = suppressWarnings(as.numeric(off)),
                    Size   = as.numeric(sizes$Size),
                    stringsAsFactors = FALSE)
  out[!is.na(out$Offset), , drop = FALSE]
}


#### . . . - .mmgc_require() - . . . ####
#### . . . - Attach a group of packages on demand, with an actionable message.

.mmgc_require <- function(PACKAGES, WHAT) {
  missing <- character(0)
  for (p in PACKAGES) {
    if (!suppressPackageStartupMessages(
          require(p, character.only = TRUE, quietly = TRUE))) {
      missing <- c(missing, p)
    }
  }
  if (length(missing)) {
    stop(WHAT, " needs these packages, which are not installed: ",
         paste(missing, collapse = ", "),
         "\n  Install them with: Rscript scripts/install_dependencies.R",
         call. = FALSE)
  }
  invisible(TRUE)
}


#### . . . - .mmgc_outfile() - . . . ####
#### . . . - Join an output directory and a file name with exactly one separator.
####         The original pasted OUT.DIR and the name directly in some places and
####         with a "/" in others, so a directory given with a trailing slash
####         produced paths with a doubled separator.

.mmgc_outfile <- function(OUT.DIR, ...) {
  file.path(sub("/+$", "", OUT.DIR), paste0(...))
}
