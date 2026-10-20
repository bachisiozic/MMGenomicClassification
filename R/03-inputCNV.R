#### . . . - inputCNV() - . . . ####

#### . . . - Normalise an allele-specific copy-number table into the single shape
####         the rest of the pipeline expects:
####
####             sample  chrom  start  end  total  minor
####
####         where chrom is an integer: autosomes keep their number, X and Y take
####         the next two codes for the active build (human X = 23, Y = 24;
####         mouse X = 20, Y = 21). Rows on contigs outside that set -- scaffolds,
####         alts, chrM -- are dropped.
####
####         DATA       = data frame of copy-number segments
####         SAMPLE.COL, CHR.COL, START.COL, END.COL, TOTAL.COL, MINOR.COL
####                    = the column in DATA holding each field
####         GENOME     = build to use for chromosome coding. Defaults to the one
####                      set by MMGC_setup().
####
####         Returns the normalised data frame, sorted by sample, chrom, start.

inputCNV <- function(DATA,
                     SAMPLE.COL, CHR.COL, START.COL, END.COL, TOTAL.COL, MINOR.COL,
                     GENOME = NULL) {

  genome <- if (is.null(GENOME)) MMGC_genome() else GENOME
  coding <- MMGC_CHROM_CODING[[genome]]
  if (is.null(coding)) {
    stop("inputCNV(): unknown genome ", sQuote(genome), call. = FALSE)
  }

  DATA <- as.data.frame(DATA, stringsAsFactors = FALSE)

  wanted <- c(sample = SAMPLE.COL, chrom = CHR.COL, start = START.COL,
              end = END.COL, total = TOTAL.COL, minor = MINOR.COL)
  absent <- wanted[!wanted %in% colnames(DATA)]
  if (length(absent)) {
    stop("inputCNV(): column(s) not found in DATA: ",
         paste(sprintf("%s = %s", names(absent), sQuote(unname(absent))),
               collapse = ", "),
         "\n  DATA has: ", paste(colnames(DATA), collapse = ", "), call. = FALSE)
  }

  # Select and rename in one step. The original rewrote colnames with a stack of
  # nested ifelse() calls, which silently produced duplicate names when the input
  # already had a column called "start", "end", "total" or "minor".
  out <- DATA[, unname(wanted), drop = FALSE]
  colnames(out) <- names(wanted)

  out$start <- as.numeric(as.character(out$start))
  out$end   <- as.numeric(as.character(out$end))
  out$total <- as.numeric(as.character(out$total))
  out$minor <- as.numeric(as.character(out$minor))

  #### . . . - Chromosome to integer code
  chrom <- gsub("^chr", "", as.character(out$chrom), ignore.case = TRUE)
  chrom[chrom %in% c("X", "x")] <- as.character(coding$X)
  chrom[chrom %in% c("Y", "y")] <- as.character(coding$Y)

  keep.codes <- as.character(c(coding$autosomes, coding$X, coding$Y))
  dropped <- !chrom %in% keep.codes
  if (any(dropped)) {
    lost <- sort(unique(as.character(out$chrom)[dropped]))
    message("inputCNV(): dropped ", sum(dropped), " segment(s) on contigs outside ",
            genome, ": ", paste(utils::head(lost, 10), collapse = ", "),
            if (length(lost) > 10) sprintf(" (+%d more)", length(lost) - 10) else "")
  }
  out <- out[!dropped, , drop = FALSE]
  out$chrom <- as.numeric(chrom[!dropped])

  if (!nrow(out)) {
    stop("inputCNV(): no segments left after filtering to ", genome,
         " chromosomes. Check CHR.COL = ", sQuote(CHR.COL), ".", call. = FALSE)
  }

  out <- out %>% dplyr::arrange(sample, chrom, start)
  out[, c("sample", "chrom", "start", "end", "total", "minor")]
}
