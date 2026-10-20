#### . . . - FocalLargeMatrix() - . . . ####

#### . . . - Score every recurrent GISTIC peak in every sample, distinguishing
####         focal from large events.
####
####         Gains  (peak lesion == "Amp"):        Losses (peak lesion == "Del"):
####           total >= 3.5, focal  -> 2             total  < 0.5, focal  -> 2  (biallelic)
####           total >= 3.5, large  -> 1.5           total  < 0.5, large  -> 1.5
####           2.5 <= total < 3.5, focal -> 1        monoallelic, focal   -> 1
####           2.5 <= total < 3.5, large -> 0.5      monoallelic, large   -> 0.5
####           otherwise -> 0                        otherwise -> 0
####
####         A segment counts as focal when it is at most FOCAL.LENGTH long.
####         Where several segments hit one peak, the peak takes the highest score.
####
####         CNV          = output of inputCNV()
####         GISTIC.FILE  = recurrent peak coordinates. Defaults to the active build's.
####         DATABASE     = cohort label, used in the output file name
####         FOCAL.LENGTH = focal/large cut-off in bp (default 5 Mb)
####         OUT.DIR      = where to write CNV_<DATABASE>matrix_focal_large_<n>Mb.txt
####         CYTOBAND.FILE, METADATA = accepted for backward compatibility; see below.
####
####         Returns the matrix invisibly as well as writing it, so it can be used
####         directly instead of being read back from disk.
####
####     NOTE on CYTOBAND.FILE and METADATA: the original signature required both,
####     but neither changes the result -- the cytoband table was merged into the
####     peaks only to eyeball centromeric peaks interactively, and the metadata
####     was read and never used. They stay in the signature so existing call
####     sites keep working, and are ignored.

FocalLargeMatrix <- function(CNV,
                             GISTIC.FILE = NULL,
                             DATABASE,
                             OUT.DIR,
                             FOCAL.LENGTH = 5000000,
                             CYTOBAND.FILE = NULL,
                             METADATA = NULL) {

  if (is.null(GISTIC.FILE)) GISTIC.FILE <- MMGC_ref("gistic_peaks")

  cnv <- as.data.frame(CNV, stringsAsFactors = FALSE)
  need <- c("sample", "chrom", "start", "end", "total", "minor")
  absent <- setdiff(need, colnames(cnv))
  if (length(absent)) {
    stop("FocalLargeMatrix(): missing column(s) in CNV: ",
         paste(absent, collapse = ", "), "\n  Run inputCNV() first.", call. = FALSE)
  }

  gistic <- read.delim(GISTIC.FILE, stringsAsFactors = FALSE)
  for (col in c("chr", "start", "end", "cytoband", "code", "lesion")) {
    if (!col %in% colnames(gistic)) {
      stop("FocalLargeMatrix(): ", GISTIC.FILE, " has no column ", sQuote(col),
           ".\n  It needs: chr, start, end, cytoband, code, lesion.", call. = FALSE)
    }
  }
  gistic$chr <- gsub("^chr", "", as.character(gistic$chr))
  gistic$gistic.seg.length <- gistic$end - gistic$start

  cnv$type       <- DATABASE
  cnv$seg.length <- cnv$end - cnv$start
  cnv$seg.type   <- ifelse(cnv$seg.length <= FOCAL.LENGTH, "focal", "large")
  cnv$chrom      <- gsub("^chr", "", as.character(cnv$chrom))

  #### . . . - Segment x peak overlaps
  hits <- .mmgc_overlap(
    A = cnv, B = gistic,
    A.CHROM = "chrom", A.START = "start", A.END = "end",
    B.CHROM = "chr",   B.START = "start", B.END = "end",
    A.KEEP = c("sample", "total", "minor", "type", "seg.length", "seg.type"),
    B.KEEP = c("cytoband", "code", "lesion", "gistic.seg.length"))

  if (!nrow(hits)) {
    stop("FocalLargeMatrix(): no copy-number segment overlaps any GISTIC peak. ",
         "Is CNV on the same build as ", sQuote(MMGC_genome()), "?", call. = FALSE)
  }
  hits <- hits %>% dplyr::arrange(sample, seqnames, startA)

  #### . . . - Gains
  amp <- hits[which(hits$lesion == "Amp"), , drop = FALSE]
  amp$anomaly <- ifelse(amp$total >= 3.5 & amp$seg.type == "focal", 2,
                 ifelse(amp$total >= 3.5 & amp$seg.type == "large", 1.5,
                 ifelse(amp$total >= 2.5 & amp$total < 3.5 & amp$seg.type == "focal", 1,
                 ifelse(amp$total >= 2.5 & amp$total < 3.5 & amp$seg.type == "large", 0.5, 0))))
  amp.mtx <- .mmgc_peak_matrix(amp, unique(cnv$sample), PREFIX = "Amp_")

  #### . . . - Losses
  del <- hits[which(hits$lesion == "Del"), , drop = FALSE]
  # Sentinel kept from the original: a segment at total >= 2 whose minor allele is
  # unknown cannot be called, and ends up NA in the matrix (see .mmgc_max_by_cell).
  del$pr <- ifelse(del$total >= 2 & is.na(del$minor), 0, "qqq")
  del$anomaly <- ifelse(del$total < 0.5 & del$seg.type == "focal", 2,
                 ifelse(del$total < 0.5 & del$seg.type == "large", 1.5,
                 ifelse(del$total >= 1.5 & del$minor < 0.5 & del$seg.type == "focal", 1,
                 ifelse(del$total < 1.5 & del$total >= 0.5 & del$seg.type == "focal", 1,
                 ifelse(del$total >= 1.5 & del$minor < 0.5 & del$seg.type == "large", 0.5,
                 ifelse(del$total < 1.5 & del$total >= 0.5 & del$seg.type == "large", 0.5, 0))))))
  del$anomaly <- ifelse(is.na(del$anomaly), del$pr, del$anomaly)
  del.mtx <- .mmgc_peak_matrix(del, unique(cnv$sample), PREFIX = "Del_")

  #### . . . - One row per sample, one column per peak
  all.peaks <- rbind(amp.mtx, del.mtx)
  out <- as.data.frame(t(all.peaks), stringsAsFactors = FALSE)
  for (j in seq_len(ncol(out))) out[[j]] <- as.numeric(as.character(out[[j]]))

  #### . . . - Peaks that nothing in this cohort touched still need a column.
  ####         The original had two attempts at this; the first looped over
  ####         `length(mis.peak)` instead of `seq_along(...)`, so it only ever
  ####         added the last missing peak.
  absent.peaks <- setdiff(gistic$code, colnames(out))
  for (p in absent.peaks) out[[p]] <- NA_real_

  out$database <- DATABASE

  focal.mb <- FOCAL.LENGTH / 1e6
  outfile <- .mmgc_outfile(
    OUT.DIR, "CNV_", DATABASE, "matrix_focal_large_",
    if (focal.mb == round(focal.mb)) sprintf("%d", round(focal.mb)) else format(focal.mb),
    "Mb.txt")

  write.table(out, outfile, col.names = TRUE, row.names = TRUE,
              quote = FALSE, sep = "\t")
  message("FocalLargeMatrix(): ", nrow(out), " samples x ",
          ncol(out) - 1L, " peaks -> ", outfile)

  invisible(out)
}


#### . . . - .mmgc_peak_matrix() - . . . ####
#### . . . - Peak x sample matrix of anomaly scores, taking the highest score when
####         several segments hit the same peak, and filling samples with no
####         qualifying event with 0.

.mmgc_peak_matrix <- function(HITS, ALL.SAMPLES, PREFIX = "") {

  if (!nrow(HITS)) {
    return(matrix(numeric(0), nrow = 0, ncol = length(ALL.SAMPLES),
                  dimnames = list(character(0), ALL.SAMPLES)))
  }

  dt <- data.table::as.data.table(unique(HITS[, c("sample", "cytoband", "anomaly")]))
  wide <- data.table::dcast(dt, cytoband ~ sample, toString, value.var = "anomaly")
  wide <- as.data.frame(wide, stringsAsFactors = FALSE)
  rownames(wide) <- wide$cytoband
  wide <- wide[, -1, drop = FALSE]

  wide <- .mmgc_max_by_cell(wide)

  #### . . . - Samples with no event in any peak of this class score 0, not NA:
  ####         their copy number was assessed, it simply was not abnormal.
  absent <- setdiff(ALL.SAMPLES, colnames(wide))
  for (s in absent) wide[[s]] <- 0

  rownames(wide) <- paste0(PREFIX, rownames(wide))
  wide
}
