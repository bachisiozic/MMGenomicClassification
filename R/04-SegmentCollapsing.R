#### . . . - SegmentCollapsing() - . . . ####

#### . . . - Merge runs of consecutive segments that carry the same rounded total
####         and minor copy number, per sample and chromosome.
####
####         COPY.NUMBER = output of inputCNV()
####         DATABASE    = cohort label, written to the `type` column
####
####         Returns  sample chrom start end total minor code Num_Segments type,
####         where `code` is "<total>_<minor>" and Num_Segments is how many input
####         segments were merged into that row.
####
####         The original walked every sample and chromosome with rle() on the
####         code, kept the first and last row of each run, re-ran rle() on that
####         reduced table, and paired starts with ends positionally. That relies
####         on the two halves lining up row for row, which only holds because the
####         input is sorted -- and it silently produced mismatched rows when it
####         was not. This does the same collapse by grouping on a run id, which
####         cannot get out of step.
####
####         COPY.NUMBER is the first argument; `copy_number` is still accepted as
####         a name so existing scripts keep working.

SegmentCollapsing <- function(COPY.NUMBER, DATABASE, copy_number = NULL) {

  if (!is.null(copy_number)) COPY.NUMBER <- copy_number
  cnv <- as.data.frame(COPY.NUMBER, stringsAsFactors = FALSE)

  need <- c("sample", "chrom", "start", "end", "total", "minor")
  absent <- setdiff(need, colnames(cnv))
  if (length(absent)) {
    stop("SegmentCollapsing(): missing column(s): ", paste(absent, collapse = ", "),
         "\n  Run the table through inputCNV() first.", call. = FALSE)
  }

  cnv$total <- round(as.numeric(cnv$total), 0)
  cnv$minor <- round(as.numeric(cnv$minor), 0)
  cnv$code  <- paste(cnv$total, cnv$minor, sep = "_")

  cnv <- cnv %>% dplyr::arrange(sample, chrom, start)

  # A run is a maximal stretch of identical code within one sample+chromosome.
  key <- paste(cnv$sample, cnv$chrom, cnv$code, sep = "\r")
  run <- cumsum(c(TRUE, key[-1L] != key[-length(key)]))

  out <- cnv %>%
    dplyr::mutate(.run = run) %>%
    dplyr::group_by(.run) %>%
    dplyr::summarise(
      sample       = dplyr::first(sample),
      chrom        = dplyr::first(chrom),
      start        = dplyr::first(start),
      end          = dplyr::last(end),
      total        = dplyr::first(total),
      minor        = dplyr::first(minor),
      code         = dplyr::first(code),
      Num_Segments = dplyr::n(),
      .groups      = "drop"
    ) %>%
    dplyr::arrange(.run) %>%
    dplyr::select(-.run) %>%
    as.data.frame()

  out$type <- DATABASE
  out
}
