#### . . . - Cytogenetic feature groups - . . . ####
#
# Which GISTIC peaks feed which downstream matrix. These were written inline as
# literal vectors in the original; naming them once means the human and mouse
# pipelines can share them, and that a cohort missing a peak fails with a clear
# message instead of a subscript error.
#
# mm10 keeps the human peak names on purpose (the peaks are relocated to the
# orthologous mouse locus but keep their human code), so these lists apply to the
# mouse build too -- minus the peaks with no confident mouse anchor, which
# MMGC_peak_groups() reports.

MMGC_PEAKS_AMP <- c("Amp_14q11.2", "Amp_14q32.33", "Amp_17q22", "Amp_2p11.2",
                    "Amp_2q24.3", "Amp_2q32.1", "Amp_4q13.2", "Amp_6p24.3",
                    "Amp_6q11.1", "Amp_8q24.21")

MMGC_PEAKS_DEL <- c("Del_10p15.3", "Del_10q24.32", "Del_10q26.3", "Del_11q22.1",
                    "Del_12p13.2", "Del_12q24.31", "Del_13q12.11", "Del_13q14.2",
                    "Del_14q24.3", "Del_14q32.33", "Del_16p11.2", "Del_16p13.3",
                    "Del_16q12.1", "Del_16q22.2", "Del_16q24.3", "Del_17p13.1",
                    "Del_1p12", "Del_1p22.1", "Del_1p32.3", "Del_20p13",
                    "Del_20q13.12", "Del_21p11.2", "Del_22q11.22", "Del_22q13.32",
                    "Del_2p11.2", "Del_2q31.1", "Del_2q37.3", "Del_6q26",
                    "Del_7p11.2", "Del_7p22.2", "Del_8p21.3", "Del_8p23.1",
                    "Del_8q24.21", "Del_9p21.3")

MMGC_PEAKS_1Q <- c("Amp_1q21.1", "Amp_1q21.3", "Amp_1q44")

#### . . . - GISTIC deletion peaks that stand in for a tumour suppressor, and the
####         gene they represent. Used to cross-check the TSG calls.
MMGC_PEAK_TSG <- c(Del_1p12    = "Del_TENT5C",
                   Del_17p13.1 = "Del_TP53",
                   Del_1p32.3  = "Del_CDKN2C",
                   Del_1p22.1  = "Del_RPL5")


#### . . . - CytogeneticMatrixGenerator() - . . . ####

#### . . . - Turn the focal/large peak matrix plus the hyperdiploidy call into the
####         three tables the final matrix is assembled from:
####
####           HRD_chromosomes<DATABASE>.txt    per-chromosome trisomy flags
####           gistic.other_All_peaks<DB>.txt   recurrent peaks, 0/1/2
####           x1q_amp_gain<DATABASE>.txt       1q gain/amplification, 0/1/2
####
####         Peak scores are collapsed from the focal/large scale to severity:
####         0.5 or 1 (gain / monoallelic loss) -> 1, and 1.5 or 2 (amplification /
####         biallelic loss) -> 2.
####
####         FocalLargeMATRIX = path written by FocalLargeMatrix(), or the data
####                            frame it returns
####         HY.FILE          = path written by HYidentification(), or the
####                            `hyperdiploid` element it returns
####         DATABASE         = cohort label
####         OUT.DIR          = output directory
####
####         Returns the three tables invisibly as a named list.

CytogeneticMatrixGenerator <- function(FocalLargeMATRIX, HY.FILE, DATABASE, OUT.DIR) {

  peaks <- if (is.character(FocalLargeMATRIX))
    read.delim(FocalLargeMATRIX, stringsAsFactors = FALSE) else
    as.data.frame(FocalLargeMATRIX, stringsAsFactors = FALSE)

  hy <- if (is.character(HY.FILE))
    read.delim(HY.FILE, stringsAsFactors = FALSE) else
    as.data.frame(HY.FILE, stringsAsFactors = FALSE)

  if (!"sample" %in% colnames(hy)) {
    stop("CytogeneticMatrixGenerator(): HY.FILE has no 'sample' column. ",
         "Pass the *hyperdiploid* table from HYidentification(), not the ",
         "no.hyperdiploid one.", call. = FALSE)
  }

  ################################################################################
  #### . . . - HRD chromosomes - . . . ####
  ################################################################################

  hyp.cols <- grep("^hyp\\.chr", colnames(hy), value = TRUE)
  if (!length(hyp.cols)) {
    stop("CytogeneticMatrixGenerator(): HY.FILE carries no hyp.chr* columns.",
         call. = FALSE)
  }

  hrd <- data.frame(row.names = rownames(peaks))
  for (cc in hyp.cols) {
    v <- hy[[cc]][match(rownames(peaks), hy$sample)]
    #### . . . - A sample absent from the hyperdiploid table has no gained
    ####         chromosome, which is a 0, not missing data.
    v[is.na(v)] <- 0
    #### . . . - Column naming follows the hyp.chr* name rather than a fixed
    ####         positional list, so a different HRD chromosome set still lines up.
    hrd[[paste0("chr", sub("^hyp\\.chr", "", cc), "+")]] <- as.numeric(v)
  }

  f.hrd <- .mmgc_outfile(OUT.DIR, "HRD_chromosomes", DATABASE, ".txt")
  write.table(hrd, f.hrd, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

  ################################################################################
  #### . . . - Other recurrent GISTIC peaks - . . . ####
  ################################################################################

  amp <- .mmgc_select_peaks(peaks, MMGC_PEAKS_AMP, "amplification")
  del <- .mmgc_select_peaks(peaks, MMGC_PEAKS_DEL, "deletion")

  other <- cbind(.mmgc_severity(amp), .mmgc_severity(del))
  other$database <- DATABASE

  f.other <- .mmgc_outfile(OUT.DIR, "gistic.other_All_peaks", DATABASE, ".txt")
  write.table(other, f.other, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

  ################################################################################
  #### . . . - 1q - . . . ####
  ################################################################################

  x1q <- .mmgc_severity(.mmgc_select_peaks(peaks, MMGC_PEAKS_1Q, "1q"))
  x1q$database <- DATABASE

  f.1q <- .mmgc_outfile(OUT.DIR, "x1q_amp_gain", DATABASE, ".txt")
  write.table(x1q, f.1q, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

  message("CytogeneticMatrixGenerator(): ", nrow(hrd), " samples -> ",
          basename(f.hrd), ", ", basename(f.other), ", ", basename(f.1q))

  invisible(list(hrd_chromosomes = hrd, gistic_other = other, x1q = x1q))
}


#### . . . - .mmgc_select_peaks() - . . . ####
#### . . . - Pull a named set of peak columns, adding NA columns for any the
####         cohort's reference did not cover. The original indexed by name
####         directly, so one missing peak aborted the run with "undefined
####         columns selected" and no indication which.

.mmgc_select_peaks <- function(X, WANTED, WHAT) {
  absent <- setdiff(WANTED, colnames(X))
  if (length(absent)) {
    message("  note: ", length(absent), " ", WHAT, " peak(s) absent from this ",
            "cohort's matrix, carried as NA: ", paste(absent, collapse = ", "))
  }
  out <- data.frame(row.names = rownames(X))
  for (p in WANTED) {
    out[[p]] <- if (p %in% colnames(X)) X[[p]] else NA_real_
  }
  out
}


#### . . . - .mmgc_severity() - . . . ####
#### . . . - Collapse focal/large scores to severity: {0.5, 1} -> 1 (single copy
####         gain or loss), {1.5, 2} -> 2 (amplification or biallelic loss),
####         0 -> 0, NA -> NA.

.mmgc_severity <- function(X) {
  for (j in seq_len(ncol(X))) {
    v <- as.numeric(X[[j]])
    X[[j]] <- ifelse(is.na(v), NA_real_,
              ifelse(v %in% c(0.5, 1), 1,
              ifelse(v %in% c(1.5, 2), 2,
              ifelse(v == 0, 0, v))))
  }
  X
}
