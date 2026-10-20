#### . . . - The published feature set - . . . ####
#
# The columns of the final matrix, in the order MMGenoCla() expects. Order is not
# cosmetic: MMGenoCla() derives its "deletions" burden by taking the span from the
# first loss column to the last tumour-suppressor column, so copy-number losses
# and TSGs have to stay contiguous and in this sequence.
#
# Two peaks that FocalLargeMatrix() scores are deliberately left out here, as in
# the published code: Amp_14q11.2 and Del_21p11.2.

MMGC_FEATURES_HRD_CHROM <- paste0("CNV_chr", c(3, 5, 7, 9, 11, 15, 19, 21, 18), ".gain")

MMGC_FEATURES_GAIN <- c("CNV_Gain_Amp1q", "CNV_Amp_17q22", "CNV_Amp_2q24.3",
                        "CNV_Amp_2q32.1", "CNV_Amp_4q13.2", "CNV_Amp_6p24.3",
                        "CNV_Amp_8q24.21")

MMGC_FEATURES_LOSS <- c("CNV_chr8p.loss", "CNV_Del_10p15.3", "CNV_Del_10q24.32",
                        "CNV_Del_10q26.3", "CNV_Del_11q22.1", "CNV_Del_12p13.2",
                        "CNV_Del_12q24.31", "CNV_Del_14q24.3", "CNV_Del_20p13",
                        "CNV_Del_20q13.12", "CNV_Del_22q13.32", "CNV_Del_2q31.1",
                        "CNV_Del_2q37.3", "CNV_Del_6q26", "CNV_Del_7p11.2",
                        "CNV_Del_7p22.2", "CNV_Del_8q24.21", "CNV_Del_9p21.3")

MMGC_FEATURES_TSG <- paste0("CNV.SNV_", c(
  "ARID1A", "ARID2", "ATM", "BTG1", "CDKN1B", "CDKN2C", "CREBBP", "CYLD",
  "DNMT3A", "EP300", "TENT5C", "FUBP1", "KMT2B", "KMT2C", "MAX", "NCOR1",
  "NF1", "NFKBIA", "POT1", "PRDM1", "RASA2", "RB1", "RPL5", "SETD2", "SP140",
  "TET2", "TGDS", "TP53", "TRAF2", "TRAF3", "UBR5", "ZNF292"))

MMGC_FEATURES_SNV <- paste0("SNV_", c(
  "BRAF", "CCND1", "DIS3", "FGFR3", "H1.4", "IDH1", "IDH2", "IRF4", "KRAS",
  "NFKB2", "NRAS", "PIK3CA", "PTPN11", "SF3B1", "ABCF1", "ACTG1", "BHLHE41",
  "DDX3X", "DTX1", "DUSP2", "H1.5", "H1.2", "H1.3", "H2BC12", "HLA.C", "HUWE1",
  "INO80E", "IRF1", "KLHL6", "LCE1D", "LTB", "MAN2C1", "MRGPRX4", "NDUFC2",
  "PABPC1", "PANK3", "PRKD2", "RFTN1", "RPRD1B", "RPS3A", "SAMHD1", "SLC35G5",
  "SYT17", "TBC1D29P", "TNFRSF11B", "USP7", "ZFP36L1", "ATRX", "KDM5C", "KDM6A",
  "RPL10", "BCL7A", "EGR1", "MAF", "MAFB", "MAML2", "PIM1", "TCL1A", "XBP1"))

MMGC_FINAL_FEATURES <- c("sample",
                         MMGC_FEATURES_HRD_CHROM,
                         MMGC_FEATURES_GAIN,
                         MMGC_FEATURES_LOSS,
                         MMGC_FEATURES_TSG,
                         MMGC_FEATURES_SNV,
                         "APOBEC")

#### . . . - Peaks that are collapsed into one feature, or dropped because a
####         gene-level call already covers them.
MMGC_1Q_PEAKS  <- c("CNV_Amp_1q21.1", "CNV_Amp_1q21.3", "CNV_Amp_1q44")
MMGC_8P_PEAKS  <- c("CNV_Del_8p21.3", "CNV_Del_8p23.1")

#### . . . - Peaks whose signal is already represented by a tumour-suppressor
####         feature, or which are immunoglobulin loci and not interpretable.
MMGC_DROP_PEAKS <- c("Del_1p32.3", "Del_1p12", "Del_1p22.1", "Del_13q14.2",
                     "Del_16p13.3", "Del_16q12.1", "Del_17p13.1",
                     "Amp_14q32.33", "Del_14q32.33", "Del_22q11.22",
                     "Del_2p11.2", "Amp_2p11.2")


#### . . . - FinalMatrixGenerator() - . . . ####

#### . . . - Assemble the per-sample feature matrix from the tables the earlier
####         steps wrote, and reduce it to the published feature set.
####
####         INPUT.FILE.FOLDER = folder holding the outputs of
####                             CytogeneticMatrixGenerator(), MutationsMatrixGenerator()
####                             and APOBECanalyzer()
####         DATABASE          = cohort label
####         OUTPUT            = path for the final matrix
####         FEATURES          = feature set to emit (default MMGC_FINAL_FEATURES)
####         STRICT            = stop when a feature is missing, instead of
####                             carrying it as NA
####
####         Returns the matrix invisibly.

FinalMatrixGenerator <- function(INPUT.FILE.FOLDER,
                                 DATABASE,
                                 OUTPUT,
                                 FEATURES = MMGC_FINAL_FEATURES,
                                 STRICT = FALSE) {

  rd <- function(...) {
    f <- .mmgc_outfile(INPUT.FILE.FOLDER, ...)
    if (!file.exists(f)) {
      stop("FinalMatrixGenerator(): expected file not found:\n  ", f,
           "\n  Run the earlier pipeline steps into this folder first.",
           call. = FALSE)
    }
    read.delim(f, stringsAsFactors = FALSE)
  }

  ################################################################################
  #### . . . - Copy number - . . . ####
  ################################################################################

  #### . . . - HRD chromosomes. The column names are written as "chr3+", which
  ####         read.delim() turns into "chr3." -- the original relied on exactly
  ####         that and pasted "gain" onto it. Strip whichever trailing character
  ####         is there, so the result is "CNV_chr3.gain" either way.
  hrd <- rd("HRD_chromosomes", DATABASE, ".txt")
  colnames(hrd) <- paste0("CNV_", sub("[.+]$", "", colnames(hrd)), ".gain")

  #### . . . - 1q
  x1q <- rd("x1q_amp_gain", DATABASE, ".txt")
  x1q <- x1q[, colnames(x1q) != "database", drop = FALSE]
  colnames(x1q) <- paste0("CNV_", colnames(x1q))

  cn <- merge(hrd, x1q, by = "row.names", all = TRUE)
  colnames(cn)[1] <- "sample"

  #### . . . - Other recurrent peaks, reduced to present/absent
  gistic <- rd("gistic.other_All_peaks", DATABASE, ".txt")
  gistic <- gistic[, colnames(gistic) != "database", drop = FALSE]
  gistic <- gistic[, !colnames(gistic) %in% MMGC_DROP_PEAKS, drop = FALSE]
  for (j in seq_len(ncol(gistic))) {
    v <- as.numeric(gistic[[j]])
    gistic[[j]] <- ifelse(is.na(v), NA_real_, ifelse(v > 1, 1, v))
  }
  colnames(gistic) <- paste0("CNV_", colnames(gistic))

  m <- merge(cn, gistic, by.x = "sample", by.y = "row.names", all = TRUE)

  ################################################################################
  #### . . . - Mutations - . . . ####
  ################################################################################

  onc <- rd("oncogene_", DATABASE, ".txt")
  for (j in seq_len(ncol(onc))) {
    v <- as.numeric(onc[[j]])
    onc[[j]] <- ifelse(is.na(v), NA_real_, ifelse(v > 1, 1, v))
  }
  colnames(onc) <- paste0("SNV_", colnames(onc))
  m <- merge(m, onc, by.x = "sample", by.y = "row.names", all = TRUE)

  tsg <- rd("tsg_", DATABASE, ".txt")
  colnames(tsg) <- paste0("CNV.SNV_", colnames(tsg))
  m <- merge(m, tsg, by.x = "sample", by.y = "row.names", all = TRUE)

  unk <- rd("unk_mutated_genes_", DATABASE, ".txt")
  colnames(unk) <- paste0("SNV_", colnames(unk))
  m <- merge(m, unk, by.x = "sample", by.y = "row.names", all = TRUE)

  ################################################################################
  #### . . . - APOBEC - . . . ####
  ################################################################################

  apo <- rd("APOBEC_relative_contribution_matrix", DATABASE, ".txt")
  apo <- apo[, "apobec", drop = FALSE]
  m <- merge(m, apo, by.x = "sample", by.y = "row.names", all = TRUE)

  #### . . . - Binned on the published cohort's deciles, not this cohort's.
  m$APOBEC <- ifelse(is.na(m$apobec), NA_real_,
              ifelse(m$apobec <  MMGC_APOBEC_CUTOFFS[["low"]],  0,
              ifelse(m$apobec <  MMGC_APOBEC_CUTOFFS[["high"]], 1, 2)))
  m <- m[, colnames(m) != "apobec", drop = FALSE]

  ################################################################################
  #### . . . - X-linked drivers: mutation only - . . . ####
  ################################################################################

  #### . . . - Drop the combined (deletion + mutation) call for these and use the
  ####         mutation-only one, because a single-copy loss on X is not a hit.
  m <- m[, !colnames(m) %in% paste0("CNV.SNV_", MMGC_CHRX_GENES), drop = FALSE]
  chrx <- rd("chrx_tsg_", DATABASE, ".txt")
  if (ncol(chrx)) {
    colnames(chrx) <- paste0("SNV_", colnames(chrx))
    #### . . . - These columns can already be present from the oncogene/unknown
    ####         matrices; keep the mutation-only version.
    m <- m[, !colnames(m) %in% colnames(chrx), drop = FALSE]
    m <- merge(m, chrx, by.x = "sample", by.y = "row.names", all.x = TRUE)
  }

  ################################################################################
  #### . . . - Collapse the multi-peak features - . . . ####
  ################################################################################

  #### . . . - 1q: gain in any of the three peaks -> 1, amplification in any -> 2.
  m$CNV_Gain_Amp1q <- .mmgc_collapse_peaks(m, MMGC_1Q_PEAKS)
  #### . . . - 8p: either peak lost -> 1. The published code never assigned 2 here.
  m$CNV_chr8p.loss <- .mmgc_collapse_peaks(m, MMGC_8P_PEAKS, MAX = 1)

  #### . . . - IGLL5 sits inside the IGL locus; its "mutations" are somatic
  ####         hypermutation, not driver events.
  ####         The original used df[, -grep("IGLL5", colnames(df))], which returns
  ####         a ZERO-column table when IGLL5 is absent -- so a cohort that had
  ####         never called it lost the entire matrix here.
  m <- .mmgc_drop_cols(m, "IGLL5")

  ################################################################################
  #### . . . - Reduce to the published feature set - . . . ####
  ################################################################################

  absent <- setdiff(FEATURES, colnames(m))
  if (length(absent)) {
    msg <- paste0("FinalMatrixGenerator(): ", length(absent),
                  " feature(s) not produced by this cohort:\n  ",
                  paste(absent, collapse = ", "))
    if (isTRUE(STRICT)) stop(msg, call. = FALSE)
    message(msg, "\n  Carried as NA. MMGenoCla() will report them again and, ",
            "because it drops samples with any NA, they must be resolved ",
            "before classification.")
    for (f in absent) m[[f]] <- NA_real_
  }

  out <- m[, FEATURES, drop = FALSE]

  write.table(out, OUTPUT, col.names = TRUE, row.names = FALSE,
              quote = FALSE, sep = "\t")
  message("FinalMatrixGenerator(): ", nrow(out), " samples x ",
          ncol(out) - 1L, " features -> ", OUTPUT)

  invisible(out)
}


#### . . . - .mmgc_collapse_peaks() - . . . ####
#### . . . - Combine several peak columns covering one region into a single
####         severity score: the highest score across the peaks, capped at MAX.
####         NA only when every contributing peak is NA.
####
####         The original spelled this out as a lookup over every observed
####         combination of three values pasted together ("1_2_0", "NA_1_1", ...),
####         which meant any combination nobody had seen yet fell through to NA.

.mmgc_collapse_peaks <- function(X, PEAKS, MAX = 2) {
  have <- intersect(PEAKS, colnames(X))
  if (!length(have)) return(rep(NA_real_, nrow(X)))

  vals <- vapply(have, function(p) as.numeric(X[[p]]), numeric(nrow(X)))
  vals <- matrix(vals, nrow = nrow(X))

  all.missing <- rowSums(!is.na(vals)) == 0
  out <- suppressWarnings(apply(vals, 1, function(v) max(v, na.rm = TRUE)))
  out[!is.finite(out)] <- NA_real_
  out[out > MAX] <- MAX
  out[all.missing] <- NA_real_
  out
}
