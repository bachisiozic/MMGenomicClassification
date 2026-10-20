#### . . . - The classification itself - . . . ####
#
# Eleven genomic groups, assigned by an ordered set of rules. The order matters:
# later rules overwrite earlier ones, which is how a sample that satisfies several
# definitions ends up in the most specific group.
#
# The published function carried this whole body TWICE -- once for the branch where
# every expected feature was present and once for the branch where some were
# missing -- 245 lines duplicated, differing only in that the second branch also
# returned the dropped features. Here it lives once.


#### . . . - Features carried in the matrix but excluded from the classification,
####         because they were called in too few of the 1,933 samples the
####         classification was derived from to be usable.
MMGC_HIGH_MISSING_FEATURES <- paste0("SNV_", c(
  "IDH1", "IDH2", "PIK3CA", "ABCF1", "BHLHE41", "DDX3X", "DTX1", "H1.5", "H1.3",
  "H2BC12", "INO80E", "IRF1", "LCE1D", "MAN2C1", "MRGPRX4", "NDUFC2", "PANK3",
  "RFTN1", "RPRD1B", "RPS3A", "SLC35G5", "SYT17", "TBC1D29P", "TNFRSF11B",
  "USP7", "ZFP36L1", "ATRX", "KDM5C", "KDM6A", "RPL10", "MAF", "MAFB", "MAML2",
  "PIM1", "TCL1A"))

#### . . . - Chromosomes whose trisomy counts towards hyperdiploidy. At least two
####         are needed. chr18 is carried as a feature but excluded here.
MMGC_HRD_GAIN_FEATURES <- paste0("CNV_chr", c(3, 5, 7, 9, 11, 15, 19, 21), ".gain")

#### . . . - Additional gains. Three or more of these marks a "gains" phenotype.
MMGC_ADD_GAIN_FEATURES <- c("CNV_Amp_14q11.2", "CNV_Amp_17q22", "CNV_Amp_2q24.3",
                            "CNV_Amp_2q32.1", "CNV_Amp_4q13.2", "CNV_Amp_6p24.3",
                            "CNV_Amp_8q24.21", "CNV_chr18.gain")

#### . . . - MAPK-pathway and FGFR3 mutations, which define the HRD_RAS group.
MMGC_RAS_FEATURES <- c("SNV_KRAS", "SNV_NRAS", "SNV_BRAF", "SNV_FGFR3")

#### . . . - Row order for the heatmap: features grouped by chromosome, then
####         mutations, then translocations and the composite features.
MMGC_HEATMAP_ROWS <- c(
  "chr3.gain", "chr5.gain", "chr7.gain", "chr9.gain", "chr11.gain",
  "chr15.gain", "chr19.gain", "chr21.gain",                                # HRD
  "Gain_Amp1q",                                                           # 1q
  "chr18.gain",
  "Amp_17q22", "Amp_2q24.3", "Amp_2q32.1", "Amp_4q13.2", "Amp_6p24.3",
  "Amp_8q24.21",                                                          # other gains
  "ARID1A", "CDKN2C", "FAM46C", "FUBP1", "RPL5",                          # chr1
  "Del_2q31.1", "Del_2q37.3", "DNMT3A", "SP140",                          # chr2
  "RASA2", "SETD2",                                                       # chr3
  "TET2",                                                                 # chr4
  "Del_6q26", "PRDM1", "ZNF292",                                          # chr6
  "Del_7p11.2", "Del_7p22.2", "KMT2C", "POT1",                            # chr7
  "chr8p.loss", "Del_8q24.21", "UBR5",                                    # chr8
  "Del_9p21.3", "TRAF2",                                                  # chr9
  "Del_10p15.3", "Del_10q24.32", "Del_10q26.3",                           # chr10
  "Del_11q22.1", "ATM",                                                   # chr11
  "Del_12p13.2", "Del_12q24.31", "ARID2", "BTG1", "CDKN1B",               # chr12
  "RB1", "TGDS",                                                          # chr13
  "Del_14q24.3", "MAX", "NFKBIA", "TRAF3",                                # chr14
  "CREBBP", "CYLD",                                                       # chr16
  "NCOR1", "NF1", "TP53",                                                 # chr17
  "KMT2B",                                                                # chr19
  "Del_20p13", "Del_20q13.12",                                            # chr20
  "Del_22q13.32", "EP300",                                                # chr22
  "BRAF", "KRAS", "NRAS",                                                 # MAPK
  "CCND1", "FGFR3", "H1.4", "IRF4", "NFKB2", "PTPN11", "SF3B1", "ACTG1",
  "DUSP2", "H1.2", "HLA.C", "HUWE1", "KLHL6", "LTB", "PABPC1", "PRKD2",
  "SAMHD1", "BCL7A", "EGR1", "XBP1",                                      # other mutations
  "HS_DIS3", "noHS_DIS3",                                                 # DIS3
  "CCND1-IGH", "NSD2-IGH", "MAF-IGH",                                     # translocations
  "APOBEC",
  "CNV.Sig")

#### . . . - The eleven published groups, in the order they are reported, with the
####         internal rule labels that map onto each.
MMGC_GROUP_MAP <- list(
  HRD_RAS                   = "A_HRD_RAS",
  HRD_Gains                 = "B_HRD_gains",
  HRD_Complex_Cytogenetic   = "C_HRD_deletions",
  CCND1_Simple              = c("E_CCND1_HRD", "F_CCND1_RAS", "FF_CCND1_RAS"),
  CCND1_Complex_Cytogenetic = "G_CCND1_complex",
  MAF_and_or_HyperAPOBEC    = c("I_MAF_APOBEC", "L_MAF_APOBEC", "M_MAF_APOBEC"),
  NSD2_GainAmp1q_Del13q     = "N_NSD2_1q_del13q",
  NSD2_HRD                  = "O_NSD2",
  NSD2_Del13q               = "P_NSD2_del13q",
  GainAmp1q_Del13q          = "Q_x1q_del13q",
  Multiple_Losses           = "R_Multiple_losses",
  Simple                    = "RR_Simple")

MMGC_GROUP_COLOURS <- c(
  HRD_RAS                   = "lightpink3",
  HRD_Gains                 = "mediumpurple2",
  HRD_Complex_Cytogenetic   = "darkorchid4",
  CCND1_Simple              = "skyblue3",
  CCND1_Complex_Cytogenetic = "navy",
  MAF_and_or_HyperAPOBEC    = "darkgreen",
  NSD2_GainAmp1q_Del13q     = "lightgoldenrod",
  NSD2_HRD                  = "gold",
  NSD2_Del13q               = "goldenrod4",
  GainAmp1q_Del13q          = "sienna3",
  Multiple_Losses           = "firebrick",
  Simple                    = "slategray")


#### . . . - MMGenoCla() - . . . ####

#### . . . - Assign each sample to one of the eleven genomic groups.
####
####         input   = final matrix (data frame, or a path to it), with the
####                   translocation calls t_CCND1 / t_MMSET / t_MAF, the
####                   chromothripsis call (CNV.Sig or Chromothripsis) and the DIS3
####                   split from dis3_hs() merged in.
####         CONFIRM = ask before classifying, as the published function did.
####                   Defaults to FALSE so the pipeline can run unattended; the
####                   original blocked on two readline() prompts, which is why
####                   scripts calling it had stray "y" lines afterwards.
####         STRICT  = stop when expected features are missing, instead of warning.
####         EXPECTED = feature set to check for, defaulting to the published one
####                   (MMGC_FINAL_FEATURES). The mouse pipeline passes its reduced
####                   set so that features with no mouse ortholog are not reported.
####         QUIET   = suppress the per-group count printed on completion.
####         LEGACY.ADD.GAINS -- see the warning below.
####
####         Returns a list:
####           $genomic     sample x feature matrix, ready for pheatmap
####           $annotation  per-sample group assignment
####           $colors      colour list for pheatmap's annotation_colors
####           $clusters    the internal rule label per sample
####           $high.missing.features  the features excluded before classifying
####           $missing.features       expected features this matrix did not have
####
####     WARNING -- LEGACY.ADD.GAINS reproduces a bug in the published code. The
####     "additional gains" burden was summed over eight features including
####     CNV_Amp_14q11.2 -- but FinalMatrixGenerator() deliberately drops that peak,
####     so the column is absent, the sum evaluates to a zero-length vector, and the
####     subset assignment that should set add_gains silently does nothing.
####     add_gains therefore stays 0 for every sample, and the HRD_Gains group can
####     never be assigned. This version sums over whichever of the eight features
####     are present (keeping the >= 3 threshold) and reports which it used. Set
####     LEGACY.ADD.GAINS = TRUE to reproduce the published behaviour, in which
####     HRD_Gains is always empty. See docs/KNOWN_ISSUES.md.

MMGenoCla <- function(input,
                      CONFIRM = FALSE,
                      STRICT = FALSE,
                      LEGACY.ADD.GAINS = FALSE,
                      EXPECTED = NULL,
                      QUIET = FALSE) {

  ################################################################################
  #### . . . - Input - . . . ####
  ################################################################################

  all <- if (is.character(input))
    read.delim(input, stringsAsFactors = FALSE, header = TRUE) else
    as.data.frame(input, stringsAsFactors = FALSE)

  if (isTRUE(CONFIRM)) {
    cat("Expected column annotations:", sep = "\n")
    cat("  CNV:        CNV_chr3.gain, CNV_Gain_Amp1q, CNV_Amp_8q24.21, CNV_chr8p.loss, ...", sep = "\n")
    cat("  TSG:        CNV.SNV_TENT5C, CNV.SNV_RPL5, CNV.SNV_RB1, CNV.SNV_TP53, ...", sep = "\n")
    cat("  Mutations:  SNV_NRAS, SNV_BRAF, SNV_KRAS, SNV_CCND1, ...", sep = "\n")
    cat("  Other:      CNV.Sig, APOBEC, t_CCND1, t_MMSET, t_MAF", sep = "\n")
    if (!tolower(substr(trimws(readline("Proceed? [y/N] ")), 1, 1)) %in% "y") {
      message("MMGenoCla(): stopped at the user's request.")
      return(invisible(NULL))
    }
  }

  ################################################################################
  #### . . . - Harmonise naming - . . . ####
  ################################################################################

  #### . . . - Translocations may arrive under cytogenetic or gene names.
  cn <- colnames(all)
  cn[cn == "t(11;14)"]                 <- "t_CCND1"
  cn[cn %in% c("t(4;14)", "t_NSD2")]   <- "t_MMSET"
  cn[cn %in% c("t(14;16)", "t(14;20)")] <- "t_MAF"
  #### . . . - The chromothripsis call is used under the name CNV.Sig.
  cn[cn == "Chromothripsis"]           <- "CNV.Sig"
  #### . . . - The classification was defined on the hg19-era symbol FAM46C; the
  ####         hg38 reference files call the same gene TENT5C.
  cn[cn == "CNV.SNV_TENT5C"]           <- "CNV.SNV_FAM46C"
  colnames(all) <- cn

  #### . . . - Translocations are scored on the same 0/2 scale as an amplification.
  for (tr in c("t_CCND1", "t_MMSET", "t_MAF")) {
    if (tr %in% colnames(all)) {
      all[[tr]] <- ifelse(all[[tr]] == 1, 2, all[[tr]])
    }
  }

  ################################################################################
  #### . . . - Keep only the genomic features - . . . ####
  ################################################################################

  if (!"sample" %in% colnames(all)) {
    stop("MMGenoCla(): input has no 'sample' column.", call. = FALSE)
  }

  copynumbers  <- grep("CNV_", colnames(all), value = TRUE)
  tumoursupp   <- grep("CNV.SNV_", colnames(all), value = TRUE)
  mutations    <- c(grep("^SNV_", colnames(all), value = TRUE),
                    grep("^SNV\\.", colnames(all), value = TRUE))
  other        <- c("CNV.Sig", "APOBEC", "t_CCND1", "t_MMSET", "t_MAF")

  absent.other <- setdiff(other, colnames(all))
  if (length(absent.other)) {
    stop("MMGenoCla(): input is missing required column(s): ",
         paste(absent.other, collapse = ", "),
         "\n  t_CCND1 / t_MMSET / t_MAF come from your translocation calls and ",
         "CNV.Sig from the chromothripsis call -- merge them into the final ",
         "matrix before classifying. APOBEC comes from APOBECanalyzer().",
         call. = FALSE)
  }

  all <- all[, c("sample", copynumbers, tumoursupp, mutations, other), drop = FALSE]
  rownames(all) <- all$sample
  all <- all[, colnames(all) != "sample", drop = FALSE]

  ################################################################################
  #### . . . - Report expected features that are absent - . . . ####
  ################################################################################

  #### . . . - EXPECTED lets a caller narrow the check to the features its build can
  ####         actually measure -- the mouse pipeline passes its reduced set so that
  ####         features with no mouse ortholog are not reported as a problem.
  expected <- if (is.null(EXPECTED)) MMGC_FINAL_FEATURES else EXPECTED
  expected <- c(setdiff(expected, "sample"), "SNV.HS_DIS3", "SNV.noHS_DIS3", other)

  #### . . . - The matrix carries FAM46C after the rename above; the feature lists
  ####         use the current symbol TENT5C.
  expected <- unique(sub("^CNV\\.SNV_TENT5C$", "CNV.SNV_FAM46C", expected))

  missing.features <- setdiff(expected, colnames(all))

  if (length(missing.features)) {
    msg <- paste0("MMGenoCla(): ", length(missing.features),
                  " expected feature(s) absent from the matrix:\n  ",
                  paste(missing.features, collapse = ", "))
    if (isTRUE(STRICT)) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
    if (isTRUE(CONFIRM) &&
        !tolower(substr(trimws(readline("Classify anyway? [y/N] ")), 1, 1)) %in% "y") {
      message("MMGenoCla(): stopped at the user's request.")
      return(invisible(NULL))
    }
  }

  ################################################################################
  #### . . . - Drop the under-called features, then require complete data - . . . ####
  ################################################################################

  high.missing <- all[, colnames(all) %in% MMGC_HIGH_MISSING_FEATURES, drop = FALSE]
  all <- all[, !colnames(all) %in% MMGC_HIGH_MISSING_FEATURES, drop = FALSE]

  n.before <- nrow(all)
  keep <- stats::complete.cases(all)
  if (!any(keep)) {
    #### . . . - Name the offending features rather than just refusing: with ~100
    ####         features, "some column is NA" is not actionable.
    na.per.feature <- sort(colSums(is.na(all)), decreasing = TRUE)
    na.per.feature <- na.per.feature[na.per.feature > 0]
    stop("MMGenoCla(): every one of the ", n.before, " samples has at least one ",
         "missing feature, so none can be classified.",
         "\n  Features with missing values (worst first): ",
         paste(sprintf("%s (%d)", names(utils::head(na.per.feature, 20)),
                       utils::head(na.per.feature, 20)), collapse = ", "),
         "\n  Check that the translocation, chromothripsis and APOBEC calls cover ",
         "every sample.", call. = FALSE)
  }
  all <- all[keep, , drop = FALSE]
  if (nrow(all) < n.before) {
    message("MMGenoCla(): ", n.before - nrow(all), " of ", n.before,
            " samples dropped for incomplete features; ", nrow(all), " classified.")
  }

  ################################################################################
  #### . . . - Composite features - . . . ####
  ################################################################################

  gen <- all
  gen$sample <- rownames(gen)
  gen <- gen[, c(ncol(gen), seq_len(ncol(gen) - 1L)), drop = FALSE]

  #### . . . - Hyperdiploidy: two or more of the trisomy chromosomes.
  gen$HRD <- as.integer(.mmgc_feature_sum(gen, MMGC_HRD_GAIN_FEATURES) >= 2)

  #### . . . - Additional gains: three or more of the secondary gain peaks.
  add.present <- intersect(MMGC_ADD_GAIN_FEATURES, colnames(gen))
  if (isTRUE(LEGACY.ADD.GAINS)) {
    warning("MMGenoCla(): LEGACY.ADD.GAINS = TRUE reproduces the published bug -- ",
            "add_gains stays 0 for every sample and the HRD_Gains group will be ",
            "empty. See docs/KNOWN_ISSUES.md.", call. = FALSE)
    gen$add_gains <- 0L
  } else {
    if (length(add.present) < length(MMGC_ADD_GAIN_FEATURES)) {
      message("MMGenoCla(): additional-gain burden computed over ",
              length(add.present), " of ", length(MMGC_ADD_GAIN_FEATURES),
              " features (absent: ",
              paste(setdiff(MMGC_ADD_GAIN_FEATURES, add.present), collapse = ", "),
              ").")
    }
    gen$add_gains <- as.integer(.mmgc_feature_sum(gen, add.present) >= 3)
  }

  #### . . . - Deletion burden: losses plus tumour suppressors, 5 or more hits.
  ####         Five is 10% of the 50 loss/TSG features.
  ####
  ####         The original took the span from the first column matching "Del",
  ####         "loss" or "CNV.SNV_" to the last, and summed everything in between,
  ####         which is only the intended set because FinalMatrixGenerator() emits
  ####         those columns contiguously. Selecting them by name instead means the
  ####         burden no longer silently depends on column order.
  loss.features <- c(intersect(MMGC_FEATURES_LOSS, colnames(gen)),
                     grep("^CNV\\.SNV_", colnames(gen), value = TRUE))
  gen$deletions <- as.integer(.mmgc_feature_sum(gen, loss.features) >= 5)

  gen$x1q <- if ("CNV_Gain_Amp1q" %in% colnames(gen)) gen$CNV_Gain_Amp1q else 0

  ################################################################################
  #### . . . - Group assignment - . . . ####
  ################################################################################

  gen$clusters <- .mmgc_assign_clusters(gen)

  ################################################################################
  #### . . . - Output - . . . ####
  ################################################################################

  gen <- gen[order(gen$clusters), , drop = FALSE]

  #### . . . - Strip the feature-class prefixes for the heatmap.
  nm <- colnames(gen)
  nm <- gsub("^CNV\\.SNV_", "", nm)
  nm <- gsub("^CNV_", "", nm)
  nm <- gsub("^SNV_", "", nm)
  nm <- gsub("^SNV\\.", "", nm)
  nm[nm == "t_CCND1"] <- "CCND1-IGH"
  nm[nm == "t_MMSET"] <- "NSD2-IGH"
  nm[nm == "t_MAF"]   <- "MAF-IGH"
  colnames(gen) <- nm

  rows <- intersect(MMGC_HEATMAP_ROWS, colnames(gen))
  absent.rows <- setdiff(MMGC_HEATMAP_ROWS, colnames(gen))
  if (length(absent.rows)) {
    message("MMGenoCla(): ", length(absent.rows),
            " heatmap row(s) not in this matrix, omitted from $genomic: ",
            paste(absent.rows, collapse = ", "))
  }

  mtx <- gen[, rows, drop = FALSE]

  groups <- .mmgc_group_names(gen$clusters)
  annotation <- data.frame(clusters.new = factor(groups, levels = names(MMGC_GROUP_MAP)),
                           row.names = gen$sample)

  present.groups <- intersect(names(MMGC_GROUP_MAP), unique(groups))
  colours <- list(clusters.new = MMGC_GROUP_COLOURS[present.groups])

  if (!isTRUE(QUIET)) {
    message("MMGenoCla(): ", nrow(mtx), " samples classified into ",
            length(present.groups), " groups.")
    print(table(annotation$clusters.new, useNA = "ifany"))
  }

  invisible(list(
    genomic                = mtx,
    annotation             = annotation,
    colors                 = colours,
    clusters               = stats::setNames(gen$clusters, gen$sample),
    high.missing.features  = high.missing,
    missing.features       = missing.features))
}


#### . . . - .mmgc_feature_sum() - . . . ####
#### . . . - Row sums over the named features that are actually present. Returns a
####         zero vector when none are, which keeps the caller's comparison the
####         right length -- the published code's arithmetic over an absent column
####         collapsed to length zero and silently disabled the feature.

.mmgc_feature_sum <- function(X, FEATURES) {
  have <- intersect(FEATURES, colnames(X))
  if (!length(have)) return(rep(0, nrow(X)))
  rowSums(as.data.frame(lapply(X[, have, drop = FALSE], as.numeric)), na.rm = TRUE)
}


#### . . . - .mmgc_assign_clusters() - . . . ####
#### . . . - The ordered rule set. Each rule may overwrite the previous assignment,
####         so the sequence is part of the definition -- it is kept exactly as
####         published.

.mmgc_assign_clusters <- function(g) {

  cl <- rep("0", nrow(g))
  f  <- function(name) if (name %in% colnames(g)) as.numeric(g[[name]]) else rep(0, nrow(g))

  HRD <- g$HRD; add_gains <- g$add_gains; deletions <- g$deletions; x1q <- g$x1q
  t_MMSET <- f("t_MMSET"); t_CCND1 <- f("t_CCND1"); t_MAF <- f("t_MAF")
  APOBEC <- f("APOBEC"); CNV.Sig <- f("CNV.Sig")
  HS_DIS3 <- f("SNV.HS_DIS3")
  TP53 <- f("CNV.SNV_TP53"); RB1 <- f("CNV.SNV_RB1"); TGDS <- f("CNV.SNV_TGDS")
  no.tr <- t_MMSET == 0 & t_CCND1 == 0 & t_MAF == 0

  #### . . . - Hyperdiploid, MAPK-mutant, otherwise quiet
  ras <- Reduce(`|`, lapply(MMGC_RAS_FEATURES, function(x) f(x) > 0))
  cl[HRD == 1 & ras & add_gains == 0 & deletions == 0 & no.tr] <- "A_HRD_RAS"

  #### . . . - Hyperdiploid with extra gains, or with a deletion burden
  cl[HRD == 1 & deletions == 0 & add_gains == 1 & cl == "0" & no.tr] <- "B_HRD_gains"
  cl[HRD == 1 & deletions == 1 & add_gains == 0 & cl == "0" & no.tr] <- "C_HRD_deletions"

  #### . . . - Remaining hyperdiploid samples fall back to HRD_RAS
  cl[HRD == 1 & cl == "0" & add_gains == 0 & no.tr] <- "A_HRD_RAS"

  #### . . . - MAF translocation
  cl[t_MAF > 0 & HS_DIS3 == 0 & APOBEC > 1] <- "I_MAF_APOBEC"
  cl[t_MAF > 0 & HS_DIS3 == 0 & APOBEC < 2] <- "I_MAF_APOBEC"

  #### . . . - CCND1 translocation
  cl[t_CCND1 > 0 & HRD == 1 & deletions == 0] <- "E_CCND1_HRD"
  for (x in MMGC_RAS_FEATURES) {
    cl[t_CCND1 > 0 & f(x) > 0 & deletions == 0] <- "F_CCND1_RAS"
  }
  cl[t_CCND1 > 0 & deletions > 0] <- "G_CCND1_complex"
  cl[t_CCND1 > 0 & cl == "0" & deletions == 0] <- "FF_CCND1_RAS"
  ccnd1.simple <- c("E_CCND1_HRD", "F_CCND1_RAS", "FF_CCND1_RAS")
  cl[t_CCND1 > 0 & cl %in% ccnd1.simple & x1q >= 1] <- "G_CCND1_complex"
  cl[t_CCND1 > 0 & cl %in% ccnd1.simple & TP53 > 1]  <- "G_CCND1_complex"

  #### . . . - NSD2 translocation and/or 1q with del(13q)
  del13q <- RB1 > 0 | TGDS > 0
  cl[t_MMSET >  0 & x1q >  0 & del13q & HS_DIS3 == 0 & cl == "0"] <- "N_NSD2_1q_del13q"
  cl[t_MMSET >  0 & x1q == 0 & del13q & HS_DIS3 == 0 & cl == "0"] <- "P_NSD2_del13q"
  cl[t_MMSET == 0 & x1q >  0 & del13q & HS_DIS3 == 0 & cl == "0"] <- "Q_x1q_del13q"

  #### . . . - Hyper-APOBEC
  cl[APOBEC > 1 & t_MAF   >  0 & cl == "0"] <- "L_MAF_APOBEC"
  cl[APOBEC > 1 & t_MMSET == 0 & cl == "0"] <- "L_MAF_APOBEC"
  cl[APOBEC > 1 & t_CCND1 == 0 & cl == "0"] <- "L_MAF_APOBEC"
  cl[APOBEC > 1 & cl == "Q_x1q_del13q"]     <- "L_MAF_APOBEC"
  cl[APOBEC > 1 & cl %in% c("A_HRD_RAS", "B_HRD_gains", "BB_HRD_gains",
                            "C_HRD_deletions", "CC_HRD_deletions",
                            "D_HRD_other", "DD_HRD_other",
                            "0_HRD_Simple")] <- "M_MAF_APOBEC"

  #### . . . - Remaining NSD2
  cl[t_MMSET > 0 & cl == "0"] <- "O_NSD2"

  #### . . . - Everything else, by deletion burden
  cl[cl == "0" & deletions == 1] <- "R_Multiple_losses"
  cl[cl == "0"]                  <- "RR_Simple"

  #### . . . - Late reassignments
  cl[cl %in% ccnd1.simple & (APOBEC > 1 | CNV.Sig == 1)] <- "G_CCND1_complex"
  cl[cl %in% c("M_MAF_APOBEC", "I_MAF_APOBEC", "L_MAF_APOBEC") &
       t_MMSET > 0 & t_MAF > 0] <- "N_NSD2_1q_del13q"
  cl[cl == "RR_Simple" & (CNV.Sig == 1 | x1q == 2)] <- "R_Multiple_losses"

  cl
}


#### . . . - .mmgc_group_names() - . . . ####
#### . . . - Internal rule labels to published group names.

.mmgc_group_names <- function(CLUSTERS) {
  out <- rep(NA_character_, length(CLUSTERS))
  for (g in names(MMGC_GROUP_MAP)) {
    out[CLUSTERS %in% MMGC_GROUP_MAP[[g]]] <- g
  }
  if (any(is.na(out))) {
    unmapped <- sort(unique(CLUSTERS[is.na(out)]))
    warning(".mmgc_group_names(): rule label(s) with no published group: ",
            paste(unmapped, collapse = ", "), call. = FALSE)
  }
  out
}
