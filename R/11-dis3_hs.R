#### . . . - DIS3 hotspots - . . . ####
#
# DIS3 mutations split into two features because they do not behave alike: the
# three RNB-domain hotspots (D479, D488, R780) sit in the catalytic site and
# abolish exonuclease activity, while other DIS3 mutations do not carry the same
# weight in the classification.

MMGC_DIS3_HOTSPOTS <- c("D479", "D488", "R780")


#### . . . - dis3_hs() - . . . ####

#### . . . - Split SNV_DIS3 into SNV.HS_DIS3 (hotspot) and SNV.noHS_DIS3 (any
####         other non-synonymous DIS3 mutation).
####
####         MATRIX     = final matrix from FinalMatrixGenerator(), or any table
####                      with `sample` and `SNV_DIS3` columns
####         CODING.MUT = the same non-synonymous mutation table given to
####                      MutationsMatrixGenerator(); needs sampleID, gene,
####                      impact and aachange
####         HOTSPOTS   = residues to treat as hotspots (default MMGC_DIS3_HOTSPOTS)
####
####         Returns sample, SNV.HS_DIS3, SNV.noHS_DIS3. Samples with no DIS3 call
####         at all stay NA in both.

dis3_hs <- function(MATRIX, CODING.MUT, HOTSPOTS = MMGC_DIS3_HOTSPOTS) {

  mtx <- as.data.frame(MATRIX, stringsAsFactors = FALSE)
  if (!"sample" %in% colnames(mtx)) {
    stop("dis3_hs(): MATRIX has no 'sample' column.", call. = FALSE)
  }
  if (!"SNV_DIS3" %in% colnames(mtx)) {
    stop("dis3_hs(): MATRIX has no 'SNV_DIS3' column. Run ",
         "FinalMatrixGenerator() first.", call. = FALSE)
  }

  annot <- read.delim(CODING.MUT, stringsAsFactors = FALSE)
  for (col in c("sampleID", "gene", "impact")) {
    if (!col %in% colnames(annot)) {
      stop("dis3_hs(): ", CODING.MUT, " has no column ", sQuote(col), call. = FALSE)
    }
  }

  dis3 <- annot[which(annot$gene == "DIS3" & !is.na(annot$sampleID) &
                        annot$impact != "Synonymous"), , drop = FALSE]

  if (!nrow(dis3)) {
    message("dis3_hs(): no DIS3 mutations in this cohort. ",
            "SNV.HS_DIS3 is 0 throughout and SNV.noHS_DIS3 copies SNV_DIS3.")
    out <- data.frame(sample        = mtx$sample,
                      SNV.HS_DIS3   = ifelse(is.na(mtx$SNV_DIS3), NA_real_, 0),
                      SNV.noHS_DIS3 = mtx$SNV_DIS3,
                      stringsAsFactors = FALSE)
    return(out)
  }

  if (!"aachange" %in% colnames(annot)) {
    stop("dis3_hs(): ", CODING.MUT, " has DIS3 mutations but no 'aachange' ",
         "column, so hotspots cannot be identified.", call. = FALSE)
  }

  #### . . . - A sample is hotspot-mutant if any of its DIS3 mutations hits one of
  ####         the listed residues.
  pattern <- paste(HOTSPOTS, collapse = "|")
  hs.samples <- unique(dis3$sampleID[grepl(pattern, dis3$aachange)])

  out <- data.frame(
    sample      = mtx$sample,
    SNV.HS_DIS3 = ifelse(mtx$sample %in% hs.samples, 1,
                         ifelse(is.na(mtx$SNV_DIS3), NA_real_, 0)),
    #### . . . - A hotspot carrier is counted once, under the hotspot feature only.
    SNV.noHS_DIS3 = ifelse(mtx$sample %in% hs.samples, 0, mtx$SNV_DIS3),
    stringsAsFactors = FALSE)

  message("dis3_hs(): ", length(hs.samples), " sample(s) with a DIS3 hotspot (",
          paste(HOTSPOTS, collapse = ", "), ")")
  out
}
