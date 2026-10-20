#####################################################################################
####                                                                             ####
####   MMGenomicClassification -- mm10 (mouse) pipeline                          ####
####                                                                             ####
####   For the vkMYC and related mouse models of multiple myeloma.               ####
####                                                                             ####
####   Load AFTER the main entry point:                                          ####
####                                                                             ####
####       source("MMGenomicClassification.R")                                   ####
####       source("R/mm10/MMGenomicClassification_mm10.R")                       ####
####       MMGC_setup("mm10")                                                    ####
####       MMGC_mm10_report()      # what transfers and what does not            ####
####                                                                             ####
#####################################################################################
#
# WHAT THIS IS
#
# The classification was derived in human myeloma. Applying it to mouse is a
# translation, not a re-derivation, and the translation is only sound for part of
# it. The design principle, inherited from how reference/mm10 was built, is that
# every feature keeps its HUMAN name and meaning while being measured at the
# ORTHOLOGOUS MOUSE LOCUS. That way the eleven groups stay comparable and
# MMGenoCla() needs no mouse-specific rules.
#
# Three things transfer cleanly, one does not:
#
#   TRANSFERS   Recurrent copy-number peaks. reference/mm10 relocates each human
#               GISTIC peak onto its target gene's mouse ortholog and keeps the
#               human peak code, so Amp_1q21.1 still means "MCL1 gain" -- measured
#               on mouse chr3, where Mcl1 lives. 55 of the 64 human peaks have a
#               confident mouse anchor.
#
#   TRANSFERS   Driver gene mutations and deletions, via one-to-one orthology.
#               All 32 tumour-suppressor features and 51 of the 59 mutation
#               features have a mouse ortholog.
#
#   TRANSFERS   Deletion burden and the 1q / del(13q) axis, since both are built
#               from the features above.
#
#   DOES NOT    Hyperdiploidy. Human HRD is trisomy of chromosomes 3, 5, 7, 9,
#               11, 15, 19 and 21. Mouse has no counterpart: synteny is broken, and
#               each of those human chromosomes scatters its genes over several
#               mouse chromosomes -- the best single proxy carries only 33-76% of
#               the orthologs, and two pairs of human chromosomes collide onto the
#               same mouse chromosome. HRD.PROXY therefore defaults to "none":
#               the HRD features are set to 0, no HRD-based group is assigned, and
#               the classification runs on the features that do transfer.
#               Run MMGC_mm10_report() to see the numbers.
#
# Also note: the APOBEC feature is scored against cut-offs from a 1,933-sample
# human cohort, and the IGH translocations (t_CCND1 / t_MMSET / t_MAF) do not occur
# in mouse models the way they do in human myeloma. Supply them as 0 unless your
# model has an engineered equivalent.
#
# Every function here is a thin mouse-specific layer over the shared
# implementation in R/. It sets up the ortholog translation, the mouse chromosome
# coding and the reduced feature set, then delegates -- rather than carrying a
# third copy of 2,500 lines that would drift out of step with the human version.


if (!exists("MMGC_setup")) {
  stop("Source MMGenomicClassification.R before this file.", call. = FALSE)
}


#### . . . - Mouse-specific constants - . . . ####

#### . . . - Human features with no one-to-one mouse ortholog: the histone cluster
####         genes and HLA-C, which are gene families rather than single orthologs.
####         Five of these are excluded from the classification anyway; H1.4, H1.2
####         and HLA.C are lost from the heatmap only -- no cluster rule uses them.
MMGC_MM10_NO_ORTHOLOG <- paste0("SNV_", c("H1.4", "H1.5", "H1.2", "H1.3",
                                          "H2BC12", "HLA.C", "LCE1D", "TBC1D29P"))

#### . . . - Human GISTIC peaks with no confident mouse anchor (see
####         reference/mm10/QC_unmapped_peaks.txt). Only Del_12q24.31 survives to
####         the final feature set; the rest are dropped by the human pipeline too.
MMGC_MM10_NO_PEAK <- c("CNV_Del_12q24.31")

#### . . . - The mouse feature set: the published one, minus what cannot be measured.
MMGC_FINAL_FEATURES_MM10 <- setdiff(MMGC_FINAL_FEATURES,
                                    c(MMGC_MM10_NO_ORTHOLOG, MMGC_MM10_NO_PEAK))


#### . . . - MMGC_mm10_orthologs() - . . . ####
#### . . . - Mouse symbol -> human symbol, from the `gene_name` / `human_gene_name`
####         columns of Driver_summary_mm10.txt. Cached after the first call.

MMGC_mm10_orthologs <- function(REFRESH = FALSE) {
  if (!REFRESH && !is.null(.MMGC$mm10.orthologs)) return(.MMGC$mm10.orthologs)

  f <- MMGC_ref("drivers")
  d <- read.delim(f, stringsAsFactors = FALSE)
  if (!all(c("gene_name", "human_gene_name") %in% colnames(d))) {
    stop("MMGC_mm10_orthologs(): ", f, " needs columns 'gene_name' (mouse) and ",
         "'human_gene_name'.", call. = FALSE)
  }
  d <- d[!is.na(d$gene_name) & !is.na(d$human_gene_name) &
           nzchar(d$gene_name) & nzchar(d$human_gene_name), , drop = FALSE]

  map <- stats::setNames(MMGC_canonical_symbols(d$human_gene_name), d$gene_name)
  map <- map[!duplicated(names(map))]

  .MMGC$mm10.orthologs <- map
  map
}


#### . . . - MMGC_mm10_to_human() - . . . ####
#### . . . - Rewrite mouse gene symbols as their human orthologs, so every feature
####         downstream carries its human name. Symbols with no ortholog are left
####         as they are -- they simply never match a feature.
####
####         Matching ignores case, because mouse symbols are title case (Trp53)
####         and the gene-model and mutation tables are not always consistent.

MMGC_mm10_to_human <- function(X) {
  if (!length(X)) return(X)
  map <- MMGC_mm10_orthologs()
  idx <- match(toupper(as.character(X)), toupper(names(map)))
  out <- as.character(X)
  out[!is.na(idx)] <- unname(map[idx[!is.na(idx)]])
  out
}


#### . . . - MMGC_mm10_report() - . . . ####
#### . . . - Print what the mouse translation can and cannot measure. Worth running
####         before any real analysis.

MMGC_mm10_report <- function() {
  if (!identical(MMGC_genome(), "mm10")) {
    stop("MMGC_mm10_report(): call MMGC_setup(\"mm10\") first.", call. = FALSE)
  }

  peaks <- read.delim(MMGC_ref("gistic_peaks"), stringsAsFactors = FALSE)
  map   <- MMGC_mm10_orthologs()
  hrd   <- MMGC_mm10_hrd_map()

  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("mm10 translation report\n")
  cat(strrep("=", 74), "\n\n", sep = "")

  cat("Copy-number peaks\n")
  cat("  relocated to a mouse ortholog: ", nrow(peaks), " of 64 human peaks\n", sep = "")
  cat("  features lost as a result:     ",
      if (length(MMGC_MM10_NO_PEAK)) paste(MMGC_MM10_NO_PEAK, collapse = ", ") else "none",
      "\n\n", sep = "")

  cat("Driver genes\n")
  cat("  mouse-human ortholog pairs:    ", length(map), "\n", sep = "")
  cat("  features with no ortholog:     ",
      paste(sub("^SNV_", "", MMGC_MM10_NO_ORTHOLOG), collapse = ", "), "\n", sep = "")
  cat("  (of those, only H1.4, H1.2 and HLA.C affect anything, and only the\n")
  cat("   heatmap -- no classification rule uses them)\n\n")

  cat("Feature set\n")
  cat("  published features:            ", length(MMGC_FINAL_FEATURES) - 1L, "\n", sep = "")
  cat("  measurable in mouse:           ", length(MMGC_FINAL_FEATURES_MM10) - 1L, "\n\n", sep = "")

  cat("Hyperdiploidy -- DOES NOT TRANSFER\n")
  for (i in seq_len(nrow(hrd))) {
    cat(sprintf("  human chr%-2s -> mouse chr%-3s %5.1f%% of orthologs%s\n",
                hrd$human_chrom[i], hrd$mouse_chrom[i],
                100 * hrd$ortholog_fraction[i],
                if (sum(hrd$mouse_chrom == hrd$mouse_chrom[i]) > 1)
                  "   <- shared with another chromosome" else ""))
  }
  cat("\n  HRD.PROXY defaults to \"none\": these features are set to 0 and no\n")
  cat("  HRD-based group is assigned. Pass HRD.PROXY = \"ortholog\" only if the\n")
  cat("  fractions above are good enough for your question.\n\n")

  cat("Other caveats\n")
  cat("  - APOBEC is binned on cut-offs from a 1,933-sample human cohort.\n")
  cat("  - IGH translocations do not occur in mouse models as they do in human\n")
  cat("    myeloma; supply t_CCND1 / t_MMSET / t_MAF as 0 unless engineered.\n")
  cat("  - DIS3 hotspots D479/D488/R780 are HUMAN residue numbers. Confirm the\n")
  cat("    mouse Dis3 alignment before relying on dis3_hs_mm10().\n")
  cat(strrep("=", 74), "\n\n", sep = "")

  invisible(list(peaks = nrow(peaks), orthologs = length(map), hrd_map = hrd,
                 features = length(MMGC_FINAL_FEATURES_MM10) - 1L))
}


#### . . . - MMGC_mm10_hrd_map() - . . . ####
#### . . . - The reviewable human-HRD-chromosome to mouse-chromosome map.

MMGC_mm10_hrd_map <- function() {
  f <- MMGC_ref("hrd_map")
  d <- read.delim(f, stringsAsFactors = FALSE, comment.char = "#")
  for (col in c("human_chrom", "mouse_chrom", "ortholog_fraction", "use")) {
    if (!col %in% colnames(d)) {
      stop("MMGC_mm10_hrd_map(): ", f, " has no column ", sQuote(col),
           ".\n  Rebuild it with: Rscript scripts/build_mm10_arms_and_hrd_map.R",
           call. = FALSE)
    }
  }
  d
}


################################################################################
####                                                                        ####
####   Pipeline                                                             ####
####                                                                        ####
################################################################################

#### . . . - inputCNV_mm10() - . . . ####
#### . . . - As inputCNV(), with mouse chromosome coding (1-19, X = 20, Y = 21).

inputCNV_mm10 <- function(DATA, SAMPLE.COL, CHR.COL, START.COL, END.COL,
                          TOTAL.COL, MINOR.COL) {
  inputCNV(DATA, SAMPLE.COL = SAMPLE.COL, CHR.COL = CHR.COL, START.COL = START.COL,
           END.COL = END.COL, TOTAL.COL = TOTAL.COL, MINOR.COL = MINOR.COL,
           GENOME = "mm10")
}


#### . . . - FocalLargeMatrix_mm10() - . . . ####
#### . . . - As FocalLargeMatrix(). The mouse GISTIC table keeps the human peak
####         codes, so no translation is needed here.

FocalLargeMatrix_mm10 <- function(CNV, DATABASE, OUT.DIR, FOCAL.LENGTH = 5000000,
                                  GISTIC.FILE = NULL) {
  FocalLargeMatrix(CNV = CNV, GISTIC.FILE = GISTIC.FILE, DATABASE = DATABASE,
                   OUT.DIR = OUT.DIR, FOCAL.LENGTH = FOCAL.LENGTH)
}


#### . . . - HYidentification_mm10() - . . . ####

#### . . . - Whole-chromosome gain in mouse, reported under the human HRD
####         chromosome labels the downstream matrix expects.
####
####         HRD.PROXY = "none"      HRD features are 0 for every sample. The
####                                 classification runs on the transferable
####                                 features only. This is the default, because
####                                 mouse has no equivalent of human trisomy.
####         HRD.PROXY = "ortholog"  Use reference/mm10/hrd_chromosome_map.txt:
####                                 a human HRD chromosome is called gained when
####                                 its mapped mouse chromosome is gained. Read
####                                 MMGC_mm10_report() first -- the best proxy
####                                 carries only a third of the orthologs for
####                                 some chromosomes, and two pairs collide.
####
####         Mouse chromosomes are telocentric, so reference/mm10/arm_mm10.txt
####         gives each one a single q arm spanning its whole length. Arm-level and
####         chromosome-level gain are therefore the same measurement in mouse.

HYidentification_mm10 <- function(CNV, DATABASE, OUT.DIR,
                                  PERCENTAGE = 0.6,
                                  HRD.PROXY = c("none", "ortholog"),
                                  ARMS = NULL, CHROM.LENGHT = NULL) {

  HRD.PROXY <- match.arg(HRD.PROXY)

  #### . . . - Score gain on the MOUSE chromosomes first.
  mouse.chroms <- 1:19
  res <- HYidentification(CNV = CNV, ARMS = ARMS, CHROM.LENGHT = CHROM.LENGHT,
                          DATABASE = DATABASE, OUT.DIR = OUT.DIR,
                          PERCENTAGE = PERCENTAGE,
                          HRD.CHROMOSOMES = mouse.chroms,
                          EXTRA.CHROMOSOMES = integer(0),
                          DOUBLE.GENOME.ARMS = 10)

  mouse.hyper <- res$hyperdiploid
  human.hrd <- c(3, 5, 7, 9, 11, 15, 19, 21, 18)

  out <- data.frame(sample = mouse.hyper$sample,
                    type   = if (nrow(mouse.hyper)) mouse.hyper$type else character(0),
                    whole.median.CN = if (nrow(mouse.hyper)) mouse.hyper$whole.median.CN
                                      else numeric(0),
                    stringsAsFactors = FALSE)

  if (identical(HRD.PROXY, "none")) {
    #### . . . - Zero, not NA: MMGenoCla() drops any sample with a missing feature,
    ####         so NA here would discard the whole cohort.
    for (h in human.hrd) out[[paste0("hyp.chr", h)]] <- rep(0L, nrow(out))
    message("HYidentification_mm10(): HRD.PROXY = \"none\" -- human HRD features ",
            "set to 0. No HRD-based group will be assigned.")
  } else {
    map <- MMGC_mm10_hrd_map()
    warning("HYidentification_mm10(): HRD.PROXY = \"ortholog\" maps human ",
            "trisomies onto single mouse chromosomes that carry only ",
            paste0(round(100 * range(map$ortholog_fraction, na.rm = TRUE)),
                   collapse = "-"),
            "% of the relevant orthologs. Run MMGC_mm10_report() and read ",
            "docs/MOUSE.md before using these HRD calls.", call. = FALSE)
    for (h in human.hrd) {
      row <- map[map$human_chrom == h, , drop = FALSE]
      col <- paste0("hyp.chr", h)
      if (!nrow(row) || is.na(row$mouse_chrom[1]) || row$use[1] != 1L) {
        out[[col]] <- rep(0L, nrow(out))
        next
      }
      mcol <- paste0("hyp.chr", row$mouse_chrom[1])
      v <- if (mcol %in% colnames(mouse.hyper)) mouse.hyper[[mcol]] else NA_integer_
      v[is.na(v)] <- 0L
      out[[col]] <- as.integer(v)
    }
  }

  #### . . . - The HRD call itself, on the human labels.
  hrd.only <- paste0("hyp.chr", setdiff(human.hrd, 18))
  out$hyp.sample <- if (nrow(out))
    ifelse(rowSums(out[, hrd.only, drop = FALSE], na.rm = TRUE) >= 2,
           "Hyperdiploidy", "No.HRD") else character(0)
  out$Cytogenetic <- if (nrow(out)) mouse.hyper$Cytogenetic else character(0)

  #### . . . - Overwrite the human-labelled table; keep the mouse-chromosome one
  ####         alongside it, since that is the measurement actually made.
  tag <- paste0(PERCENTAGE * 100, "%_ARMS_", DATABASE)
  f.human <- .mmgc_outfile(OUT.DIR, tag, "hyperdiploid.txt")
  f.mouse <- .mmgc_outfile(OUT.DIR, tag, "hyperdiploid_mouse_chromosomes.txt")

  write.table(mouse.hyper, f.mouse, col.names = TRUE, row.names = FALSE,
              quote = FALSE, sep = "\t")
  write.table(out, f.human, col.names = TRUE, row.names = FALSE,
              quote = FALSE, sep = "\t")

  message("HYidentification_mm10(): mouse-chromosome gains in ", basename(f.mouse),
          "; human-labelled features in ", basename(f.human))

  invisible(list(hyperdiploid = out,
                 mouse_chromosomes = mouse.hyper,
                 no.hyperdiploid = res$no.hyperdiploid))
}


#### . . . - CytogeneticMatrixGenerator_mm10() - . . . ####
#### . . . - Unchanged: the mouse peak codes are the human ones.

CytogeneticMatrixGenerator_mm10 <- function(FocalLargeMATRIX, HY.FILE, DATABASE, OUT.DIR) {
  CytogeneticMatrixGenerator(FocalLargeMATRIX = FocalLargeMATRIX, HY.FILE = HY.FILE,
                             DATABASE = DATABASE, OUT.DIR = OUT.DIR)
}


#### . . . - MutationsMatrixGenerator_mm10() - . . . ####

#### . . . - As MutationsMatrixGenerator(), but every mouse gene symbol -- in the
####         mutation calls, the driver list and the gene models -- is rewritten to
####         its human ortholog first. Coordinates stay mouse; only the labels
####         become human, so the feature names match the published set.
####
####         CODING.MUT must use mouse symbols in its `gene` column (as dndscv
####         against a mouse reference produces).

MutationsMatrixGenerator_mm10 <- function(CODING.MUT, METADATA, CNV, DATABASE, OUT.DIR,
                                          UCSC.gene = NULL, DRIVERs = NULL,
                                          GisticOtherAllPeaks = NULL) {

  if (is.null(UCSC.gene)) UCSC.gene <- MMGC_ref("gene_models")
  if (is.null(DRIVERs))   DRIVERs   <- MMGC_ref("drivers")

  work <- file.path(tempdir(), paste0("mmgc_mm10_", DATABASE))
  dir.create(work, showWarnings = FALSE, recursive = TRUE)

  #### . . . - Driver list: swap `gene_name` (mouse) for the human ortholog, so the
  ####         shared code produces human-named features.
  drv <- read.delim(DRIVERs, stringsAsFactors = FALSE)
  drv$gene_name <- MMGC_canonical_symbols(drv$human_gene_name)
  drv <- drv[!is.na(drv$gene_name) & nzchar(drv$gene_name), , drop = FALSE]
  drv <- drv[!duplicated(drv$gene_name), , drop = FALSE]
  drv.file <- file.path(work, "drivers_humanised.txt")
  write.table(drv, drv.file, sep = "\t", quote = FALSE, row.names = FALSE)

  #### . . . - Mutation calls: mouse gene symbol -> human.
  mut <- read.delim(CODING.MUT, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(mut)) {
    stop("MutationsMatrixGenerator_mm10(): ", CODING.MUT, " has no 'gene' column.",
         call. = FALSE)
  }
  n.before <- length(unique(mut$gene))
  mut$gene <- MMGC_mm10_to_human(mut$gene)
  mut.file <- file.path(work, "mutations_humanised.txt")
  write.table(mut, mut.file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("MutationsMatrixGenerator_mm10(): ", n.before,
          " mouse gene symbols translated to human orthologs where known.")

  #### . . . - Gene models: keep the mouse coordinates, relabel with human symbols.
  gm <- read.delim(UCSC.gene, stringsAsFactors = FALSE)
  gm$geneSymbol <- MMGC_mm10_to_human(gm$geneSymbol)
  gm.file <- file.path(work, "gene_models_humanised.txt")
  write.table(gm, gm.file, sep = "\t", quote = FALSE, row.names = FALSE)

  MutationsMatrixGenerator(CODING.MUT = mut.file, METADATA = METADATA,
                           UCSC.gene = gm.file, DRIVERs = drv.file,
                           GisticOtherAllPeaks = GisticOtherAllPeaks,
                           CNV = CNV, DATABASE = DATABASE, OUT.DIR = OUT.DIR,
                           CHRX.GENES = MMGC_CHRX_GENES)
}


#### . . . - APOBECanalyzer_mm10() - . . . ####

#### . . . - Signature fitting against mm10 trinucleotide contexts.
####
####         The SBS reference signatures are human. Fitting them to mouse mutations
####         is standard practice for APOBEC, whose CpG-independent TCW context is
####         conserved, but the absolute contributions are not directly comparable
####         with human cohorts -- and MMGC_APOBEC_CUTOFFS are human deciles. Treat
####         the resulting APOBEC feature as a rank within your mouse cohort.

APOBECanalyzer_mm10 <- function(SNV.OBJECT, SAMPLE, CHROM, POS, REF, ALT,
                                ITER = 10, DATABASE, OUT.DIR,
                                STRAND.BIAS = TRUE, PLOT = TRUE) {
  message("APOBECanalyzer_mm10(): fitting human SBS reference signatures to mouse ",
          "contexts. The APOBEC bins use human cohort cut-offs -- see docs/MOUSE.md.")
  APOBECanalyzer(SNV.OBJECT = SNV.OBJECT, SAMPLE = SAMPLE, CHROM = CHROM,
                 POS = POS, REF = REF, ALT = ALT, ITER = ITER,
                 DATABASE = DATABASE, OUT.DIR = OUT.DIR, GENOME = "mm10",
                 STRAND.BIAS = STRAND.BIAS, PLOT = PLOT)
}


#### . . . - FinalMatrixGenerator_mm10() - . . . ####
#### . . . - As FinalMatrixGenerator(), with the reduced mouse feature set.

FinalMatrixGenerator_mm10 <- function(INPUT.FILE.FOLDER, DATABASE, OUTPUT,
                                      FEATURES = MMGC_FINAL_FEATURES_MM10,
                                      STRICT = FALSE) {
  FinalMatrixGenerator(INPUT.FILE.FOLDER = INPUT.FILE.FOLDER, DATABASE = DATABASE,
                       OUTPUT = OUTPUT, FEATURES = FEATURES, STRICT = STRICT)
}


#### . . . - dis3_hs_mm10() - . . . ####

#### . . . - DIS3 hotspots in mouse.
####
####     WARNING: D479, D488 and R780 are HUMAN DIS3 residue numbers. Mouse Dis3 is
####     not guaranteed to number the same residues identically. Check the
####     alignment for your annotation and pass HOTSPOTS explicitly if it differs.

dis3_hs_mm10 <- function(MATRIX, CODING.MUT, HOTSPOTS = MMGC_DIS3_HOTSPOTS) {

  warning("dis3_hs_mm10(): hotspot residues ", paste(HOTSPOTS, collapse = ", "),
          " are HUMAN DIS3 positions. Confirm they correspond to the same ",
          "residues in your mouse Dis3 annotation.", call. = FALSE)

  #### . . . - Translate the gene column so "Dis3" matches the human "DIS3" the
  ####         shared function looks for.
  mut <- read.delim(CODING.MUT, stringsAsFactors = FALSE)
  mut$gene <- MMGC_mm10_to_human(mut$gene)
  work <- file.path(tempdir(), paste0("mmgc_mm10_dis3_", basename(CODING.MUT)))
  write.table(mut, work, sep = "\t", quote = FALSE, row.names = FALSE)

  dis3_hs(MATRIX = MATRIX, CODING.MUT = work, HOTSPOTS = HOTSPOTS)
}


#### . . . - MMGenoCla_mm10() - . . . ####

#### . . . - The classification, on the mouse feature set.
####
####         Identical rules to MMGenoCla() -- that is the point of keeping human
####         feature names. With HRD.PROXY = "none" upstream, the HRD features are
####         all 0, so the HRD_RAS / HRD_Gains / HRD_Complex_Cytogenetic groups
####         are unreachable and samples resolve into the translocation, 1q /
####         del(13q), Multiple_Losses and Simple groups instead.

MMGenoCla_mm10 <- function(input, CONFIRM = FALSE, STRICT = FALSE,
                           LEGACY.ADD.GAINS = FALSE,
                           EXPECTED = MMGC_FINAL_FEATURES_MM10,
                           QUIET = FALSE) {

  res <- MMGenoCla(input = input, CONFIRM = CONFIRM, STRICT = STRICT,
                   LEGACY.ADD.GAINS = LEGACY.ADD.GAINS, EXPECTED = EXPECTED,
                   QUIET = QUIET)
  if (is.null(res)) return(invisible(NULL))

  hrd.assigned <- res$annotation$clusters.new %in%
    c("HRD_RAS", "HRD_Gains", "HRD_Complex_Cytogenetic")
  if (any(hrd.assigned)) {
    message("MMGenoCla_mm10(): ", sum(hrd.assigned), " sample(s) in an HRD group. ",
            "These rest on the ortholog HRD proxy -- see MMGC_mm10_report().")
  } else {
    message("MMGenoCla_mm10(): no HRD groups assigned, as expected with ",
            "HRD.PROXY = \"none\".")
  }
  invisible(res)
}


message("mm10 pipeline loaded. Run MMGC_mm10_report() after MMGC_setup(\"mm10\") ",
        "to see what transfers.")
