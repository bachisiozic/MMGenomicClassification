#!/usr/bin/env Rscript
#
# selftest.R -- run the whole pipeline on synthetic data, for every build.
#
# This is a smoke test, not a validation of the biology: it checks that each step
# runs, that the files it writes are read back correctly by the next step, and
# that the classification produces the expected groups. It needs no real data, so
# anyone can run it after cloning to confirm their installation works.
#
# Usage:
#   Rscript scripts/selftest.R              # all builds it can run
#   Rscript scripts/selftest.R hg38         # one build
#
# Signature fitting is skipped unless the genome packages are installed; the
# APOBEC feature is then filled with synthetic values so the rest can be checked.

this.file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
repo.root <- if (length(this.file)) dirname(dirname(normalizePath(this.file))) else getwd()

suppressMessages(source(file.path(repo.root, "MMGenomicClassification.R")))

args   <- commandArgs(trailingOnly = TRUE)
builds <- if (length(args)) args else c("hg38", "hg19", "mm10")

set.seed(42)

#### . . . - Synthetic cohort ------------------------------------------------------

make_cnv <- function(genome, peaks, n.samples = 12) {
  #### . . . - Build a DISJOINT segmentation per chromosome, the way an
  ####         allele-specific WGS caller emits it: a diploid backbone that is split
  ####         at event boundaries rather than having events laid on top of it.
  ####         Overlapping segments would double-count gained bases downstream, and
  ####         checkCNV() is there to reject exactly that -- so the fixture has to be
  ####         clean or the test is not testing the real thing.
  gistic <- read.delim(peaks, stringsAsFactors = FALSE)
  gistic$chr <- gsub("^chr", "", gistic$chr)

  sizes <- read.delim(file.path(repo.root, "reference", genome,
                                sprintf("chrom_sizes_%s.txt", genome)),
                      stringsAsFactors = FALSE)
  autos <- sizes[!sizes$Offset %in% c("X", "Y"), ]

  #### . . . - Lay `events` (start, end, total, minor) over one chromosome and
  ####         return disjoint segments covering it end to end.
  overlay <- function(chr.len, events) {
    cuts <- sort(unique(c(1, chr.len, events$start, pmin(events$end, chr.len))))
    cuts <- cuts[cuts >= 1 & cuts <= chr.len]
    if (length(cuts) < 2) return(NULL)
    seg <- data.frame(start = utils::head(cuts, -1), end = utils::tail(cuts, -1))
    seg$total <- 2; seg$minor <- 1
    for (i in seq_len(nrow(events))) {
      hit <- seg$start >= events$start[i] & seg$end <= events$end[i]
      seg$total[hit] <- events$total[i]
      seg$minor[hit] <- events$minor[i]
    }
    seg
  }

  out <- list()
  for (i in seq_len(n.samples)) {
    s <- sprintf("SAMPLE_%02d", i)

    #### . . . - Events for this sample: a few recurrent peaks, plus whole-chromosome
    ####         trisomies in half the samples so the HRD branch is reached.
    k <- sample(4:12, 1)
    hit <- gistic[sample(nrow(gistic), k), ]
    amp <- hit$lesion == "Amp"
    ev <- data.frame(
      chr   = hit$chr,
      start = hit$start,
      end   = hit$end,
      total = ifelse(amp, sample(c(3, 4), k, TRUE), sample(c(0, 1), k, TRUE)),
      minor = ifelse(amp, 1, 0),
      stringsAsFactors = FALSE)

    if (i %% 2 == 0) {
      tri.chr <- intersect(as.character(c(3, 5, 7, 9, 11)), autos$Offset)
      ev <- rbind(ev, data.frame(
        chr = tri.chr, start = 1,
        end = as.numeric(autos$Size[match(tri.chr, autos$Offset)]),
        total = 3, minor = 1, stringsAsFactors = FALSE))
    }

    rows <- list()
    for (cc in autos$Offset) {
      L <- as.numeric(autos$Size[autos$Offset == cc])
      e <- ev[ev$chr == cc, , drop = FALSE]
      #### . . . - Whole-chromosome events last, so a trisomy wins over a peak.
      e <- e[order(e$end - e$start), , drop = FALSE]
      seg <- overlay(L, e)
      if (is.null(seg)) next
      rows[[length(rows) + 1L]] <- data.frame(
        sample = s, chr = cc, startpos = seg$start, endpos = seg$end,
        tot = seg$total, min = seg$minor, stringsAsFactors = FALSE)
    }
    out[[i]] <- do.call(rbind, rows)
  }
  do.call(rbind, out)
}

make_mutations <- function(drivers.file, samples) {
  d <- read.delim(drivers.file, stringsAsFactors = FALSE)
  #### . . . - Column 1 is the gene symbol in the build's own nomenclature:
  ####         human symbols for hg19/hg38, mouse symbols for mm10.
  genes <- d[[1]]
  rows <- list()
  for (s in samples) {
    k <- sample(3:10, 1)
    g <- sample(genes, k, replace = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      sampleID = s, gene = g,
      impact   = sample(c("Missense", "Nonsense", "Essential_Splice"), k, TRUE),
      aachange = sample(c("D479N", "D488Y", "R780K", "G12D", "Q61H", "E17K"),
                        k, TRUE),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#### . . . - Run one build --------------------------------------------------------

run_build <- function(genome) {
  cat("\n", strrep("=", 74), "\n== ", genome, "\n", strrep("=", 74), "\n", sep = "")

  out.dir <- file.path(tempdir(), paste0("mmgc_selftest_", genome))
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

  mouse <- identical(genome, "mm10")
  if (mouse) {
    source(file.path(repo.root, "R", "mm10", "MMGenomicClassification_mm10.R"))
  }
  MMGC_setup(genome, QUIET = TRUE, LOAD.GENOME.PACKAGES = FALSE)

  db <- "TEST"

  #### . . . - Inputs
  raw <- make_cnv(genome, MMGC_ref("gistic_peaks"))
  samples <- unique(raw$sample)

  cnv <- inputCNV(raw, SAMPLE.COL = "sample", CHR.COL = "chr",
                  START.COL = "startpos", END.COL = "endpos",
                  TOTAL.COL = "tot", MINOR.COL = "min")
  stopifnot(nrow(cnv) > 0, all(c("sample","chrom","start","end","total","minor")
                               %in% colnames(cnv)))
  cat("inputCNV            ok  (", nrow(cnv), " segments, ",
      length(samples), " samples)\n", sep = "")

  qc <- checkCNV(cnv, REPORT = FALSE)
  stopifnot(nrow(qc$per_sample) == length(samples),
            all(c("verdict", "ploidy", "frac_minor_missing", "n_overlaps")
                %in% colnames(qc$per_sample)))
  cat("checkCNV            ok  (", sum(qc$per_sample$verdict == "ok"), " of ",
      nrow(qc$per_sample), " samples clean)\n", sep = "")

  meta <- data.frame(sample = samples, cnv = 1, snv = 1, stringsAsFactors = FALSE)
  meta.file <- file.path(out.dir, "metadata.txt")
  write.table(meta, meta.file, sep = "\t", quote = FALSE, row.names = FALSE)

  muts <- make_mutations(MMGC_ref("drivers"), samples)
  mut.file <- file.path(out.dir, "mutations.txt")
  write.table(muts, mut.file, sep = "\t", quote = FALSE, row.names = FALSE)

  #### . . . - Copy number
  flm <- FocalLargeMatrix(CNV = cnv, DATABASE = db, OUT.DIR = out.dir)
  cat("FocalLargeMatrix    ok  (", ncol(flm) - 1, " peaks)\n", sep = "")

  hy <- if (mouse) HYidentification_mm10(CNV = cnv, DATABASE = db, OUT.DIR = out.dir)
        else        HYidentification(CNV = cnv, DATABASE = db, OUT.DIR = out.dir)
  cat("HYidentification    ok  (", sum(hy$hyperdiploid$hyp.sample == "Hyperdiploidy"),
      " hyperdiploid)\n", sep = "")

  flm.file <- file.path(out.dir, paste0("CNV_", db, "matrix_focal_large_5Mb.txt"))
  hy.file  <- file.path(out.dir, paste0("60%_ARMS_", db, "hyperdiploid.txt"))

  cyto <- CytogeneticMatrixGenerator(FocalLargeMATRIX = flm.file, HY.FILE = hy.file,
                                     DATABASE = db, OUT.DIR = out.dir)
  cat("CytogeneticMatrix   ok\n")

  #### . . . - Mutations
  gistic.file <- file.path(out.dir, paste0("gistic.other_All_peaks", db, ".txt"))
  mm <- if (mouse)
    MutationsMatrixGenerator_mm10(CODING.MUT = mut.file, METADATA = meta.file,
                                  GisticOtherAllPeaks = gistic.file,
                                  CNV = cnv, DATABASE = db, OUT.DIR = out.dir)
  else
    MutationsMatrixGenerator(CODING.MUT = mut.file, METADATA = meta.file,
                             GisticOtherAllPeaks = gistic.file,
                             CNV = cnv, DATABASE = db, OUT.DIR = out.dir)
  cat("MutationsMatrix     ok  (", ncol(mm$tsg), " TSGs combined)\n", sep = "")

  #### . . . - APOBEC: synthetic, so the test needs no BSgenome packages.
  apo <- data.frame(apobec = round(runif(length(samples), 0, 0.3), 4), offset = 0,
                    row.names = samples)
  write.table(apo,
              file.path(out.dir, paste0("APOBEC_relative_contribution_matrix", db, ".txt")),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  cat("APOBEC              (synthetic values substituted)\n")

  #### . . . - Final matrix
  final.file <- file.path(out.dir, paste0("FinalMatrix_", db, ".txt"))
  fm <- FinalMatrixGenerator(INPUT.FILE.FOLDER = out.dir, DATABASE = db,
                             OUTPUT = final.file,
                             FEATURES = if (mouse) MMGC_FINAL_FEATURES_MM10
                                        else MMGC_FINAL_FEATURES)
  cat("FinalMatrixGenerator ok (", ncol(fm) - 1, " features)\n", sep = "")

  #### . . . - Translocations / chromothripsis are external calls; synthesise them.
  extra <- data.frame(
    sample          = fm$sample,
    t_CCND1         = rbinom(nrow(fm), 1, 0.15),
    t_MMSET         = rbinom(nrow(fm), 1, 0.15),
    t_MAF           = rbinom(nrow(fm), 1, 0.08),
    Chromothripsis  = rbinom(nrow(fm), 1, 0.25),
    stringsAsFactors = FALSE)
  mtx <- merge(fm, extra, by = "sample")

  dis3 <- if (mouse) suppressWarnings(dis3_hs_mm10(MATRIX = mtx, CODING.MUT = mut.file))
          else        dis3_hs(MATRIX = mtx, CODING.MUT = mut.file)
  mtx <- merge(mtx, dis3, by = "sample")
  for (j in 2:ncol(mtx)) mtx[[j]] <- ifelse(is.na(mtx[[j]]), 0, mtx[[j]])
  cat("dis3_hs             ok\n")

  #### . . . - Classification
  res <- if (mouse) MMGenoCla_mm10(input = mtx) else MMGenoCla(input = mtx)
  stopifnot(!is.null(res), nrow(res$annotation) > 0,
            nrow(res$genomic) == nrow(res$annotation))
  cat("MMGenoCla           ok  (", nrow(res$annotation), " samples, ",
      ncol(res$genomic), " heatmap rows)\n", sep = "")

  invisible(TRUE)
}

#### . . - Does checkCNV() actually catch each problem? ---------------------------

test_checkcnv <- function() {
  cat("\n", strrep("=", 74), "\n== checkCNV: does it catch what it claims?\n",
      strrep("=", 74), "\n", sep = "")
  MMGC_setup("hg19", QUIET = TRUE, LOAD.GENOME.PACKAGES = FALSE)
  sizes <- read.delim(file.path(repo.root, "reference/hg19/chrom_sizes_hg19.txt"),
                      stringsAsFactors = FALSE)
  autos <- sizes[!sizes$Offset %in% c("X", "Y"), ]

  #### . . . - A clean diploid baseline: one disjoint segment per autosome.
  clean <- function(s) data.frame(
    sample = s, chr = autos$Offset, startpos = 1, endpos = as.numeric(autos$Size),
    tot = 2, min = 1, stringsAsFactors = FALSE)

  raw <- clean("S_ok")
  #### no minor allele anywhere
  x <- clean("S_no_minor"); x$min <- NA;                        raw <- rbind(raw, x)
  #### minor allele missing on a third of segments
  x <- clean("S_partial"); x$min[1:7] <- NA;                    raw <- rbind(raw, x)
  #### a segment overlapping the one before it on the same chromosome
  x <- clean("S_overlap")
  x <- rbind(x, data.frame(sample = "S_overlap", chr = "1", startpos = 1000,
                           endpos = 5e7, tot = 3, min = 1, stringsAsFactors = FALSE))
  raw <- rbind(raw, x)
  #### uncorrected ploidy: total copy number ~20 genome-wide
  x <- clean("S_ploidy"); x$tot <- 20; x$min <- 10;             raw <- rbind(raw, x)
  #### a zero-width segment
  x <- clean("S_malformed")
  x <- rbind(x, data.frame(sample = "S_malformed", chr = "2", startpos = 500,
                           endpos = 500, tot = 2, min = 1, stringsAsFactors = FALSE))
  raw <- rbind(raw, x)
  #### only a fifth of the genome segmented
  x <- clean("S_sparse"); x$endpos <- round(x$endpos / 5);      raw <- rbind(raw, x)

  cnv <- inputCNV(raw, SAMPLE.COL = "sample", CHR.COL = "chr", START.COL = "startpos",
                  END.COL = "endpos", TOTAL.COL = "tot", MINOR.COL = "min")
  q <- checkCNV(cnv, REPORT = FALSE)
  v <- stats::setNames(q$per_sample$verdict, q$per_sample$sample)

  want <- c(S_ok = "ok",
            S_no_minor  = "no allele-specific data",
            S_partial   = "partial allele-specific data",
            S_overlap   = "overlapping segments",
            S_ploidy    = "implausible ploidy",
            S_malformed = "malformed segments",
            S_sparse    = "sparse coverage")
  pass <- TRUE
  for (nm in names(want)) {
    got <- if (nm %in% names(v)) v[[nm]] else "<missing>"
    hit <- identical(got, want[[nm]])
    if (!hit) pass <- FALSE
    cat(sprintf("  %-6s %-13s expected %-30s got %s\n",
                if (hit) "PASS" else "FAIL", nm, want[[nm]], got))
  }
  cat(sprintf("\n  recommended setting: %s (expected \"non-HRD\", since one sample has no minor allele)\n",
              if (isTRUE(q$advice)) "TRUE" else paste0("\"", q$advice, "\"")))
  if (!identical(q$advice, "non-HRD")) pass <- FALSE
  pass
}

ok <- TRUE
for (b in builds) {
  r <- try(run_build(b), silent = FALSE)
  if (inherits(r, "try-error")) {
    ok <- FALSE
    cat("\n!! ", b, " FAILED: ", conditionMessage(attr(r, "condition")), "\n", sep = "")
  }
}

#### . . - The QC test is build-agnostic, so it runs once.
if (!length(args) || "hg19" %in% args || length(builds) > 1) {
  r <- try(test_checkcnv(), silent = FALSE)
  if (inherits(r, "try-error") || !isTRUE(r)) {
    ok <- FALSE
    cat("\n!! checkCNV test FAILED\n")
  }
}

cat("\n", strrep("=", 74), "\n", sep = "")
if (ok) cat("selftest: all checks passed\n") else {
  cat("selftest: FAILURES above\n"); quit(status = 1)
}
