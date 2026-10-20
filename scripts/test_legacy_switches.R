#!/usr/bin/env Rscript
#
# test_legacy_switches.R
#
# Regression test for the defects documented in docs/KNOWN_ISSUES.md, and for
# the REQUIRE.MINOR.ALLELE contract. It
# demonstrates each one on a minimal synthetic case and checks that:
#
#   * the default (corrected) behaviour is right, and
#   * the LEGACY.* switch still reproduces the published behaviour.
#
# Run it after changing anything in HYidentification() or MMGenoCla().
#
# Usage:
#   Rscript scripts/test_legacy_switches.R

this.file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
repo.root <- if (length(this.file)) dirname(dirname(normalizePath(this.file))) else getwd()
setwd(repo.root)

suppressMessages(source("MMGenomicClassification.R"))
MMGC_setup("hg38", QUIET = TRUE, LOAD.GENOME.PACKAGES = FALSE)

out <- file.path(tempdir(), "mmgc_legacy_test")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

pass <- TRUE
check <- function(ok, what) {
  cat(if (ok) "  PASS  " else "  FAIL  ", what, "\n", sep = "")
  if (!ok) pass <<- FALSE
}

################################################################################
cat("\n== Defect 1: gain fractions compared as text, not numbers ==\n\n")
################################################################################

#### . . . - Two samples, both with gains on chromosomes 3, 5 and 7.
####         S_TINY's gains are 8 kb -- about 4e-5 of a chromosome, which is where
####         R switches to scientific notation and the text comparison breaks.
####         S_REAL's gains span whole chromosomes.
sizes <- read.delim("reference/hg38/chrom_sizes_hg38.txt", stringsAsFactors = FALSE)
autos <- sizes[!sizes$Offset %in% c("X", "Y"), ]

make_sample <- function(s, span) {
  backbone <- data.frame(sample = s, chr = autos$Offset, start = 1,
                         end = as.numeric(autos$Size), tot = 2, min = 1,
                         stringsAsFactors = FALSE)
  gains <- do.call(rbind, lapply(c("3", "5", "7"), function(c) {
    L <- as.numeric(autos$Size[autos$Offset == c])
    data.frame(sample = s, chr = c, start = 1e6,
               end = 1e6 + if (span == "tiny") 8000 else L - 2e6,
               tot = 3, min = 1, stringsAsFactors = FALSE)
  }))
  rbind(backbone, gains)
}

cnv <- inputCNV(rbind(make_sample("S_TINY", "tiny"), make_sample("S_REAL", "full")),
                SAMPLE.COL = "sample", CHR.COL = "chr", START.COL = "start",
                END.COL = "end", TOTAL.COL = "tot", MINOR.COL = "min")

fixed  <- suppressMessages(
  HYidentification(CNV = cnv, DATABASE = "FIX", OUT.DIR = out))
legacy <- suppressWarnings(suppressMessages(
  HYidentification(CNV = cnv, DATABASE = "LEG", OUT.DIR = out,
                   LEGACY.STRING.COMPARISON = TRUE)))

call_for <- function(res, s) {
  v <- res$hyperdiploid$hyp.sample[res$hyperdiploid$sample == s]
  if (!length(v)) "not called" else v
}

cat(sprintf("  %-8s  8 kb gains on chr3/5/7    default: %-14s legacy: %s\n",
            "S_TINY", call_for(fixed, "S_TINY"), call_for(legacy, "S_TINY")))
cat(sprintf("  %-8s  whole-chromosome gains    default: %-14s legacy: %s\n\n",
            "S_REAL", call_for(fixed, "S_REAL"), call_for(legacy, "S_REAL")))

check(call_for(fixed, "S_TINY")  != "Hyperdiploidy",
      "8 kb gains are not hyperdiploidy by default")
check(call_for(fixed, "S_REAL")  == "Hyperdiploidy",
      "real trisomies are still called")
check(call_for(legacy, "S_TINY") == "Hyperdiploidy",
      "LEGACY.STRING.COMPARISON reproduces the false positive")

################################################################################
cat("\n== Defect 2: add_gains always zero, HRD_Gains unreachable ==\n\n")
################################################################################

#### . . . - Hyperdiploid samples with four secondary gains and nothing else: the
####         definition of B_HRD_gains / HRD_Gains.
m <- data.frame(sample = sprintf("S%02d", 1:20), stringsAsFactors = FALSE)
for (f in setdiff(MMGC_FINAL_FEATURES, "sample")) m[[f]] <- 0
for (f in MMGC_FEATURES_HRD_CHROM) m[[f]] <- 1
for (f in c("CNV_Amp_17q22", "CNV_Amp_2q24.3", "CNV_Amp_2q32.1",
            "CNV_Amp_4q13.2")) m[[f]] <- 1
m$t_CCND1 <- 0; m$t_MMSET <- 0; m$t_MAF <- 0; m$CNV.Sig <- 0
m$SNV.HS_DIS3 <- 0; m$SNV.noHS_DIS3 <- 0; m$APOBEC <- 0

groups <- function(res) {
  tb <- table(res$annotation$clusters.new)
  paste(names(tb)[tb > 0], collapse = ", ")
}

fx <- suppressWarnings(suppressMessages(MMGenoCla(m, QUIET = TRUE)))
lg <- suppressWarnings(suppressMessages(
  MMGenoCla(m, LEGACY.ADD.GAINS = TRUE, QUIET = TRUE)))

cat("  hyperdiploid + 4 secondary gains   default: ", groups(fx), "\n", sep = "")
cat("                                     legacy:  ", groups(lg), "\n\n", sep = "")

check(any(fx$annotation$clusters.new == "HRD_Gains"),
      "HRD_Gains is reachable by default")
check(!any(lg$annotation$clusters.new == "HRD_Gains"),
      "LEGACY.ADD.GAINS reproduces the empty HRD_Gains group")

################################################################################
cat("\n== REQUIRE.MINOR.ALLELE: the three settings ==\n\n")
################################################################################

#### . . . - Three samples, each with whole-chromosome gains on 3, 5 and 7, but
####         differing in the minor allele: retained (3+1), lost (3+0), or unknown.
mk2 <- function(s, minor) rbind(
  data.frame(sample = s, chr = autos$Offset, start = 1,
             end = as.numeric(autos$Size), tot = 2, min = 1, stringsAsFactors = FALSE),
  data.frame(sample = s, chr = c("3", "5", "7"), start = 1,
             end = as.numeric(autos$Size[match(c("3", "5", "7"), autos$Offset)]),
             tot = 3, min = minor, stringsAsFactors = FALSE))

cnv2 <- inputCNV(rbind(mk2("S_3plus1", 1), mk2("S_3plus0", 0), mk2("S_NA", NA)),
                 SAMPLE.COL = "sample", CHR.COL = "chr", START.COL = "start",
                 END.COL = "end", TOTAL.COL = "tot", MINOR.COL = "min")

called <- function(setting) {
  r <- suppressWarnings(suppressMessages(
    HYidentification(CNV = cnv2, DATABASE = "M", OUT.DIR = out,
                     REQUIRE.MINOR.ALLELE = setting)))
  h <- r$hyperdiploid
  vapply(c("S_3plus1", "S_3plus0", "S_NA"), function(x) {
    v <- h$hyp.sample[h$sample == x]
    if (!length(v)) "no" else if (v == "Hyperdiploidy") "yes" else "no"
  }, character(1))
}
g.true <- called(TRUE); g.non <- called("non-HRD"); g.false <- called(FALSE)

cat(sprintf("  %-10s  3+1: %-4s 3+0: %-4s minor=NA: %s\n", "TRUE",
            g.true[1], g.true[2], g.true[3]))
cat(sprintf("  %-10s  3+1: %-4s 3+0: %-4s minor=NA: %s\n", "non-HRD",
            g.non[1], g.non[2], g.non[3]))
cat(sprintf("  %-10s  3+1: %-4s 3+0: %-4s minor=NA: %s\n\n", "FALSE",
            g.false[1], g.false[2], g.false[3]))

check(identical(unname(g.true), c("yes", "no", "no")),
      "TRUE (default, as published) counts only 3+1")
check(identical(unname(g.non), c("yes", "yes", "yes")),
      "\"non-HRD\" counts 3+0 and unknown-minor trisomies")
check(identical(unname(g.false), c("yes", "yes", "yes")),
      "FALSE counts them too")
check(inherits(try(HYidentification(CNV = cnv2, DATABASE = "M", OUT.DIR = out,
                                    REQUIRE.MINOR.ALLELE = "maybe"), silent = TRUE),
               "try-error"),
      "an invalid REQUIRE.MINOR.ALLELE is rejected")

################################################################################

cat("\n", strrep("=", 70), "\n", sep = "")
if (pass) {
  cat("test_legacy_switches: all checks passed\n")
} else {
  cat("test_legacy_switches: FAILURES above\n"); quit(status = 1)
}
