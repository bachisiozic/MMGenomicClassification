#### . . . - HYidentification() - . . . ####

#### . . . - Call hyperdiploidy from arm- and chromosome-level gain.
####
####         A segment counts towards gain when its collapsed total copy number is
####         >= 2.5 with minor allele >= 0.5. Gained bases are summed per arm and
####         per chromosome and expressed as a fraction of the arm / chromosome
####         length. An arm or chromosome is "gained" above PERCENTAGE. A sample is
####         hyperdiploid when at least two of HRD.CHROMOSOMES are gained.
####
####         CNV             = output of inputCNV()
####         ARMS            = arm coordinate table. Defaults to the active build's.
####         CHROM.LENGHT    = chromosome size table. Defaults to the active build's.
####                           (Spelling kept from the original signature.)
####         DATABASE        = cohort label
####         OUT.DIR         = where to write the two output tables
####         PERCENTAGE      = fraction of an arm / chromosome that must be gained
####                           (default 0.6)
####         HRD.PERCENTAGE  = threshold for the per-chromosome HRD call. The
####                           original hardcoded 0.6 here even when PERCENTAGE was
####                           set to something else; it is exposed so the two can
####                           be kept in step, and defaults to 0.6 so that
####                           PERCENTAGE = 0.6 behaves exactly as before.
####         HRD.CHROMOSOMES = chromosomes whose trisomy defines hyperdiploidy
####                           (default 3, 5, 7, 9, 11, 15, 19, 21)
####         EXTRA.CHROMOSOMES = reported but excluded from the HRD call (default 18)
####         DOUBLE.GENOME.ARMS = above this many gained arms the sample is called
####                           Double_Genomes instead (default 20)
####         REQUIRE.MINOR.ALLELE = whether a gained segment must also retain its
####                           minor allele. TRUE (default) always requires it,
####                           which is the published code path; "non-HRD" requires
####                           it everywhere except the trisomy chromosomes; FALSE
####                           never requires it. See the note below -- this choice
####                           changes trisomy calls.
####         CYTOBAND.FILE, GISTIC.ARMS = accepted and ignored; see note below.
####
####         LEGACY.STRING.COMPARISON -- see the warning below. Leave it FALSE.
####
####         Writes  <PERCENTAGE*100>%_ARMS_<DATABASE>hyperdiploid.txt  and
####                 <PERCENTAGE*100>%_ARMS_<DATABASE>no.hyperdiploid.txt
####         and returns both invisibly as a list.
####
####     NOTE on CYTOBAND.FILE and GISTIC.ARMS: the original required both and used
####     neither -- it read them into local variables that nothing else touched.
####     They stay in the signature so existing call sites keep working.
####
####     NOTE on REQUIRE.MINOR.ALLELE -- this changes trisomy calls.
####
####     The published function selected gained segments with
####         total >= 2.5 & minor >= 0.5
####     i.e. it required the second allele to be retained, so a chromosome at three
####     copies with the minor allele lost (3+0) did not count as a trisomy.
####     Checked against the published feature matrix, that is NOT how the matrix
####     was made: scoring gain on total copy number alone reproduces its
####     CNV_chr*.gain columns for 15,235 of 15,237 sample-chromosome calls
####     (99.99%), while requiring the minor allele reproduces 96.49% and misses 535
####     trisomies the matrix contains -- every one of them either a 3+0 gain or a
####     sample with no allele-specific data. It also matters practically:
####     allele-specific copy number is absent for whole cohorts, and requiring it
####     makes every one of those samples look diploid.
####
####     The default is nevertheless TRUE -- the published code path -- so that the
####     function behaves as published unless asked otherwise. Be aware of what that
####     costs: on the paper's own cohort it reproduces 96.49% of the published
####     trisomy calls rather than 99.99%, and any sample whose copy-number calls
####     carry no minor allele is silently reported as diploid. The function warns
####     when it finds such samples.
####
####     The two alternatives:
####
####       "non-HRD"  total copy number alone decides a trisomy; a gain anywhere
####                  else must still retain its minor allele. This is the setting
####                  that reproduces the published matrix (99.99% of trisomy
####                  calls), and it keeps the whole-genome-doubling call strict,
####                  since doubling means gaining BOTH alleles and a 3+0 gain
####                  should not count towards it.
####
####       FALSE      never requires the minor allele. Produces an identical feature
####                  matrix to "non-HRD" -- verified at 222,592 of 222,592 cells on
####                  the published cohort -- and differs only in that
####                  whole-genome-doubling label, so "non-HRD" is preferable.
####
####     Why "non-HRD" and FALSE tie on the matrix: the only features this function
####     contributes are the hyp.chr* trisomy flags, which exist solely for
####     HRD.CHROMOSOMES and EXTRA.CHROMOSOMES. Gains on any other chromosome are
####     scored but feed only `Cytogenetic` (Trisomie vs Double_Genomes) and
####     `double.chrs` here, neither of which CytogeneticMatrixGenerator() reads.
####
####     Note also that the rest of the pipeline already makes this split without a
####     switch: FocalLargeMatrix() scores peak-level GAINS on total copy number
####     alone (the allele-aware variant is commented out in the published source)
####     while its LOSSES and the gene-level deletions in MutationsMatrixGenerator()
####     are allele-aware. REQUIRE.MINOR.ALLELE therefore governs trisomy calls
####     only.
####
####     WARNING -- LEGACY.STRING.COMPARISON reproduces a bug in the published
####     code. There, the gain fractions were carried in a character matrix, so
####     `seg.gain.prop > PERCENTAGE` and `perc.gain.chr > PERCENTAGE` compared
####     TEXT, not numbers. Lexicographic order agrees with numeric order for
####     values R prints in full ("0.05" < "0.6"), but not for values small enough
####     to be printed in scientific notation: "2.5e-05" > "0.6" is TRUE. So a
####     chromosome carrying only a tiny gain -- under roughly 10-25 kb, which is
####     a very ordinary focal event -- was scored as gained along its whole
####     length, and two such chromosomes were enough to call a sample
####     hyperdiploid. This version compares numbers. Set
####     LEGACY.STRING.COMPARISON = TRUE only to reproduce output from the
####     original scripts; see docs/KNOWN_ISSUES.md.

HYidentification <- function(CNV,
                             ARMS = NULL,
                             CHROM.LENGHT = NULL,
                             DATABASE,
                             OUT.DIR,
                             PERCENTAGE = 0.6,
                             HRD.PERCENTAGE = 0.6,
                             HRD.CHROMOSOMES = c(3, 5, 7, 9, 11, 15, 19, 21),
                             EXTRA.CHROMOSOMES = 18,
                             DOUBLE.GENOME.ARMS = 20,
                             REQUIRE.MINOR.ALLELE = TRUE,
                             LEGACY.STRING.COMPARISON = FALSE,
                             CYTOBAND.FILE = NULL,
                             GISTIC.ARMS = NULL) {

  if (is.null(ARMS))         ARMS         <- MMGC_ref("arms")
  if (is.null(CHROM.LENGHT)) CHROM.LENGHT <- MMGC_ref("chrom_sizes")

  if (!(isTRUE(REQUIRE.MINOR.ALLELE) || isFALSE(REQUIRE.MINOR.ALLELE) ||
        identical(REQUIRE.MINOR.ALLELE, "non-HRD"))) {
    stop("HYidentification(): REQUIRE.MINOR.ALLELE must be TRUE, FALSE or ",
         "\"non-HRD\".", call. = FALSE)
  }

  if (isTRUE(LEGACY.STRING.COMPARISON)) {
    warning("HYidentification(): LEGACY.STRING.COMPARISON = TRUE reproduces the ",
            "published text-comparison bug -- gains below ~1e-4 of a chromosome ",
            "are scored as whole-chromosome gains. See docs/KNOWN_ISSUES.md.",
            call. = FALSE)
  }

  #### . . . - Gained segments, after collapsing consecutive equal-CN runs
  cnv.n <- SegmentCollapsing(CNV, DATABASE = DATABASE)
  all.samples <- unique(cnv.n$sample)

  #### . . . - Which segments count as gained. See REQUIRE.MINOR.ALLELE above.
  ####
  ####         Indexed through which(), because a logical index carrying NA does not
  ####         skip those rows -- it inserts an all-NA row for each one, and an
  ####         unknown minor allele is common in real calls. (The published code hit
  ####         the same thing and swept the phantom rows up afterwards with
  ####         complete.cases(); selecting only the TRUE rows is the same set.)
  gain.ok <- if (identical(REQUIRE.MINOR.ALLELE, "non-HRD")) {
    #### . . . - Total copy number alone on the trisomy chromosomes, allele-aware
    ####         everywhere else. NA propagates correctly: on an HRD chromosome the
    ####         first term is TRUE so an unknown minor allele does not matter; off
    ####         one, FALSE | NA is NA and which() drops the row.
    hrd.set <- c(HRD.CHROMOSOMES, EXTRA.CHROMOSOMES)
    cnv.n$total >= 2.5 & (cnv.n$chrom %in% hrd.set | cnv.n$minor >= 0.5)
  } else if (isTRUE(REQUIRE.MINOR.ALLELE)) {
    cnv.n$total >= 2.5 & cnv.n$minor >= 0.5
  } else {
    cnv.n$total >= 2.5
  }
  gain <- cnv.n[which(!is.na(cnv.n$chrom) & gain.ok), , drop = FALSE]

  #### . . . - Warn when the minor allele is required but not available: those
  ####         samples cannot register a gain at all, and would silently come back
  ####         as non-hyperdiploid rather than as unassessable.
  if (isTRUE(REQUIRE.MINOR.ALLELE)) {
    no.allelic <- vapply(split(cnv.n$minor, cnv.n$sample),
                         function(v) all(is.na(v)), logical(1))
    if (any(no.allelic)) {
      warning("HYidentification(): ", sum(no.allelic), " of ", length(no.allelic),
              " sample(s) have no minor-allele values at all, so with ",
              "REQUIRE.MINOR.ALLELE = TRUE no gain can be detected for them and ",
              "they will all be reported as non-hyperdiploid. Pass ",
              "REQUIRE.MINOR.ALLELE = FALSE for total-copy-number-only calls.",
              call. = FALSE)
    }
  }

  #### . . . - Arm coordinates
  arms <- read.delim(ARMS, stringsAsFactors = FALSE)
  for (col in c("Arm", "Chromosome", "Start_Position", "End_Position")) {
    if (!col %in% colnames(arms)) {
      stop("HYidentification(): ", ARMS, " has no column ", sQuote(col),
           ".\n  It needs: Arm, Chromosome, Start_Position, End_Position.",
           call. = FALSE)
    }
  }
  arms$arm.length <- arms$End_Position - arms$Start_Position

  chrom.size <- .mmgc_chrom_sizes(CHROM.LENGHT)

  #### . . . - Empty-cohort short circuit: nothing gained anywhere.
  if (!nrow(gain)) {
    return(.mmgc_write_hy(hyper = NULL, no.hyper = all.samples,
                          DATABASE = DATABASE, OUT.DIR = OUT.DIR,
                          PERCENTAGE = PERCENTAGE,
                          HRD.CHROMOSOMES = HRD.CHROMOSOMES,
                          EXTRA.CHROMOSOMES = EXTRA.CHROMOSOMES))
  }

  #### . . . - Gained segment x arm
  hits <- .mmgc_overlap(
    A = gain, B = arms,
    A.CHROM = "chrom",      A.START = "start",          A.END = "end",
    B.CHROM = "Chromosome", B.START = "Start_Position", B.END = "End_Position",
    A.KEEP = c("sample", "total", "minor", "code", "Num_Segments", "type"),
    B.KEEP = c("Arm", "arm.length"))

  if (!nrow(hits)) {
    return(.mmgc_write_hy(hyper = NULL, no.hyper = all.samples,
                          DATABASE = DATABASE, OUT.DIR = OUT.DIR,
                          PERCENTAGE = PERCENTAGE,
                          HRD.CHROMOSOMES = HRD.CHROMOSOMES,
                          EXTRA.CHROMOSOMES = EXTRA.CHROMOSOMES))
  }

  #### . . . - How much of the arm each gained segment covers, clipping the
  ####         segment to the arm boundaries.
  hits$overlap <- pmin(hits$endA, hits$endB) - pmax(hits$startA, hits$startB)
  hits$overlap[hits$overlap < 0] <- 0
  hits$chrom <- as.numeric(as.character(hits$seqnames))

  #### . . . - Per sample and arm
  per.arm <- hits %>%
    dplyr::group_by(sample, chrom, Arm, arm.length, type) %>%
    dplyr::summarise(CN.median  = stats::median(as.numeric(total), na.rm = TRUE),
                     CN.max     = max(as.numeric(total), na.rm = TRUE),
                     CN.min     = min(as.numeric(total), na.rm = TRUE),
                     length.seg = sum(overlap, na.rm = TRUE),
                     .groups    = "drop") %>%
    as.data.frame()
  per.arm$seg.gain.prop <- per.arm$length.seg / per.arm$arm.length

  #### . . . - Per sample and chromosome: gained bases over the whole chromosome
  per.chrom <- per.arm %>%
    dplyr::group_by(sample, chrom) %>%
    dplyr::summarise(gain.length = sum(length.seg, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
  per.chrom$Size <- chrom.size$Size[match(per.chrom$chrom, chrom.size$Offset)]
  per.chrom$perc.gain.chr <- per.chrom$gain.length / per.chrom$Size

  per.arm <- merge(per.arm,
                   per.chrom[, c("sample", "chrom", "gain.length", "perc.gain.chr")],
                   by = c("sample", "chrom"), all.x = TRUE)

  #### . . . - Thresholds. `.mmgc_gt` is where the legacy text comparison lives.
  per.arm$arm.gain <- as.integer(.mmgc_gt(per.arm$seg.gain.prop, PERCENTAGE,
                                          LEGACY.STRING.COMPARISON))
  per.arm$chr.gain <- as.integer(.mmgc_gt(per.arm$perc.gain.chr, PERCENTAGE,
                                          LEGACY.STRING.COMPARISON))

  kept <- per.arm[!is.na(per.arm$chr.gain) & per.arm$chr.gain == 1, , drop = FALSE]

  #### . . . - Per-sample call
  hrd.all <- c(HRD.CHROMOSOMES, EXTRA.CHROMOSOMES)
  rows <- list()
  no.hyper <- setdiff(all.samples, unique(kept$sample))

  for (s in unique(kept$sample)) {
    ks <- kept[which(kept$sample == s), , drop = FALSE]

    n.gained.arms  <- nrow(ks)
    on.hrd         <- ks$chrom %in% c(HRD.CHROMOSOMES, EXTRA.CHROMOSOMES)

    if (n.gained.arms > DOUBLE.GENOME.ARMS) {
      cytogenetic <- "Double_Genomes"
      double.chrs <- sum(ks$arm.gain, na.rm = TRUE)
    } else if (any(on.hrd)) {
      cytogenetic <- "Trisomie"
      double.chrs <- sum(ks$arm.gain[on.hrd], na.rm = TRUE)
    } else {
      no.hyper <- c(no.hyper, s)
      next
    }

    #### . . . - Per-chromosome HRD flags, in the fixed order the downstream
    ####         matrix expects: HRD chromosomes first, then the extras.
    hyp <- vapply(hrd.all, function(cc) {
      p <- unique(ks$perc.gain.chr[ks$chrom == cc])
      if (!length(p) || all(is.na(p))) return(NA_integer_)
      as.integer(.mmgc_gt(p[1], HRD.PERCENTAGE, LEGACY.STRING.COMPARISON,
                          OR.EQUAL = TRUE))
    }, integer(1))
    names(hyp) <- paste0("hyp.chr", hrd.all)

    #### . . . - chr18 (EXTRA.CHROMOSOMES) is reported but does not count here.
    hyp.hrd.only <- hyp[paste0("hyp.chr", HRD.CHROMOSOMES)]
    hyp.sample <- if (sum(hyp.hrd.only, na.rm = TRUE) >= 2) "Hyperdiploidy" else "No.HRD"

    rows[[length(rows) + 1L]] <- data.frame(
      sample          = s,
      type            = ks$type[1],
      whole.median.CN = stats::median(as.numeric(ks$CN.median), na.rm = TRUE),
      as.list(hyp),
      hyp.sample      = hyp.sample,
      Cytogenetic     = cytogenetic,
      double.chrs     = double.chrs,
      stringsAsFactors = FALSE, check.names = FALSE)
  }

  hyper <- if (length(rows)) do.call(rbind, rows) else NULL

  .mmgc_write_hy(hyper = hyper, no.hyper = unique(no.hyper),
                 DATABASE = DATABASE, OUT.DIR = OUT.DIR, PERCENTAGE = PERCENTAGE,
                 HRD.CHROMOSOMES = HRD.CHROMOSOMES,
                 EXTRA.CHROMOSOMES = EXTRA.CHROMOSOMES)
}


#### . . . - .mmgc_gt() - . . . ####
#### . . . - Threshold comparison. LEGACY = TRUE reproduces the published
####         character comparison (see the warning on HYidentification()).

.mmgc_gt <- function(X, THRESHOLD, LEGACY = FALSE, OR.EQUAL = FALSE) {
  if (isTRUE(LEGACY)) {
    x <- as.character(X)
    t <- as.character(THRESHOLD)
    return(if (OR.EQUAL) x >= t else x > t)
  }
  x <- as.numeric(X)
  if (OR.EQUAL) x >= THRESHOLD else x > THRESHOLD
}


#### . . . - .mmgc_write_hy() - . . . ####
#### . . . - Write the hyperdiploid / non-hyperdiploid tables in the layout the
####         rest of the pipeline reads, including the empty-cohort case (which
####         the original did not handle -- it indexed into NULL and failed).

.mmgc_write_hy <- function(hyper, no.hyper, DATABASE, OUT.DIR, PERCENTAGE,
                           HRD.CHROMOSOMES, EXTRA.CHROMOSOMES) {

  hyp.cols <- paste0("hyp.chr", c(HRD.CHROMOSOMES, EXTRA.CHROMOSOMES))

  if (is.null(hyper) || !nrow(hyper)) {
    hyper <- data.frame(sample = character(0), type = character(0),
                        whole.median.CN = numeric(0), stringsAsFactors = FALSE)
    for (cc in hyp.cols) hyper[[cc]] <- integer(0)
    hyper$hyp.sample  <- character(0)
    hyper$Cytogenetic <- character(0)
    hyper$double.chrs <- integer(0)
  }

  no.hyper.df <- data.frame(sample = no.hyper, stringsAsFactors = FALSE)
  if (nrow(no.hyper.df)) no.hyper.df$Cytogenetic <- "No_HRD"
  else no.hyper.df$Cytogenetic <- character(0)

  tag <- paste0(PERCENTAGE * 100, "%_ARMS_", DATABASE)
  f1 <- .mmgc_outfile(OUT.DIR, tag, "hyperdiploid.txt")
  f2 <- .mmgc_outfile(OUT.DIR, tag, "no.hyperdiploid.txt")

  write.table(hyper, f1, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(no.hyper.df, f2, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  n.hrd <- sum(hyper$hyp.sample == "Hyperdiploidy")
  message("HYidentification(): ", nrow(hyper), " sample(s) with gained chromosomes (",
          n.hrd, " hyperdiploid), ", nrow(no.hyper.df), " without -> ", f1)

  invisible(list(hyperdiploid = hyper, no.hyperdiploid = no.hyper.df))
}
