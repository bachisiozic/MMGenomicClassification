#### . . . - APOBEC / mutational signature fitting - . . . ####
#
# The APOBEC feature is the combined relative contribution of SBS2 and SBS13,
# fitted with mmsig.
#
# mmsig hardcodes BSgenome.Hsapiens.UCSC.hg19 in two places: the trinucleotide
# context lookup in mm_fit_signatures(), and the transcript annotation in its
# internal getStrandBias(). Fitting hg38 or mm10 data against hg19 contexts is
# wrong -- the reference base at a given coordinate differs between builds.
#
# The original hg38 script worked around this by source()-ing a privately patched
# copy of mmsig from an absolute path on one machine:
#
#     source("~/Google Drive/.../mmsig-master/R/mmsig_functions_hg38.R")
#
# which of course fails for everyone else, and leaves no record of what was
# patched. Instead, the two genome-dependent steps are reimplemented here on top
# of mmsig's own unexported internals, taking the genome as an argument. The
# fitting itself -- vcfCleaner(), fit_signatures(), bootstrap_fit_signatures() --
# is mmsig's, untouched.

#### . . . - Signatures fitted by default, in the reference's column order.
MMGC_SIGNATURES <- c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS18")

#### . . . - APOBEC cut-offs: the 6th and 10th deciles of the combined SBS2+SBS13
####         relative contribution across the 1,933-sample cohort the
####         classification was built on. Hardcoded so that a new cohort is scored
####         on the published scale rather than on its own quantiles.
MMGC_APOBEC_CUTOFFS <- c(low = 0.02891755, high = 0.11203844)


#### . . . - APOBECanalyzer() - . . . ####

#### . . . - Fit mutational signatures and write the APOBEC contribution per sample.
####
####         SNV.OBJECT = data frame of single-base substitutions
####         SAMPLE, CHROM, POS, REF, ALT = the column in SNV.OBJECT holding each field
####         ITER       = bootstrap iterations (default 10; the published analysis
####                      used 1000 for the cohort-level fit)
####         DATABASE   = cohort label
####         OUT.DIR    = output directory
####         GENOME     = build to take trinucleotide contexts from. Defaults to
####                      the one set by MMGC_setup().
####         STRAND.BIAS = also test transcriptional strand bias (needs the TxDb
####                      package for the build)
####         PLOT       = draw mmsig's bootstrap plot, as the original did
####
####         Writes APOBEC_relative_contribution_matrix<DATABASE>.txt and returns
####         the full fit invisibly, so the per-signature estimates are available
####         and not just APOBEC.

APOBECanalyzer <- function(SNV.OBJECT,
                           SAMPLE, CHROM, POS, REF, ALT,
                           ITER = 10,
                           DATABASE,
                           OUT.DIR,
                           GENOME = NULL,
                           STRAND.BIAS = TRUE,
                           PLOT = TRUE) {

  .mmgc_require(MMGC_PACKAGES_SIGNATURES, "APOBECanalyzer()")

  genome <- if (is.null(GENOME)) MMGC_genome() else GENOME
  set.seed(1)

  ################################################################################
  #### . . . - Input - . . . ####
  ################################################################################

  snv <- as.data.frame(SNV.OBJECT, stringsAsFactors = FALSE)
  wanted <- c(sample = SAMPLE, chr = CHROM, pos = POS, ref = REF, alt = ALT)
  absent <- wanted[!wanted %in% colnames(snv)]
  if (length(absent)) {
    stop("APOBECanalyzer(): column(s) not found in SNV.OBJECT: ",
         paste(sprintf("%s = %s", names(absent), sQuote(unname(absent))),
               collapse = ", "), call. = FALSE)
  }
  snv <- snv[, unname(wanted), drop = FALSE]
  colnames(snv) <- names(wanted)

  #### . . . - Single-base substitutions only.
  snv <- snv[(nchar(snv$ref) + nchar(snv$alt)) == 2, , drop = FALSE]
  if (!nrow(snv)) {
    stop("APOBECanalyzer(): no single-base substitutions left in SNV.OBJECT.",
         call. = FALSE)
  }
  snv$chr <- paste0("chr", gsub("^chr", "", as.character(snv$chr)))

  ################################################################################
  #### . . . - Signature reference - . . . ####
  ################################################################################

  sig.ref <- mmsig::signature_ref[, c("class", "sub", "tri", MMGC_SIGNATURES)]
  sig <- sig.ref[, c("sub", "tri", "SBS1", "SBS2", "SBS13", "SBS5",
                     "SBS8", "SBS9", "SBS18")]

  ################################################################################
  #### . . . - Fit - . . . ####
  ################################################################################

  fit <- .mmgc_fit_signatures(muts.input = snv,
                              sig.input  = sig,
                              genome     = genome,
                              iterations = ITER,
                              strandbias = STRAND.BIAS)

  if (isTRUE(PLOT)) {
    try(print(mmsig::bootSigsPlot(fit$bootstrap)), silent = TRUE)
  }

  ################################################################################
  #### . . . - APOBEC = SBS2 + SBS13 - . . . ####
  ################################################################################

  est <- fit$estimate
  apobec <- data.frame(
    sample = rownames(est),
    apobec = rowSums(est[, intersect(c("SBS2", "SBS13"), colnames(est)), drop = FALSE],
                     na.rm = TRUE),
    stringsAsFactors = FALSE)

  #### . . . - `offset` is a filler column: the original wrote it so that the file
  ####         had two columns and FinalMatrixGenerator() could drop the last one.
  ####         Kept so the file layout is unchanged.
  out <- data.frame(apobec = apobec$apobec, offset = 0,
                    row.names = apobec$sample)

  f <- .mmgc_outfile(OUT.DIR, "APOBEC_relative_contribution_matrix", DATABASE, ".txt")
  write.table(out, f, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

  message("APOBECanalyzer(): ", nrow(out), " samples, ",
          sum(apobec$apobec > MMGC_APOBEC_CUTOFFS[["low"]]),
          " above the APOBEC threshold -> ", f)

  invisible(list(apobec = out, fit = fit))
}


#### . . . - .mmgc_fit_signatures() - . . . ####
#### . . . - mmsig::mm_fit_signatures(), with the reference genome as a parameter.
####         Only the two genome-dependent steps are reimplemented; the fitting
####         itself calls mmsig's own internals.

.mmgc_fit_signatures <- function(muts.input, sig.input, genome,
                                 iterations = 10, strandbias = TRUE,
                                 cos_sim_threshold = 0.01,
                                 force_include = c("SBS1", "SBS5")) {

  pkgs <- MMGC_PACKAGES_GENOME[[genome]]
  if (is.null(pkgs)) stop(".mmgc_fit_signatures(): unknown genome ", genome, call. = FALSE)
  .mmgc_require(unname(pkgs["bsgenome"]), paste0("Signature fitting on ", genome))
  bsg <- get(pkgs[["bsgenome"]], envir = as.environment(paste0("package:", pkgs[["bsgenome"]])))

  options(scipen = 999)

  #### . . . - Clean the variant table with mmsig's own checker.
  muts <- mmsig:::vcfCleaner(muts.input)
  if (!nrow(muts)) {
    stop(".mmgc_fit_signatures(): mmsig rejected every variant. Check chromosome ",
         "naming and that ref/alt are single bases.", call. = FALSE)
  }

  #### . . . - Trinucleotide contexts, from the requested build. This is the step
  ####         mmsig hardcodes to hg19.
  samples.muts <- deconstructSigs::mut.to.sigs.input(
    mut.ref = muts, sample.id = "sample", chr = "chr", pos = "pos",
    ref = "ref", alt = "alt", bsg = bsg)
  samples.muts <- as.data.frame(t(samples.muts))

  #### . . . - Reference signature definitions, reshaped as mmsig expects.
  consigts.defn <- sig.input
  mutlist <- paste(consigts.defn[, 2],
                   paste(substr(consigts.defn[, 2], 1, 1),
                         substr(consigts.defn[, 1], 3, 3),
                         substr(consigts.defn[, 2], 3, 3), sep = ""),
                   sep = ">")
  consigts.defn <- sapply(consigts.defn[, 3:ncol(consigts.defn)], as.numeric)
  rownames(consigts.defn) <- mutlist

  samples <- colnames(samples.muts)
  sigt.profs <- stats::setNames(
    replicate(length(samples), colnames(consigts.defn), simplify = FALSE), samples)

  output <- list()
  output$estimate <- mmsig:::fit_signatures(
    samples.muts = samples.muts, consigts.defn = consigts.defn,
    sigt.profs = sigt.profs, cos_sim_threshold = cos_sim_threshold,
    force_include = force_include, dbg = FALSE)

  #### . . . - Bootstrap, so downstream code can see the confidence intervals.
  sig.est <- output$estimate
  sig.est$sample <- rownames(sig.est)
  sig.est <- reshape2::melt(sig.est[names(sig.est) != "mutations"],
                            id.vars = "sample", variable.name = "signature",
                            value.name = "estimate")
  sig.est$signature <- as.character(sig.est$signature)

  sigboot <- mmsig:::bootstrap_fit_signatures(
    samples.muts = samples.muts, consigts.defn = consigts.defn,
    sigt.profs = sigt.profs, iterations = iterations,
    cos_sim_threshold = cos_sim_threshold, force_include = force_include)
  output$bootstrap <- dplyr::left_join(sigboot, sig.est,
                                       by = c("sample", "signature"))

  if (isTRUE(strandbias)) {
    sb <- try(.mmgc_strand_bias(muts, genome), silent = TRUE)
    if (inherits(sb, "try-error")) {
      warning("Strand-bias test skipped: ", conditionMessage(attr(sb, "condition")),
              call. = FALSE)
    } else {
      output$strand_bias_all_3nt <- sb$all_3nt
      output$strand_bias_mm1     <- sb$mm1
      output$strand_bias_SBS35   <- sb$SBS35
    }
  }

  output$mutmatrix <- samples.muts
  output
}


#### . . . - .mmgc_strand_bias() - . . . ####
#### . . . - mmsig's transcriptional strand-bias test, with the genome and
####         transcript annotation taken from the active build instead of hg19.

.mmgc_strand_bias <- function(MUTS, GENOME) {

  pkgs <- MMGC_PACKAGES_GENOME[[GENOME]]
  .mmgc_require(unname(pkgs), paste0("Strand-bias testing on ", GENOME))
  txdb <- get(pkgs[["txdb"]], envir = as.environment(paste0("package:", pkgs[["txdb"]])))

  samples <- unique(MUTS$sample)
  vcf.list <- lapply(samples, function(s) {
    d <- MUTS[which(MUTS$sample == s), , drop = FALSE]
    gr <- GenomicRanges::GRanges(d$chr,
                                 IRanges::IRanges(start = d$pos, end = d$pos),
                                 REF = d$ref, ALT = d$alt)
    GenomeInfoDb::genome(gr) <- GENOME
    gr
  })
  names(vcf.list) <- samples

  genes <- suppressMessages(GenomicFeatures::genes(txdb))

  mm <- MutationalPatterns::mut_matrix_stranded(
    vcf_list = vcf.list, ref_genome = pkgs[["bsgenome"]],
    ranges = genes, mode = "transcription")

  stranded <- as.data.frame(mm)
  stranded <- tibble::rownames_to_column(stranded, var = "type")
  stranded <- reshape2::melt(stranded, id.var = "type",
                             variable.name = "group", value.name = "no_mutations")
  stranded <- tidyr::separate(stranded, col = "type",
                              into = c("type", "strand"), sep = "-")

  all.3nt <- suppressMessages(MutationalPatterns::strand_bias_test(stranded))

  #### . . . - Two signature-specific contexts mmsig reports separately:
  ####         SBS-MM1 (melphalan) and SBS35 (platinum).
  ratio.test <- function(contexts, flag.name) {
    sub <- stranded[stranded$type %in% contexts, , drop = FALSE]
    agg <- stats::aggregate(no_mutations ~ group + strand, data = sub, FUN = sum)
    wide <- reshape2::dcast(agg, group ~ strand, value.var = "no_mutations")
    wide$ratio <- wide$transcribed / wide$untranscribed
    wide$p_poisson <- vapply(seq_len(nrow(wide)), function(i) {
      p <- try(stats::poisson.test(c(wide$transcribed[i], wide$untranscribed[i]),
                                   r = 1)$p.value, silent = TRUE)
      if (inherits(p, "try-error") || !is.numeric(p)) NA_real_ else p
    }, numeric(1))
    wide[[flag.name]] <- ifelse(!is.na(wide$p_poisson) & wide$p_poisson < 0.05 &
                                  wide$ratio > 1, "*", "")
    wide
  }

  list(all_3nt = all.3nt,
       mm1   = ratio.test(c("C[C>T]A", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T"),
                          "MM1_flag"),
       SBS35 = ratio.test(c("C[C>A]C", "C[C>T]C"), "SBS35_flag"))
}
