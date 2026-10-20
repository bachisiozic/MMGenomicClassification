#### . . . - X-linked drivers - . . . ####
#
# These are carried as mutation-only features. They sit on chromosome X, where a
# single-copy loss is not interpretable the way it is on an autosome, so the
# deletion half of the call is dropped and only the mutation is kept.

MMGC_CHRX_GENES <- c("ATRX", "KDM5C", "KDM6A", "RPL10")

#### . . . - How far either side of a gene's transcript to look for an overlapping
####         copy-number segment (bp).
MMGC_GENE_FLANK <- 50000


#### . . . - MutationsMatrixGenerator() - . . . ####

#### . . . - Build the per-gene driver matrices from coding mutations and, for
####         tumour suppressors, from copy number as well.
####
####         Writes six tables into OUT.DIR:
####           oncogene_<DB>.txt           oncogene mutations, 0/1
####           tsg_mutated_genes_<DB>.txt  TSG mutations, 0/1
####           tsg_deleted_genes_<DB>.txt  TSG deletions, 0/1/2
####           tsg_<DB>.txt                the two combined, 0/1/2
####           unk_mutated_genes_<DB>.txt  drivers of unknown role, 0/1
####           chrx_tsg_<DB>.txt           X-linked drivers, mutation only
####
####         A gene is scored 0 when the sample was assessed and is wild type, and
####         NA when the relevant assay is missing for that sample -- which is what
####         METADATA is for.
####
####         CODING.MUT  = non-synonymous coding mutations. Needs columns
####                       sampleID, gene, impact (the dndscv annotmuts layout).
####         METADATA    = per-sample table with columns sample, snv, cnv, each 0/1
####                       saying whether that data type exists for the sample.
####         UCSC.gene   = gene model table (chrom, txStart, txEnd, geneSymbol).
####                       Defaults to the active build's.
####         DRIVERs     = driver summary table. Defaults to the active build's.
####         GisticOtherAllPeaks = gistic.other_All_peaks<DB>.txt from
####                       CytogeneticMatrixGenerator(); used for the peak-vs-gene
####                       cross-check returned in `qc`.
####         CNV         = output of inputCNV()
####         DATABASE    = cohort label
####         OUT.DIR     = output directory
####         CHRX.GENES  = X-linked drivers (default MMGC_CHRX_GENES)
####         UCSC.HG19   = old name for UCSC.gene, still accepted.
####
####         Returns the tables invisibly as a named list.

MutationsMatrixGenerator <- function(CODING.MUT,
                                     METADATA,
                                     UCSC.gene = NULL,
                                     DRIVERs = NULL,
                                     GisticOtherAllPeaks = NULL,
                                     CNV,
                                     DATABASE,
                                     OUT.DIR,
                                     CHRX.GENES = MMGC_CHRX_GENES,
                                     UCSC.HG19 = NULL) {

  #### . . . - UCSC.HG19 was the hg19-only argument name; keep it working.
  if (!is.null(UCSC.HG19) && is.null(UCSC.gene)) UCSC.gene <- UCSC.HG19
  if (is.null(UCSC.gene)) UCSC.gene <- MMGC_ref("gene_models")
  if (is.null(DRIVERs))   DRIVERs   <- MMGC_ref("drivers")

  ################################################################################
  #### . . . - Drivers, split by role - . . . ####
  ################################################################################

  drivers <- read.delim(DRIVERs, stringsAsFactors = FALSE)
  for (col in c("Role", "MM.driver", "DNDS", "MUTSIG", "ONCODRIVER", "FISHHOOK")) {
    if (!col %in% colnames(drivers)) {
      stop("MutationsMatrixGenerator(): ", DRIVERs, " has no column ", sQuote(col),
           call. = FALSE)
    }
  }

  #### . . . - A gene is kept when at least two callers found it, or when it is a
  ####         known myeloma driver.
  drivers$n.CallBy.new <- drivers$DNDS + drivers$MUTSIG +
                          drivers$ONCODRIVER + drivers$FISHHOOK
  drivers <- drivers[which(drivers$n.CallBy.new > 1 | drivers$MM.driver == "MM"), , drop = FALSE]
  colnames(drivers)[1] <- "Gene.Symbol"
  drivers <- drivers[!is.na(drivers$Gene.Symbol), , drop = FALSE]

  #### . . . - hg19-era symbols are rewritten to current names so that one set of
  ####         feature vectors serves both human builds.
  drivers$Gene.Symbol <- MMGC_canonical_symbols(drivers$Gene.Symbol)

  #### . . . - "No role assigned" is spelled NA in the human driver tables and as an
  ####         empty string in the mouse one. The original tested only for NA and
  ####         "Unknown", so every gene with a blank role silently vanished from all
  ####         three matrices -- which is why nine mouse mutation features came out
  ####         missing. Normalise first.
  role <- trimws(as.character(drivers$Role))
  role[is.na(role) | role == "" | role == "NA"] <- "Unknown"

  oncogene <- unique(drivers$Gene.Symbol[role == "Oncogene"])
  tsg      <- unique(drivers$Gene.Symbol[role == "TSG"])
  unknown  <- unique(drivers$Gene.Symbol[role == "Unknown"])

  ################################################################################
  #### . . . - Mutations and per-sample data availability - . . . ####
  ################################################################################

  mut <- read.delim(CODING.MUT, stringsAsFactors = FALSE)
  for (col in c("sampleID", "gene", "impact")) {
    if (!col %in% colnames(mut)) {
      stop("MutationsMatrixGenerator(): ", CODING.MUT, " has no column ",
           sQuote(col), ".\n  Expected the dndscv annotmuts layout ",
           "(sampleID, gene, impact, aachange).", call. = FALSE)
    }
  }
  mut <- mut[which(mut$impact != "Synonymous"), , drop = FALSE]
  mut$gene <- MMGC_canonical_symbols(mut$gene)

  meta <- read.delim(METADATA, stringsAsFactors = FALSE)
  for (col in c("sample", "snv", "cnv")) {
    if (!col %in% colnames(meta)) {
      stop("MutationsMatrixGenerator(): ", METADATA, " has no column ",
           sQuote(col), ".\n  It needs one row per sample with 0/1 flags ",
           "'snv' and 'cnv'.", call. = FALSE)
    }
  }
  snv.assessed <- meta$sample[meta$snv == 1]
  snv.missing  <- meta$sample[meta$snv == 0]
  cnv.missing  <- meta$sample[meta$cnv == 0]

  cnv <- as.data.frame(CNV, stringsAsFactors = FALSE)

  ################################################################################
  #### . . . - Mutation matrices - . . . ####
  ################################################################################

  matrix.oncogene <- .mmgc_mutation_matrix(mut, oncogene, snv.assessed, snv.missing)
  matrix.tsg.mut  <- .mmgc_mutation_matrix(mut, tsg,      snv.assessed, snv.missing)
  matrix.unknown  <- .mmgc_mutation_matrix(mut, unknown,  snv.assessed, snv.missing)

  f.onc <- .mmgc_outfile(OUT.DIR, "oncogene_", DATABASE, ".txt")
  f.tsm <- .mmgc_outfile(OUT.DIR, "tsg_mutated_genes_", DATABASE, ".txt")
  f.unk <- .mmgc_outfile(OUT.DIR, "unk_mutated_genes_", DATABASE, ".txt")
  .mmgc_write(matrix.oncogene, f.onc)
  .mmgc_write(matrix.tsg.mut,  f.tsm)
  .mmgc_write(matrix.unknown,  f.unk)

  ################################################################################
  #### . . . - TSG deletions from copy number - . . . ####
  ################################################################################

  ref.genes <- read.delim(UCSC.gene, stringsAsFactors = FALSE)
  for (col in c("chrom", "txStart", "txEnd", "geneSymbol")) {
    if (!col %in% colnames(ref.genes)) {
      stop("MutationsMatrixGenerator(): ", UCSC.gene, " has no column ",
           sQuote(col), ".\n  Build the slim gene models with ",
           "scripts/build_slim_gene_references.R.", call. = FALSE)
    }
  }
  ref.genes$geneSymbol <- MMGC_canonical_symbols(ref.genes$geneSymbol)

  #### . . . - One transcript per gene: the first listed, as before.
  ts.genes <- ref.genes[ref.genes$geneSymbol %in% tsg,
                        c("chrom", "txStart", "txEnd", "geneSymbol")]
  ts.genes <- ts.genes %>% dplyr::distinct(geneSymbol, .keep_all = TRUE)
  colnames(ts.genes) <- c("chrom", "start", "end", "Symbol")
  ts.genes$chrom <- gsub("^chr", "", ts.genes$chrom)

  unplaced <- setdiff(tsg, ts.genes$Symbol)
  if (length(unplaced)) {
    message("  note: ", length(unplaced), " TSG(s) absent from the gene models, ",
            "so they carry mutation calls only: ",
            paste(unplaced, collapse = ", "))
  }

  matrix.tsg.del <- .mmgc_gene_cn_matrix(cnv, ts.genes, cnv.missing,
                                         FLANK = MMGC_GENE_FLANK)
  colnames(matrix.tsg.del) <- paste0("Del_", colnames(matrix.tsg.del))

  f.tsd <- .mmgc_outfile(OUT.DIR, "tsg_deleted_genes_", DATABASE, ".txt")
  .mmgc_write(matrix.tsg.del, f.tsd)

  ################################################################################
  #### . . . - Mutation and deletion combined, per TSG - . . . ####
  ################################################################################

  #### . . . - 0 wild type, 1 one hit, 2 two hits (biallelic loss, or deletion plus
  ####         mutation). NA only when BOTH assays are missing for the sample; when
  ####         just one is missing, the other one's call stands.
  ####
  ####         The original did this with a set of magic numbers (NA -> 999, then
  ####         1998 -> NA, 1001 -> 2, 1000 -> 1, ...) over hardcoded column
  ####         positions 1:36, 37:72 and 73:108. That silently produced garbage for
  ####         any driver list that did not contain exactly 36 tumour suppressors.
  ####         Here the two matrices are paired by gene symbol, so the number of
  ####         TSGs no longer has to be 36 -- which is what makes the mouse driver
  ####         list usable.
  matrix.tsg <- .mmgc_combine_hits(
    DELETION = matrix.tsg.del,
    MUTATION = matrix.tsg.mut,
    GENES    = sort(tsg))

  f.tsg <- .mmgc_outfile(OUT.DIR, "tsg_", DATABASE, ".txt")
  .mmgc_write(matrix.tsg, f.tsg)

  ################################################################################
  #### . . . - X-linked drivers: mutation only - . . . ####
  ################################################################################

  chrx <- matrix.tsg.mut[, colnames(matrix.tsg.mut) %in% CHRX.GENES, drop = FALSE]
  f.chrx <- .mmgc_outfile(OUT.DIR, "chrx_tsg_", DATABASE, ".txt")
  .mmgc_write(chrx, f.chrx)

  ################################################################################
  #### . . . - QC: GISTIC peak vs the gene it stands for - . . . ####
  ################################################################################

  #### . . . - The original built this table and then never looked at it. It is
  ####         genuinely useful -- it says whether the peak-level and gene-level
  ####         calls agree -- so it is returned rather than discarded.
  qc.peak.vs.gene <- NULL
  if (!is.null(GisticOtherAllPeaks) && file.exists(GisticOtherAllPeaks)) {
    gistic <- read.delim(GisticOtherAllPeaks, stringsAsFactors = FALSE)
    have <- intersect(names(MMGC_PEAK_TSG), colnames(gistic))
    if (length(have)) {
      qc.peak.vs.gene <- gistic[, have, drop = FALSE]
      colnames(qc.peak.vs.gene) <- unname(MMGC_PEAK_TSG[have])
    }
  }

  message("MutationsMatrixGenerator(): ", length(oncogene), " oncogenes, ",
          length(tsg), " TSGs, ", length(unknown), " unknown-role drivers -> ",
          OUT.DIR)

  invisible(list(oncogene        = matrix.oncogene,
                 tsg_mutated     = matrix.tsg.mut,
                 tsg_deleted     = matrix.tsg.del,
                 tsg             = matrix.tsg,
                 unknown         = matrix.unknown,
                 chrx            = chrx,
                 qc_peak_vs_gene = qc.peak.vs.gene))
}


#### . . . - .mmgc_mutation_matrix() - . . . ####
#### . . . - sample x gene 0/1 matrix of non-synonymous mutations.
####         Samples with no mutation calls at all are NA; samples that were
####         sequenced but carry no mutation in a gene are 0. Every gene in GENES
####         gets a column even if nothing in the cohort is mutated there.

.mmgc_mutation_matrix <- function(MUT, GENES, ASSESSED, MISSING) {
  GENES <- sort(unique(GENES))
  samples <- unique(c(ASSESSED, MISSING))

  out <- matrix(0, nrow = length(samples), ncol = length(GENES),
                dimnames = list(samples, GENES))

  hit <- MUT[MUT$gene %in% GENES & MUT$sampleID %in% samples, , drop = FALSE]
  if (nrow(hit)) {
    tab <- table(factor(hit$sampleID, levels = samples),
                 factor(hit$gene, levels = GENES))
    out[] <- ifelse(as.matrix(tab) > 0, 1, 0)
  }

  out[rownames(out) %in% MISSING, ] <- NA
  as.data.frame(out)
}


#### . . . - .mmgc_gene_cn_matrix() - . . . ####
#### . . . - sample x gene copy-number call over each gene's locus:
####           0  retained (total >= 1.5 with minor allele present)
####           1  one copy lost (total < 1.5, or minor allele lost with major kept)
####           2  both copies lost (total < 0.5)
####         Where several segments cover a gene, the worst call wins.
####         Samples without copy-number data are NA.

.mmgc_gene_cn_matrix <- function(CNV, GENES, CNV.MISSING, FLANK = 50000) {

  genes <- GENES
  genes$start.flank <- genes$start
  genes$end.flank   <- genes$end

  cnv <- CNV
  cnv$chrom <- gsub("^chr", "", as.character(cnv$chrom))
  #### . . . - The flank is applied to the SEGMENT, as in the original, so a
  ####         segment ending just short of the gene still counts.
  cnv$start.flank <- cnv$start - FLANK
  cnv$end.flank   <- cnv$end + FLANK

  hits <- .mmgc_overlap(
    A = cnv, B = genes,
    A.CHROM = "chrom", A.START = "start.flank", A.END = "end.flank",
    B.CHROM = "chrom", B.START = "start",       B.END = "end",
    A.KEEP = c("sample", "total", "minor"),
    B.KEEP = "Symbol")

  samples <- unique(c(CNV$sample, CNV.MISSING))
  out <- matrix(0, nrow = length(samples), ncol = nrow(genes),
                dimnames = list(samples, genes$Symbol))

  if (nrow(hits)) {
    #### . . . - Round the copy number onto integer states first, as before.
    total <- as.numeric(hits$total)
    total <- ifelse(total < 2.5 & total >= 1.5, 2,
             ifelse(total < 1.5 & total >= 0.5, 1,
             ifelse(total < 3.5 & total >= 2.5, 3,
             ifelse(total < 4.5 & total >= 3.5, 4, total))))
    minor <- as.numeric(hits$minor)
    minor <- ifelse(minor < 0.5, 0, minor)
    major <- total - minor

    code <- ifelse(major + minor >= 1.5 & minor >= 0.5, 0,
            ifelse(major + minor < 1.5 | (minor < 0.5 & major >= 1.5), 1,
            ifelse(total < 0.5, 2, NA)))
    #### . . . - Fall back on total alone where the allelic split is unusable.
    code <- ifelse(!is.na(code), code,
            ifelse(total < 0.5, 2,
            ifelse(total < 1.5, 1, 0)))

    #### . . . - Worst call wins where several segments cover one gene. A gene
    ####         covered only by segments with an unknown copy number stays NA
    ####         rather than collapsing to -Inf.
    agg <- stats::aggregate(list(code = code),
                            by = list(sample = hits$sample, Symbol = hits$Symbol),
                            FUN = function(v) if (all(is.na(v))) NA_real_ else
                                              max(v, na.rm = TRUE))
    idx <- cbind(match(agg$sample, rownames(out)), match(agg$Symbol, colnames(out)))
    keep <- !is.na(idx[, 1]) & !is.na(idx[, 2])
    out[idx[keep, , drop = FALSE]] <- agg$code[keep]
  }

  out[rownames(out) %in% CNV.MISSING, ] <- NA
  as.data.frame(out)
}


#### . . . - .mmgc_combine_hits() - . . . ####
#### . . . - Pair a deletion matrix (Del_<gene>, 0/1/2) with a mutation matrix
####         (<gene>, 0/1) by GENE SYMBOL and add the hits, capped at 2.
####         NA only where both are NA.

.mmgc_combine_hits <- function(DELETION, MUTATION, GENES) {

  samples <- union(rownames(DELETION), rownames(MUTATION))
  out <- data.frame(row.names = samples)

  for (g in GENES) {
    d.col <- paste0("Del_", g)
    d <- if (d.col %in% colnames(DELETION))
      as.numeric(DELETION[[d.col]])[match(samples, rownames(DELETION))] else
      rep(NA_real_, length(samples))
    m <- if (g %in% colnames(MUTATION))
      as.numeric(MUTATION[[g]])[match(samples, rownames(MUTATION))] else
      rep(NA_real_, length(samples))

    both.missing <- is.na(d) & is.na(m)
    total <- rowSums(cbind(d, m), na.rm = TRUE)
    total[total > 2] <- 2
    total[both.missing] <- NA_real_
    out[[g]] <- total
  }
  out
}


#### . . . - .mmgc_write() - . . . ####

.mmgc_write <- function(X, FILE) {
  write.table(X, FILE, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
  invisible(FILE)
}
