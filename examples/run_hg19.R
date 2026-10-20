################################################################################
####                                                                        ####
####   MMGenomicClassification -- complete run, hg19                        ####
####                                                                        ####
####   Copy this file, change the CONFIGURE block, run it.                  ####
####   For hg38, change GENOME to "hg38" and use an hg38 dndscv reference.   ####
####   For mouse, start from run_mm10.R instead.                            ####
####                                                                        ####
################################################################################

################################################################################
#### . . . - CONFIGURE - . . . ####
################################################################################

REPO    <- "~/MMGenomicClassification"   # where you cloned this repository
GENOME  <- "hg19"
COHORT  <- "MyCohort"                    # label used in every output file name
OUT.DIR <- "~/analysis/MyCohort"         # working folder; created if absent

#### . . . - Your data
CNV.FILE <- "~/data/MyCohort_cnv.txt"    # allele-specific copy number
SNV.FILE <- "~/data/MyCohort_snv.txt"    # single-base substitutions

#### . . . - Your translocation and chromothripsis calls. See docs/INPUT_FORMATS.md.
####         sample, t_CCND1, t_MMSET, t_MAF   /   sample, CNV.Sig
TRANSLOCATION.FILE  <- "~/data/MyCohort_translocations.txt"
CHROMOTHRIPSIS.FILE <- "~/data/MyCohort_chromothripsis.txt"

#### . . . - dndscv ships with the hg19 reference built in, so "hg19" works here.
####         For a specific RefCDS build, give the .rda path instead.
DNDS.REFDB <- "hg19"

#### . . . - Column names in YOUR copy-number file
CNV.COLS <- c(sample = "sample", chr = "chr", start = "startpos",
              end = "endpos", total = "tot", minor = "min")

#### . . . - Column names in YOUR SNV file
SNV.COLS <- c(sample = "sample", chrom = "chrom", pos = "pos",
              ref = "ref", alt = "alt")

#### . . . - Bootstrap iterations for signature fitting. The published analysis
####         used 1000. Lower it only for a quick check.
APOBEC.ITER <- 1000

################################################################################
#### . . . - SETUP - . . . ####
################################################################################

source(file.path(REPO, "MMGenomicClassification.R"))
MMGC_setup(GENOME)

library(dndscv)

dir.create(OUT.DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
#### . . . - 0. Metadata: which assays exist per sample - . . . ####
################################################################################

#### . . . - This is what tells the pipeline apart "wild type" from "not
####         assessed". Without it, missing data reads as absence of a lesion.

cnv.raw <- read.delim(CNV.FILE, stringsAsFactors = FALSE)
snv.raw <- read.delim(SNV.FILE, stringsAsFactors = FALSE)

metadata <- merge(
  data.frame(sample = unique(cnv.raw[[CNV.COLS[["sample"]]]]), cnv = 1),
  data.frame(sample = unique(snv.raw[[SNV.COLS[["sample"]]]]), snv = 1),
  by = "sample", all = TRUE)
metadata$cnv[is.na(metadata$cnv)] <- 0
metadata$snv[is.na(metadata$snv)] <- 0

#### . . . - Samples missing either assay cannot be classified. Check now rather
####         than wondering later why they vanished.
incomplete <- metadata$sample[metadata$cnv == 0 | metadata$snv == 0]
if (length(incomplete)) {
  message(length(incomplete), " sample(s) lack copy number or mutations and will ",
          "not be classified: ", paste(utils::head(incomplete, 10), collapse = ", "))
}

METADATA.FILE <- file.path(OUT.DIR, paste0(COHORT, "_metadata.txt"))
write.table(metadata, METADATA.FILE, sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
#### . . . - 1. Copy number in - . . . ####
################################################################################

cnv <- inputCNV(cnv.raw,
                SAMPLE.COL = CNV.COLS[["sample"]],
                CHR.COL    = CNV.COLS[["chr"]],
                START.COL  = CNV.COLS[["start"]],
                END.COL    = CNV.COLS[["end"]],
                TOTAL.COL  = CNV.COLS[["total"]],
                MINOR.COL  = CNV.COLS[["minor"]])

################################################################################
#### . . . - 1b. Check the copy-number input BEFORE running anything - . . . ####
################################################################################

#### . . . - This pipeline expects allele-specific, purity- and ploidy-corrected
####         WGS segments (ASCAT, Battenberg, FACETS). checkCNV() says whether the
####         input meets that, and which samples do not. The two failures it exists
####         to catch are silent: segments with no minor allele can never register a
####         gain, so those samples come back non-hyperdiploid rather than flagged;
####         and overlapping segments double-count gained bases, inflating the
####         gained fraction. Read the report before going further.
qc <- checkCNV(cnv)

if (nrow(qc$problems)) {
  write.table(qc$per_sample, file.path(OUT.DIR, paste0(COHORT, "_cnv_qc.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Per-sample QC written to ", COHORT, "_cnv_qc.txt")
}

#### . . . - To exclude the samples the check flagged, uncomment:
# cnv <- cnv[cnv$sample %in% qc$per_sample$sample[qc$per_sample$verdict == "ok"], ]

################################################################################
#### . . . - 2. Recurrent peaks, focal vs large - . . . ####
################################################################################

#### . . . - Reference files come from MMGC_setup(); pass GISTIC.FILE only to
####         override.
FocalLargeMatrix(CNV = cnv, DATABASE = COHORT, OUT.DIR = OUT.DIR,
                 FOCAL.LENGTH = 5000000)

################################################################################
#### . . . - 3. Hyperdiploidy - . . . ####
################################################################################

HYidentification(CNV = cnv, DATABASE = COHORT, OUT.DIR = OUT.DIR,
                 PERCENTAGE = 0.6)

################################################################################
#### . . . - 4. Cytogenetic feature tables - . . . ####
################################################################################

CytogeneticMatrixGenerator(
  FocalLargeMATRIX = file.path(OUT.DIR,
                               paste0("CNV_", COHORT, "matrix_focal_large_5Mb.txt")),
  HY.FILE          = file.path(OUT.DIR,
                               paste0("60%_ARMS_", COHORT, "hyperdiploid.txt")),
  DATABASE = COHORT, OUT.DIR = OUT.DIR)

################################################################################
#### . . . - 5. Coding mutations - . . . ####
################################################################################

mutations <- snv.raw[, unname(SNV.COLS[c("sample", "chrom", "pos", "ref", "alt")])]
colnames(mutations) <- c("sampleID", "chr", "pos", "ref", "mut")

dnds <- dndscv(mutations, refdb = DNDS.REFDB, cv = NULL)

#### . . . - Genes significantly selected in THIS cohort, for your own reference.
####         The pipeline uses the published driver list, not this one.
message("Significant drivers in this cohort: ",
        paste(dnds$sel_cv$gene_name[dnds$sel_cv$qglobal_cv < 0.1], collapse = ", "))

coding <- dnds$annotmuts[dnds$annotmuts$impact != "Synonymous", ]
CODING.MUT.FILE <- file.path(OUT.DIR, paste0(COHORT, "_dnds_ns_mutations.txt"))
write.table(coding, CODING.MUT.FILE, sep = "\t", quote = FALSE, row.names = FALSE)

MutationsMatrixGenerator(
  CODING.MUT          = CODING.MUT.FILE,
  METADATA            = METADATA.FILE,
  GisticOtherAllPeaks = file.path(OUT.DIR,
                                  paste0("gistic.other_All_peaks", COHORT, ".txt")),
  CNV      = cnv,
  DATABASE = COHORT,
  OUT.DIR  = OUT.DIR)

################################################################################
#### . . . - 6. APOBEC - . . . ####
################################################################################

#### . . . - The slow step: scales with APOBEC.ITER x samples.
APOBECanalyzer(SNV.OBJECT = snv.raw,
               SAMPLE = SNV.COLS[["sample"]], CHROM = SNV.COLS[["chrom"]],
               POS    = SNV.COLS[["pos"]],    REF   = SNV.COLS[["ref"]],
               ALT    = SNV.COLS[["alt"]],
               ITER = APOBEC.ITER, DATABASE = COHORT, OUT.DIR = OUT.DIR)

################################################################################
#### . . . - 7. Feature matrix - . . . ####
################################################################################

FINAL.FILE <- file.path(OUT.DIR, paste0("FinalMatrix_", COHORT, ".txt"))
matrix <- FinalMatrixGenerator(INPUT.FILE.FOLDER = OUT.DIR,
                               DATABASE = COHORT,
                               OUTPUT   = FINAL.FILE)

################################################################################
#### . . . - 8. Merge the external calls - . . . ####
################################################################################

#### . . . - Translocations and chromothripsis are not computed here.
translocations <- read.delim(TRANSLOCATION.FILE, stringsAsFactors = FALSE)
chromothripsis <- read.delim(CHROMOTHRIPSIS.FILE, stringsAsFactors = FALSE)

matrix <- merge(matrix, translocations, by = "sample", all.x = TRUE)
matrix <- merge(matrix, chromothripsis, by = "sample", all.x = TRUE)

#### . . . - DIS3 hotspots
dis3 <- dis3_hs(MATRIX = matrix, CODING.MUT = CODING.MUT.FILE)
matrix <- merge(matrix, dis3, by = "sample")

#### . . . - A sample absent from the translocation or chromothripsis file does
####         not carry that lesion, so 0. Do NOT blanket-fill everything else:
####         an NA in an assay-driven feature means "not assessed" and should stay
####         NA so the sample is excluded rather than silently called wild type.
for (col in c("t_CCND1", "t_MMSET", "t_MAF", "CNV.Sig", "Chromothripsis")) {
  if (col %in% colnames(matrix)) matrix[[col]][is.na(matrix[[col]])] <- 0
}

################################################################################
#### . . . - 9. Classify - . . . ####
################################################################################

res <- MMGenoCla(input = matrix)

write.table(res$annotation,
            file.path(OUT.DIR, paste0(COHORT, "_GenomicClassification.txt")),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(res$genomic,
            file.path(OUT.DIR, paste0(COHORT, "_matrix.txt")),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

################################################################################
#### . . . - 10. Heatmap - . . . ####
################################################################################

library(pheatmap)

#### . . . - Order samples by group, then put a gap between groups.
ord  <- order(res$annotation$clusters.new)
mtx  <- t(res$genomic[ord, , drop = FALSE])
ann  <- res$annotation[ord, , drop = FALSE]
gaps <- cumsum(table(droplevels(ann$clusters.new)))

pdf(file.path(OUT.DIR, paste0(COHORT, "_heatmap.pdf")), width = 14, height = 10)
pheatmap(mtx,
         annotation_col    = ann,
         annotation_colors = res$colors,
         cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, fontsize_row = 5,
         border_color = NA, legend = FALSE,
         col = c("grey85", "salmon1", "darkred"),
         gaps_col = gaps[-length(gaps)])
dev.off()

message("Done. Results in ", OUT.DIR)
