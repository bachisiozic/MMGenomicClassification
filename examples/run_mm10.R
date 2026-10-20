################################################################################
####                                                                        ####
####   MMGenomicClassification -- complete run, mm10 (mouse)                ####
####                                                                        ####
####   READ docs/MOUSE.md FIRST. Hyperdiploidy has no mouse equivalent and   ####
####   is switched off by default; the copy-number peaks and driver genes    ####
####   do transfer, via orthology.                                          ####
####                                                                        ####
################################################################################

################################################################################
#### . . . - CONFIGURE - . . . ####
################################################################################

REPO    <- "~/MMGenomicClassification"
COHORT  <- "vkMYC"
OUT.DIR <- "~/analysis/vkMYC"

CNV.FILE <- "~/data/vkMYC_cnv.txt"
SNV.FILE <- "~/data/vkMYC_snv.txt"

#### . . . - dndscv reference for MOUSE. The `gene` column of the output must hold
####         mouse symbols; the pipeline translates them to human orthologs itself.
DNDS.REFDB <- "~/ref/RefCDS_mouse_GRCm38.rda"

#### . . . - "none"      HRD features set to 0, no HRD group assigned (default).
#### . . . - "ortholog"  use reference/mm10/hrd_chromosome_map.txt. Only after
####                     reading MMGC_mm10_report() -- the best proxy carries as
####                     little as a third of the relevant orthologs.
HRD.PROXY <- "none"

CNV.COLS <- c(sample = "sample", chr = "chr", start = "startpos",
              end = "endpos", total = "tot", minor = "min")
SNV.COLS <- c(sample = "sample", chrom = "chrom", pos = "pos",
              ref = "ref", alt = "alt")

APOBEC.ITER <- 1000

################################################################################
#### . . . - SETUP - . . . ####
################################################################################

source(file.path(REPO, "MMGenomicClassification.R"))
source(file.path(REPO, "R", "mm10", "MMGenomicClassification_mm10.R"))
MMGC_setup("mm10")

library(dndscv)

dir.create(OUT.DIR, recursive = TRUE, showWarnings = FALSE)

#### . . . - What can and cannot be measured in mouse. Worth reading the output.
MMGC_mm10_report()

################################################################################
#### . . . - 0. Metadata - . . . ####
################################################################################

cnv.raw <- read.delim(CNV.FILE, stringsAsFactors = FALSE)
snv.raw <- read.delim(SNV.FILE, stringsAsFactors = FALSE)

metadata <- merge(
  data.frame(sample = unique(cnv.raw[[CNV.COLS[["sample"]]]]), cnv = 1),
  data.frame(sample = unique(snv.raw[[SNV.COLS[["sample"]]]]), snv = 1),
  by = "sample", all = TRUE)
metadata$cnv[is.na(metadata$cnv)] <- 0
metadata$snv[is.na(metadata$snv)] <- 0

METADATA.FILE <- file.path(OUT.DIR, paste0(COHORT, "_metadata.txt"))
write.table(metadata, METADATA.FILE, sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
#### . . . - 1-4. Copy number - . . . ####
################################################################################

#### . . . - Mouse chromosome coding: 1-19, X = 20, Y = 21.
cnv <- inputCNV_mm10(cnv.raw,
                     SAMPLE.COL = CNV.COLS[["sample"]],
                     CHR.COL    = CNV.COLS[["chr"]],
                     START.COL  = CNV.COLS[["start"]],
                     END.COL    = CNV.COLS[["end"]],
                     TOTAL.COL  = CNV.COLS[["total"]],
                     MINOR.COL  = CNV.COLS[["minor"]])

#### . . . - Check the input before running anything. See run_hg38.R for what this
####         catches and why it matters.
qc <- checkCNV(cnv)

#### . . . - The mouse GISTIC peaks keep the human peak codes, relocated to the
####         orthologous mouse locus -- so this step is unchanged.
FocalLargeMatrix_mm10(CNV = cnv, DATABASE = COHORT, OUT.DIR = OUT.DIR)

#### . . . - Scores gain on the MOUSE chromosomes, then applies HRD.PROXY. The
####         mouse-chromosome table is written alongside, and is the measurement
####         actually made -- use it if you want mouse aneuploidy directly.
HYidentification_mm10(CNV = cnv, DATABASE = COHORT, OUT.DIR = OUT.DIR,
                      HRD.PROXY = HRD.PROXY)

CytogeneticMatrixGenerator_mm10(
  FocalLargeMATRIX = file.path(OUT.DIR,
                               paste0("CNV_", COHORT, "matrix_focal_large_5Mb.txt")),
  HY.FILE          = file.path(OUT.DIR,
                               paste0("60%_ARMS_", COHORT, "hyperdiploid.txt")),
  DATABASE = COHORT, OUT.DIR = OUT.DIR)

################################################################################
#### . . . - 5. Mutations - . . . ####
################################################################################

mutations <- snv.raw[, unname(SNV.COLS[c("sample", "chrom", "pos", "ref", "alt")])]
colnames(mutations) <- c("sampleID", "chr", "pos", "ref", "mut")

dnds <- dndscv(mutations, refdb = DNDS.REFDB, cv = NULL)

coding <- dnds$annotmuts[dnds$annotmuts$impact != "Synonymous", ]
CODING.MUT.FILE <- file.path(OUT.DIR, paste0(COHORT, "_dnds_ns_mutations.txt"))
write.table(coding, CODING.MUT.FILE, sep = "\t", quote = FALSE, row.names = FALSE)

#### . . . - Mouse gene symbols are translated to their human orthologs here, so
####         every feature downstream carries its human name.
MutationsMatrixGenerator_mm10(
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

#### . . . - Human SBS signatures fitted to mm10 contexts. The APOBEC bins use
####         human cohort cut-offs, so treat the result as a rank within this
####         cohort -- see docs/MOUSE.md.
APOBECanalyzer_mm10(SNV.OBJECT = snv.raw,
                    SAMPLE = SNV.COLS[["sample"]], CHROM = SNV.COLS[["chrom"]],
                    POS    = SNV.COLS[["pos"]],    REF   = SNV.COLS[["ref"]],
                    ALT    = SNV.COLS[["alt"]],
                    ITER = APOBEC.ITER, DATABASE = COHORT, OUT.DIR = OUT.DIR)

################################################################################
#### . . . - 7. Feature matrix (117 features, not 126) - . . . ####
################################################################################

FINAL.FILE <- file.path(OUT.DIR, paste0("FinalMatrix_", COHORT, ".txt"))
matrix <- FinalMatrixGenerator_mm10(INPUT.FILE.FOLDER = OUT.DIR,
                                    DATABASE = COHORT,
                                    OUTPUT   = FINAL.FILE)

################################################################################
#### . . . - 8. External calls - . . . ####
################################################################################

#### . . . - IGH translocations do not occur in mouse models as they do in human
####         myeloma. Set them to 0 unless your model has an engineered
####         equivalent, in which case read them from a file as run_hg38.R does.
matrix$t_CCND1 <- 0
matrix$t_MMSET <- 0
matrix$t_MAF   <- 0

#### . . . - Chromothripsis, if you have called it. 0 is "not present".
matrix$CNV.Sig <- 0

#### . . . - Hotspot residues are HUMAN DIS3 positions -- confirm the mouse Dis3
####         alignment, and pass HOTSPOTS explicitly if it differs.
dis3 <- dis3_hs_mm10(MATRIX = matrix, CODING.MUT = CODING.MUT.FILE)
matrix <- merge(matrix, dis3, by = "sample")

################################################################################
#### . . . - 9. Classify - . . . ####
################################################################################

res <- MMGenoCla_mm10(input = matrix)

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
message("Remember: with HRD.PROXY = \"", HRD.PROXY, "\", the HRD groups are ",
        if (identical(HRD.PROXY, "none")) "unreachable by design."
        else "based on a weak ortholog proxy -- see docs/MOUSE.md.")
