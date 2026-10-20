#!/usr/bin/env Rscript
#
# build_mm10_arms_and_hrd_map.R
#
# Produces the two mm10 reference tables that the mouse pipeline needs and that
# the existing mm10 reference set did not include:
#
#   reference/mm10/arm_mm10.txt             -- arm coordinates (see note 1)
#   reference/mm10/hrd_chromosome_map.txt   -- human HRD chromosome -> mouse
#                                              chromosome, by gene orthology
#                                              (see note 2)
#
# NOTE 1 -- arms.
# Mouse chromosomes are telocentric: the centromere sits at the very start and
# there is effectively no p arm. mm10 cytoBand names reflect this (every band is
# "qA1", "qB2", ... -- there are no p bands and no acen rows). Each mouse
# chromosome therefore gets a single q arm spanning the whole chromosome. The
# table keeps the column names used by the human arm files so the same code path
# reads both.
#
# NOTE 2 -- the HRD chromosome map. READ THIS BEFORE USING THE MOUSE PIPELINE.
# Hyperdiploid myeloma is defined by trisomies of human chromosomes 3, 5, 7, 9,
# 11, 15, 19 and 21 (plus 18, tracked separately). There is no mouse equivalent
# of those trisomies: synteny is broken, and each human chromosome's genes are
# scattered over several mouse chromosomes. This script therefore reports, for
# each human HRD chromosome, which mouse chromosome carries the largest share of
# its orthologs, together with that share, so the mapping can be judged rather
# than assumed. `ortholog_fraction` is usually well below 1 -- that number is the
# honest measure of how poor the analogy is.
#
# The resulting file is EDITABLE, in the same spirit as peak_target_override.txt.
# Set `use` to 0 for any human HRD chromosome that should not be called at all in
# mouse, or change `mouse_chrom` if a different chromosome is the better proxy
# for your model.
#
# Usage:
#   Rscript scripts/build_mm10_arms_and_hrd_map.R
#
# Needs: data.table, babelgene, and the slim gene models built by
# scripts/build_slim_gene_references.R.

suppressPackageStartupMessages({
  library(data.table)
  library(babelgene)
})

this.file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
repo.root <- if (length(this.file)) dirname(dirname(normalizePath(this.file))) else getwd()

ref <- function(...) file.path(repo.root, "reference", ...)

## ---------------------------------------------------------------------------
## 1. arm_mm10.txt
## ---------------------------------------------------------------------------

sizes <- fread(ref("mm10", "chrom_sizes_mm10.txt"), sep = "\t", header = TRUE,
               colClasses = "character")

# The mouse chromosome coding used throughout the mm10 pipeline: 1-19 autosomes,
# X = 20, Y = 21. (Human uses 1-22, X = 23, Y = 24.)
mouse.code <- function(x) {
  x <- gsub("^chr", "", x)
  ifelse(x == "X", "20", ifelse(x == "Y", "21", x))
}

arms <- data.table(
  Arm            = paste0(mouse.code(sizes$Offset), "q"),
  Chromosome     = mouse.code(sizes$Offset),
  Start_Position = 1L,
  End_Position   = as.numeric(sizes$Size)
)
arms <- arms[order(as.integer(Chromosome))]

fwrite(arms, ref("mm10", "arm_mm10.txt"), sep = "\t", quote = FALSE)
message(sprintf("-- wrote reference/mm10/arm_mm10.txt (%d arms, one q arm per chromosome)",
                nrow(arms)))

## ---------------------------------------------------------------------------
## 2. hrd_chromosome_map.txt
## ---------------------------------------------------------------------------

# Human HRD chromosomes. 18 is included because the pipeline carries chr18 gain
# as a separate feature (it is deliberately excluded from the HRD call itself).
hrd.human <- c(3, 5, 7, 9, 11, 15, 19, 21, 18)

hs <- unique(fread(ref("hg38", "UCSC_hg38_gene_models.txt"), sep = "\t",
                   header = TRUE, select = c("chrom", "geneSymbol"),
                   colClasses = "character"))
hs[, chrom := gsub("^chr", "", chrom)]
hs <- hs[chrom %in% as.character(1:22)]

mm <- unique(fread(ref("mm10", "UCSC_mm10_gene_models.txt"), sep = "\t",
                   header = TRUE, select = c("chrom", "geneSymbol"),
                   colClasses = "character"))
mm[, chrom := gsub("^chr", "", chrom)]
mm <- mm[chrom %in% c(as.character(1:19), "X", "Y")]
setnames(mm, c("mouse_chrom", "mouse_symbol"))

# babelgene: human symbol -> mouse symbol, one-to-one orthologs only, so that a
# gene family cannot inflate one chromosome's share.
orth <- as.data.table(babelgene::orthologs(genes = unique(hs$geneSymbol),
                                           species = "mouse", human = TRUE))
orth <- unique(orth[, .(geneSymbol = human_symbol, mouse_symbol = symbol)])

pairs <- merge(hs, orth, by = "geneSymbol", allow.cartesian = TRUE)
pairs <- merge(pairs, mm, by = "mouse_symbol", allow.cartesian = TRUE)
pairs <- unique(pairs[, .(geneSymbol, chrom, mouse_chrom)])

out <- rbindlist(lapply(hrd.human, function(h) {
  sub <- pairs[chrom == as.character(h)]
  if (!nrow(sub)) {
    return(data.table(human_chrom = h, mouse_chrom = NA_character_,
                      n_orthologs_total = 0L, n_orthologs_on_mouse_chrom = 0L,
                      ortholog_fraction = NA_real_, runner_up = NA_character_,
                      use = 0L))
  }
  tab <- sort(table(sub$mouse_chrom), decreasing = TRUE)
  data.table(
    human_chrom                = h,
    mouse_chrom                = names(tab)[1],
    n_orthologs_total          = nrow(sub),
    n_orthologs_on_mouse_chrom = as.integer(tab[[1]]),
    ortholog_fraction          = round(as.numeric(tab[[1]]) / nrow(sub), 3),
    runner_up                  = if (length(tab) > 1)
      sprintf("chr%s (%.3f)", names(tab)[2], as.numeric(tab[[2]]) / nrow(sub)) else NA_character_,
    use                        = 1L
  )
}))

header <- c(
  "# Human HRD chromosome -> mm10 proxy chromosome, by one-to-one gene orthology.",
  "# Generated by scripts/build_mm10_arms_and_hrd_map.R -- EDIT FREELY.",
  "#",
  "# ortholog_fraction is the share of that human chromosome's one-to-one mouse",
  "# orthologs that sit on `mouse_chrom`. It is the measure of how good (or bad)",
  "# the proxy is. Values well under 1 mean the human chromosome's genes are",
  "# scattered across several mouse chromosomes and a whole-chromosome gain in",
  "# mouse is NOT equivalent to the human trisomy. Review every row.",
  "#",
  "# Set use = 0 to stop calling a given human HRD chromosome in mouse. That",
  "# feature is then reported as 0 for every sample (not NA -- NA would make",
  "# MMGenoCla's complete.cases() filter drop the sample entirely).",
  "#",
  "# The whole map is ignored unless the mouse pipeline is called with",
  "# HRD.PROXY = \"ortholog\"; the default, \"none\", zeroes all of these features.",
  "#"
)
outfile <- ref("mm10", "hrd_chromosome_map.txt")
writeLines(header, outfile)
fwrite(out, outfile, sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE)

message("-- wrote reference/mm10/hrd_chromosome_map.txt")
print(out)

## --- quality report on the mapping itself ------------------------------------

message("\n== How good is this mapping? ==")

weak <- out[!is.na(ortholog_fraction) & ortholog_fraction < 0.5]
if (nrow(weak)) {
  message(sprintf("!! %d of %d human HRD chromosomes have under half their orthologs",
                  nrow(weak), nrow(out)))
  message(sprintf("!! on the chosen mouse chromosome: %s",
                  paste(sprintf("hs%s->mm%s (%.2f)", weak$human_chrom,
                                weak$mouse_chrom, weak$ortholog_fraction),
                        collapse = ", ")))
}

dups <- out[!is.na(mouse_chrom)][, .N, by = mouse_chrom][N > 1]
if (nrow(dups)) {
  for (i in seq_len(nrow(dups))) {
    who <- out[mouse_chrom == dups$mouse_chrom[i], human_chrom]
    message(sprintf("!! COLLISION: human chromosomes %s both map to mouse chr%s --",
                    paste(who, collapse = " and "), dups$mouse_chrom[i]))
    message("!!            a single gain there cannot mean both trisomies.")
  }
}

message("")
message("!! A whole-chromosome gain in mouse is NOT equivalent to the human")
message("!! trisomy it is mapped to here. This is why the mouse pipeline defaults")
message("!! to HRD.PROXY = \"none\": the HRD features are reported as 0 and no")
message("!! HRD-based group is assigned, while the ortholog-anchored GISTIC peak")
message("!! features -- which do transfer -- still drive the classification.")
message("!! Pass HRD.PROXY = \"ortholog\" only if you have reviewed the rows above.")
