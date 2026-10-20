# Input formats

Four inputs. Column *names* are yours to choose where a function takes them as
arguments; where a file is read directly, the names are fixed.

---

## 1. Copy number — passed to `inputCNV()`

Allele-specific segments, one row per segment. You name the columns:

```r
cnv <- inputCNV(my.table,
                SAMPLE.COL = "sample",
                CHR.COL    = "chr",
                START.COL  = "startpos",
                END.COL    = "endpos",
                TOTAL.COL  = "tot",     # total copy number
                MINOR.COL  = "min")     # minor allele copy number
```

| Field | Meaning |
|---|---|
| sample | sample identifier; must match the other inputs exactly |
| chr | `1`–`22`, `X`, `Y`, with or without a `chr` prefix |
| start, end | 1-based segment bounds |
| total | total copy number, may be fractional |
| minor | minor allele copy number, may be fractional; `NA` allowed |

`inputCNV()` recodes chromosomes to integers (human `X` = 23, `Y` = 24; mouse
`X` = 20, `Y` = 21) and drops anything outside that set — scaffolds, alts, chrM —
reporting what it dropped.

**Fractional values are used as-is against fixed thresholds** (gain at
`total >= 2.5`, amplification at `>= 3.5`, loss at `< 1.5`, biallelic at `< 0.5`).
They are not rescaled for purity or ploidy, so pass purity-corrected,
ploidy-aware calls — ASCAT, FACETS, Battenberg or equivalent. Purity that is
left uncorrected shifts every segment towards 2 and loses events; a sample left on
raw read-depth ratios can come out "gained" genome-wide.

**Segments must be disjoint** within a sample and chromosome. Gained bases are
summed per arm, so overlapping segments double-count them and inflate the gained
fraction. Allele-specific WGS callers emit disjoint segments; overlap usually means
the table was merged from more than one source.

**Run `checkCNV()` on the output of `inputCNV()`** before anything else — it checks
all of the above per sample and says which ones are safe to analyse.

**`minor` matters.** Loss severity depends on it: a segment with
`total >= 1.5` and `minor < 0.5` is a copy-neutral LOH-style single-copy loss.
Where `minor` is `NA` and `total >= 2`, the segment cannot be called and the peak
becomes `NA` rather than 0.

---

## 2. Coding mutations — passed to `MutationsMatrixGenerator()` and `dis3_hs()`

The [dndscv](https://github.com/im3sanger/dndscv) `annotmuts` layout, filtered to
non-synonymous:

```r
dnds <- dndscv(mutations, refdb = "hg19")      # or an hg38/mm10 RefCDS
out  <- dnds$annotmuts[dnds$annotmuts$impact != "Synonymous", ]
write.table(out, "cohort_dnds_ns_mutations.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

| Column | Required | Notes |
|---|---|---|
| `sampleID` | yes | must match the copy-number sample ids |
| `gene` | yes | HGNC symbol; mouse symbols for mm10 |
| `impact` | yes | rows equal to `Synonymous` are dropped |
| `aachange` | for `dis3_hs()` | e.g. `D479N`; needed for DIS3 hotspots |

`dndscv` input itself is `sampleID`, `chr`, `pos`, `ref`, `mut`.

hg19-era symbols are translated automatically, so `FAM46C`, `HIST1H1E`,
`HIST1H1C`, `HIST1H1D`, `HIST1H1B`, `HIST1H2BK` and `TBC1D29` are accepted
alongside `TENT5C`, `H1.4`, `H1.2`, `H1.3`, `H1.5`, `H2BC12` and `TBC1D29P`.

---

## 3. Metadata — passed to `MutationsMatrixGenerator()`

Which assays exist per sample. This is what separates "wild type" from "not
assessed": without it every gap becomes a 0 and absent data reads as absent
mutation.

```
sample      cnv   snv
SAMPLE_01   1     1
SAMPLE_02   1     0
```

| Column | Meaning |
|---|---|
| `sample` | sample identifier |
| `cnv` | 1 if copy number was called for this sample, else 0 |
| `snv` | 1 if mutations were called, else 0 |

A sample with `snv = 0` gets `NA` across the mutation features, not 0. Since
`MMGenoCla()` drops samples with any missing feature, such a sample will not be
classified — which is the intended behaviour, not a failure.

Build it from whatever you have:

```r
meta <- merge(
  unique(data.frame(sample = unique(cnv.calls$sample), cnv = 1)),
  unique(data.frame(sample = unique(snv.calls$sample), snv = 1)),
  by = "sample", all = TRUE)
meta[is.na(meta)] <- 0
```

---

## 4. SNVs — passed to `APOBECanalyzer()`

Single-base substitutions, columns named by you:

```r
APOBECanalyzer(SNV.OBJECT = snv, SAMPLE = "sample", CHROM = "chrom",
               POS = "pos", REF = "ref", ALT = "alt",
               ITER = 1000, DATABASE = "MyCohort", OUT.DIR = out.dir)
```

Rows where `ref` and `alt` are not both single bases are dropped — indels do not
contribute to SBS signatures. Coordinates must be on the build passed to
`MMGC_setup()`: the trinucleotide context is read from the reference genome, so a
build mismatch silently produces wrong contexts.

`ITER` is the bootstrap count. The published analysis used 1000; the default here
is 10, which is fast enough for a check but too few for reported confidence
intervals.

---

## 5. What the pipeline does not compute

Merge these into the final matrix before calling `MMGenoCla()`.

### Translocations

| Column | Lesion | Accepted aliases |
|---|---|---|
| `t_CCND1` | t(11;14) | `t(11;14)` |
| `t_MMSET` | t(4;14) | `t_NSD2`, `t(4;14)` |
| `t_MAF` | t(14;16) or t(14;20) | `t(14;16)`, `t(14;20)` |

Code 0 or 1; `MMGenoCla()` rescales 1 to 2 internally. A sample with no
translocation is 0, not `NA`.

### Chromothripsis

| Column | Values | Accepted aliases |
|---|---|---|
| `CNV.Sig` | 0 or 1 | `Chromothripsis` |

---

## Putting it together

```r
matrix <- read.delim("FinalMatrix_MyCohort.txt", stringsAsFactors = FALSE)

matrix <- merge(matrix, translocations, by = "sample", all.x = TRUE)
matrix <- merge(matrix, chromothripsis,  by = "sample", all.x = TRUE)

dis3   <- dis3_hs(MATRIX = matrix, CODING.MUT = "cohort_dnds_ns_mutations.txt")
matrix <- merge(matrix, dis3, by = "sample")

# Absent translocation / chromothripsis calls mean "not present", so fill with 0.
# Do NOT blanket-fill the assay-driven features -- an NA there means not assessed
# and should stay NA.
for (col in c("t_CCND1", "t_MMSET", "t_MAF", "CNV.Sig")) {
  matrix[[col]][is.na(matrix[[col]])] <- 0
}

res <- MMGenoCla(input = matrix)
```

---

## Sample identifiers

The single most common cause of empty output is identifiers that do not match
between inputs — a `IID_` prefix on one, a `_1_BM` suffix on another. Check
before you start:

```r
Reduce(intersect, list(unique(cnv$sample),
                       unique(mutations$sampleID),
                       unique(snv$sample),
                       meta$sample))
```

If that is shorter than you expect, normalise the ids first. The pipeline joins
on exact string equality throughout.
