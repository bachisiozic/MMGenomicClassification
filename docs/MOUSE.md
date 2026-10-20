# The mouse (mm10) pipeline

For vkMYC and related mouse models of multiple myeloma.

**Read this before using it.** The classification was derived in human myeloma.
Applying it to mouse is a translation, and the translation is sound for some
features and not for others. This page says which is which, with numbers.

```r
source("MMGenomicClassification.R")
source("R/mm10/MMGenomicClassification_mm10.R")
MMGC_setup("mm10")
MMGC_mm10_report()          # prints everything below, for your reference files
```

---

## The design principle

Every feature keeps its **human name and meaning** while being measured at the
**orthologous mouse locus**. `Amp_1q21.1` still means "MCL1 gain" — but is
detected on mouse chr3, where *Mcl1* lives. `CNV.SNV_TP53` is measured at mouse
*Trp53*.

This is inherited from how `reference/mm10` was built (see
`reference/mm10/README.md`) and extended here to gene-level features. The point is
that the eleven groups stay comparable with human cohorts, and `MMGenoCla()` needs
no mouse-specific rules — the mouse classification is the same function on the
same feature names.

---

## What transfers

**Recurrent copy-number peaks.** 55 of the 64 human GISTIC peaks have a confident
mouse anchor. Each was relocated by orthologous gene rather than by lifting
coordinates, because synteny is broken. The nine without an anchor are mostly
immunoglobulin loci (`reference/mm10/QC_unmapped_peaks.txt`); of those only
`Del_12q24.31` survives to the feature set, so the mouse matrix has **117
features against the human 126**.

**Driver genes.** All 32 tumour-suppressor features and 51 of the 59 mutation
features have a one-to-one mouse ortholog. The eight without are the histone
cluster genes (`H1.4`, `H1.5`, `H1.2`, `H1.3`, `H2BC12`), `HLA.C`, `LCE1D` and
`TBC1D29P` — gene families and an MHC locus, none of which has a clean single
ortholog. Five of the eight are excluded from the classification anyway; the
remaining three (`H1.4`, `H1.2`, `HLA.C`) are lost from the heatmap only. **No
classification rule uses any of them**, so the group assignments are unaffected.

**Deletion burden and the 1q / del(13q) axis.** Both are built from the features
above, so both carry over.

---

## What does not transfer: hyperdiploidy

Human HRD is trisomy of chromosomes 3, 5, 7, 9, 11, 15, 19 and 21. There is no
mouse equivalent. Each of those human chromosomes scatters its genes across
several mouse chromosomes, and the best single proxy is poor:

| Human | Best mouse proxy | Share of orthologs there | Runner-up |
|---|---|---|---|
| chr3 | mouse chr9 | 33% | chr16 (28%) |
| chr5 | mouse chr13 | 36% | chr18 (29%) |
| chr7 | mouse chr5 | 44% | chr6 (38%) |
| chr9 | mouse chr4 | 39% | chr2 (37%) |
| chr11 | mouse chr7 | 33% | chr9 (27%) |
| chr15 | mouse chr9 | 42% | chr2 (33%) |
| chr19 | mouse chr7 | 52% | chr8 (17%) |
| chr21 | mouse chr16 | 61% | chr10 (26%) |
| chr18 | mouse chr18 | 76% | chr17 (12%) |

Six of nine are below half. Worse, the mapping is **degenerate**: human chr3 and
chr15 both land on mouse chr9, and human chr11 and chr19 both on mouse chr7. A
single gain of mouse chr9 cannot mean both trisomy 3 and trisomy 15.

So `HYidentification_mm10()` defaults to **`HRD.PROXY = "none"`**: the human HRD
features are set to 0, no HRD-based group is assigned, and the classification runs
on the features that do transfer. Samples resolve into the translocation, 1q /
del(13q), `Multiple_Losses` and `Simple` groups.

The features are set to **0, not `NA`**, deliberately: `MMGenoCla()` drops any
sample with a missing feature, so `NA` here would discard the entire cohort.

### If you want the proxy anyway

```r
HYidentification_mm10(CNV = cnv, DATABASE = db, OUT.DIR = out.dir,
                      HRD.PROXY = "ortholog")
```

It warns, and you should have read the table above first. The map lives in
`reference/mm10/hrd_chromosome_map.txt` and is **meant to be edited** — set
`use = 0` for a chromosome you do not want called, or change `mouse_chrom` if a
different chromosome is the better proxy for your model. Regenerate the defaults
with `Rscript scripts/build_mm10_arms_and_hrd_map.R`.

Mouse-chromosome gains are always written to
`<pct>%_ARMS_<DATABASE>hyperdiploid_mouse_chromosomes.txt` regardless of
`HRD.PROXY`, because that is the measurement actually made. Use it directly if
you want mouse aneuploidy rather than a human-labelled proxy.

---

## Three more caveats

**Arms.** Mouse chromosomes are telocentric — the centromere is at the start and
there is effectively no p arm, which is why every mm10 cytoBand is `qA1`, `qB2`
and so on with no `acen` rows. `reference/mm10/arm_mm10.txt` therefore gives each
chromosome one q arm spanning its whole length, so arm-level and
chromosome-level gain are the same measurement in mouse.

**APOBEC.** The SBS reference signatures are human, and `MMGC_APOBEC_CUTOFFS` are
deciles from a 1,933-sample human cohort. Fitting them to mouse is reasonable for
APOBEC — the TCW context is conserved — but the absolute contributions are not
comparable with human cohorts, and the bin boundaries are not calibrated for
mouse. Treat the APOBEC feature as a rank within your own cohort, or recompute
the cut-offs from it.

**DIS3 hotspots.** `D479`, `D488` and `R780` are **human** DIS3 residue numbers.
Mouse *Dis3* is not guaranteed to number the same residues identically. Check the
alignment for your annotation and pass `HOTSPOTS` explicitly if it differs:

```r
dis3_hs_mm10(MATRIX = matrix, CODING.MUT = muts, HOTSPOTS = c("D477", "D486", "R778"))
```

**Translocations.** IGH translocations do not occur in mouse models the way they
do in human myeloma. Supply `t_CCND1`, `t_MMSET` and `t_MAF` as 0 unless your
model has an engineered equivalent. With all three at 0 and HRD off, the
classification rests on the 1q / del(13q) axis, deletion burden, APOBEC and
chromothripsis.

---

## Running it

`examples/run_mm10.R` is a complete script. The functions mirror the human ones
with an `_mm10` suffix; each sets up the mouse specifics and delegates to the
shared implementation, so there is no third copy of the pipeline to drift out of
step.

| Mouse | Differs from the human function by |
|---|---|
| `inputCNV_mm10()` | mouse chromosome coding (1–19, X = 20, Y = 21) |
| `FocalLargeMatrix_mm10()` | nothing — mouse peaks carry human codes |
| `HYidentification_mm10()` | scores mouse chromosomes, then applies `HRD.PROXY` |
| `CytogeneticMatrixGenerator_mm10()` | nothing |
| `MutationsMatrixGenerator_mm10()` | translates mouse gene symbols to human orthologs |
| `APOBECanalyzer_mm10()` | mm10 trinucleotide contexts |
| `FinalMatrixGenerator_mm10()` | the reduced 117-feature set |
| `dis3_hs_mm10()` | translates `Dis3`, warns about residue numbering |
| `MMGenoCla_mm10()` | reports whether any HRD group was assigned |

### Prerequisites

```bash
Rscript scripts/install_dependencies.R mm10
```

Mutation calls need `dndscv` against a **mouse** RefCDS, with mouse gene symbols
in the `gene` column — the pipeline translates them to human orthologs itself.

---

## What is missing from the reference set

`reference/mm10/README.md` refers to `../script/build_mm10_references.R`, which
built the relocated GISTIC peaks and the mouse driver table. **That script is not
in this repository** — the reference files it produced are here, but not the code
that made them. `scripts/build_mm10_arms_and_hrd_map.R` covers only the two tables
added here (`arm_mm10.txt` and `hrd_chromosome_map.txt`). If you have the original
builder, adding it would make the mouse reference set fully reproducible.
