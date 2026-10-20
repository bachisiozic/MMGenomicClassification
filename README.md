# MMGenomicClassification

Genomic classification of multiple myeloma into eleven groups, from copy number,
driver mutations, translocations and mutational signatures.

Built for **whole-genome sequencing**: allele-specific, purity- and
ploidy-corrected copy number (ASCAT, Battenberg, FACETS) plus somatic SNVs.
`checkCNV()` verifies your input meets that before you run anything — see
[Input QC](#input-qc).

This is the reference implementation of the classification published in:

> Maura F, Rajanna AR, Ziccheddu B, et al. **Genomic Classification and Individualized
> Prognosis in Multiple Myeloma.** *Journal of Clinical Oncology*.
> doi:[10.1200/JCO.23.01277](https://doi.org/10.1200/JCO.23.01277)

Supports **hg19**, **hg38** and **mm10** (mouse models — read
[docs/MOUSE.md](docs/MOUSE.md) first, the translation is only partly sound).

---

## Quick start

```bash
git clone https://github.com/<your-account>/MMGenomicClassification.git
cd MMGenomicClassification
Rscript scripts/install_dependencies.R hg38
Rscript scripts/selftest.R
```

`selftest.R` runs the whole pipeline on synthetic data and needs no input of your
own. If it passes, your installation works.

Then, in R:

```r
source("MMGenomicClassification.R")
MMGC_setup("hg38")            # or "hg19"
```

That one call resolves every reference file for the build, so you never type a
reference path again.

---

## The pipeline

Eight steps. Each writes its output into `OUT.DIR` and the next one reads it from
there, so `OUT.DIR` and `DATABASE` must stay the same throughout.

| Step | Function | Input | Produces |
|---|---|---|---|
| 1 | `inputCNV()` | your allele-specific copy-number calls | normalised segments |
| 1b | `checkCNV()` | step 1 | input QC — run this first |
| 2 | `FocalLargeMatrix()` | step 1 | recurrent peak scores, focal vs large |
| 3 | `HYidentification()` | step 1 | hyperdiploidy call |
| 4 | `CytogeneticMatrixGenerator()` | steps 2, 3 | HRD chromosomes, peaks, 1q |
| 5 | `MutationsMatrixGenerator()` | coding mutations + step 1 | per-gene driver matrices |
| 6 | `APOBECanalyzer()` | your SNVs | APOBEC (SBS2+SBS13) contribution |
| 7 | `FinalMatrixGenerator()` | steps 4–6 | the feature matrix |
| 8 | `dis3_hs()` → `MMGenoCla()` | step 7 + translocations + chromothripsis | **the classification** |

A complete, runnable script is in
[examples/run_hg38.R](examples/run_hg38.R). See
[docs/PIPELINE.md](docs/PIPELINE.md) for what each step actually computes and
[docs/INPUT_FORMATS.md](docs/INPUT_FORMATS.md) for the exact columns expected.


### What you must supply yourself

The pipeline computes most features but not all. Three come from elsewhere and
have to be merged into the matrix before `MMGenoCla()`:

- **`t_CCND1`, `t_MMSET`, `t_MAF`** — IGH translocation calls, from your own SV
  caller. `t_NSD2` and the cytogenetic names `t(11;14)`, `t(4;14)`, `t(14;16)`,
  `t(14;20)` are accepted as aliases.
- **`CNV.Sig`** — chromothripsis, 0/1. `Chromothripsis` is accepted as an alias.
- **Coding mutations** — the pipeline expects
  [dndscv](https://github.com/im3sanger/dndscv) `annotmuts` output, filtered to
  non-synonymous.

`MMGenoCla()` drops any sample with a missing feature.

---

## The eleven groups

| Group | Defined by |
|---|---|
| `HRD_RAS` | hyperdiploid, RAS mutations, otherwise quiet |
| `HRD_Gains` | hyperdiploid with additional recurrent gains |
| `HRD_Complex_Cytogenetic` | hyperdiploid with a high deletion burden |
| `CCND1_Simple` | t(11;14), low complexity |
| `CCND1_Complex_Cytogenetic` | t(11;14) with deletions, 1q, TP53 loss or APOBEC |
| `MAF_and_or_HyperAPOBEC` | MAF translocation and/or high APOBEC |
| `NSD2_GainAmp1q_Del13q` | t(4;14) with 1q gain and del(13q) |
| `NSD2_HRD` | t(4;14), other |
| `NSD2_Del13q` | t(4;14) with del(13q), no 1q |
| `GainAmp1q_Del13q` | 1q gain with del(13q), no canonical translocation |
| `Multiple_Losses` | high deletion burden or chromothripsis |
| `Simple` | none of the above |

`MMGenoCla()` returns everything needed to draw the published heatmap:

```r
res <- MMGenoCla(input = matrix)

pheatmap::pheatmap(t(res$genomic),
                   annotation_col    = res$annotation,
                   annotation_colors = res$colors,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   col = c("grey85", "salmon1", "darkred"),
                   show_colnames = FALSE, border_color = NA, legend = FALSE)
```

---

## Repository layout

```
MMGenomicClassification.R     entry point -- source this
R/                            the pipeline, one file per step (13-checkCNV.R = input QC)
R/mm10/                       mouse-specific layer
reference/                    reference data, per build
  hg19/ hg38/ mm10/
scripts/
  install_dependencies.R      install what you need
  selftest.R                  end-to-end check on synthetic data
  test_legacy_switches.R      regression test for the fixed defects
  compare_matrices.R          cell-by-cell diff of two feature matrices
  build_slim_gene_references.R      regenerate the gene models
  build_mm10_arms_and_hrd_map.R     regenerate the mouse-specific tables
examples/                     runnable scripts, one per build
docs/                         formats, pipeline detail, validation, known issues, mouse
legacy/                       the original scripts, exactly as published
```

### Reference data

Everything needed ships in `reference/` — about 11 MB, no external downloads and
no Git LFS. The UCSC gene tables are shipped as the four columns the pipeline
reads (`chrom`, `txStart`, `txEnd`, `geneSymbol`). Row order is preserved, so the transcript chosen
per gene is unchanged; `scripts/build_slim_gene_references.R` regenerates them
from the full UCSC dumps.

---

## Relationship to the published scripts

`legacy/` holds the original `JCO_GenomicClassification.R` (hg19),
`JCO_GenomicClassification_hg38.R` and `RUN_JCO_Genomic_Classification.R`
unmodified, so published results stay reproducible.

The version in `R/` is the same analysis, reorganised.


## Citation

If you use this, please cite the paper:

```bibtex
@article{Maura2023MMGenomicClassification,
  author  = {Maura, Francesco and Rajanna, Arjun R and Ziccheddu, Bachisio and others},
  title   = {Genomic Classification and Individualized Prognosis in Multiple Myeloma},
  journal = {Journal of Clinical Oncology},
  doi     = {10.1200/JCO.23.01277}
}
```

See [CITATION.cff](CITATION.cff).

## License

[MIT](LICENSE). The reference data in `reference/` derives from UCSC, COSMIC and
GISTIC outputs and carries the terms of those sources.
