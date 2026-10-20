# Validation against the published feature matrix

The pipeline was run end to end on the cohort behind
[doi:10.1200/JCO.23.01277](https://doi.org/10.1200/JCO.23.01277) — 1,834
copy-number profiles, 2,046 mutation profiles, 511,972 substitutions — and its
output compared cell by cell with the published feature matrix
(`PMMM_matrix_04052023_NEW_xArjun_COL_SELECTED_ORIGINAL.txt`, 1,933 samples ×
133 columns).

**Headline: 98.78% of comparable cells are identical out of the box** (209,292 of
211,878), rising to **99.03%** with `REQUIRE.MINOR.ALLELE = "non-HRD"` — see
[the minor-allele question](#what-the-comparison-found-the-minor-allele-question),
which is the one setting that matters here.

A "comparable" cell is one where both matrices carry a value. Cells the published
matrix leaves as `NA` record assays that were not done, which the pipeline cannot
and should not invent — see [Coverage](#coverage) below.

## By feature group

Measured with the defaults, i.e. `REQUIRE.MINOR.ALLELE = TRUE` (the published code
path). The last column is the same measurement with `"non-HRD"`, which is the only
setting that changes any of these numbers.

| Feature group | Cells | Identical | Agreement | with `"non-HRD"` |
|---|---|---|---|---|
| Driver mutations (`SNV_*`) | 95,680 | 95,673 | **99.99%** | 99.99% |
| APOBEC (`0/1/2` bins) | 1,526 | 1,524 | **99.87%** | 99.87% |
| Loss peaks (`CNV_Del_*`, `chr8p.loss`) | 31,718 | 31,532 | 99.41% | 99.41% |
| Gain peaks (`CNV_Amp_*`, `Gain_Amp1q`) | 13,541 | 13,269 | 97.99% | 97.99% |
| Tumour suppressors (`CNV.SNV_*`) | 54,176 | 52,592 | 97.08% | 97.08% |
| HRD trisomies (`CNV_chr*.gain`) | 15,237 | 14,702 | 96.49% | **99.99%** |
| **All** | **211,878** | **209,292** | **98.78%** | **99.03%** |

988 of 1,739 samples (57%) match on every single comparable feature; the median
per-sample agreement is 100%, and 43 samples fall below 90%. With `"non-HRD"` those
become 1,045 samples (60%) and 25 respectively.

Two results are worth calling out because they test the parts that were rewritten
most heavily:

- **APOBEC at 99.87%.** This feature comes from mutational-signature fitting, which
  was reimplemented from scratch to take the genome as a parameter instead of
  depending on a privately patched copy of mmsig (see
  [KNOWN_ISSUES.md](KNOWN_ISSUES.md#3-hg38-mutational-signatures-were-fitted-against-hg19-contexts)).
  Fitting 437,723 substitutions across 1,727 samples reproduces the published
  `APOBEC` bins for 1,524 of 1,526 samples — the two disagreements are off by one
  bin, none by two.
- **Mutations at 99.99%** (7 differing cells out of 95,680), which confirms the
  driver-role handling, the hg19→current symbol translation, and the
  gene-by-gene pairing that replaced the hardcoded `1:36` / `73:108` column
  positions.

## What the comparison found: the minor-allele question

The single largest discrepancy, before it was resolved, was hyperdiploidy. The
published `HYidentification()` selects gained segments with

```r
cnv.gain = cnv.n[cnv.n$total >= 2.5 & cnv.n$minor >= 0.5, ]
```

— it requires the **second allele to be retained**, so a chromosome at three
copies with the minor allele lost (`3+0`) does not count as a trisomy.

Measured against the published matrix, that is not how the matrix was made:

| Gain definition | Trisomy calls | Agreement with the matrix |
|---|---|---|
| `total >= 2.5 & minor >= 0.5` (published code) | 4,878 | 96.49% — misses 535 trisomies |
| `total >= 2.5` (total copy number only) | 5,413 | **99.99%** — 1 cell either way |

The matrix contains exactly 5,413 positive trisomy calls; total-copy-number-only
scoring reproduces 15,235 of 15,237 sample-chromosome calls. Every one of the 535
trisomies the allele-aware rule misses is either a `3+0` gain or a sample with no
allele-specific data at all. For example `MMRF_1719` is at `total = 3, minor = 0`
across the whole of chromosome 9 — three copies, one allele — which the matrix
counts and the published filter does not.

This also matters practically. Allele-specific copy number is simply absent for
some cohorts: **all 159 Moffitt samples in this dataset have `minor = NA` on every
segment**, and requiring the minor allele makes every one of them look diploid
with no warning.

`HYidentification()` therefore takes `REQUIRE.MINOR.ALLELE`, with three settings:

| Setting | Trisomie | Double_Genomes | Hyperdiploid | Trisomy calls | Overall |
|---|---|---|---|---|---|
| `TRUE` — **default**, published code path | 913 | 56 | 820 | 96.49% | 98.78% |
| `"non-HRD"` — total CN for trisomies only | 1,017 | **62** | 905 | **99.99%** | **99.03%** |
| `FALSE` — never require it | 1,016 | 63 | 905 | 99.99% | 99.03% |

**The default is `TRUE`, so the function behaves as published unless asked
otherwise.** Know what that costs: 535 of the published matrix's trisomies are not
reproduced, and any sample whose copy-number calls carry no minor allele is
reported as diploid. `HYidentification()` warns when it finds such samples — on
this cohort it reports all 159 Moffitt samples.

`"non-HRD"` is what reproduces the published matrix. It scores a trisomy on total
copy number alone but keeps gains elsewhere allele-aware, which also keeps the
whole-genome-doubling call strict: doubling means gaining *both* alleles, so a
`3+0` gain should not count towards it.

`"non-HRD"` and `FALSE` produce an **identical feature matrix** — verified at
222,592 of 222,592 cells — differing only in that one Double_Genomes label. The
reason they tie is that `HYidentification()`'s only contribution to the matrix is
the `hyp.chr*` trisomy flags, which exist solely for the HRD chromosomes; gains
elsewhere are scored but feed only `Cytogenetic` and `double.chrs`, which nothing
downstream reads.

Note also that the rest of the pipeline already makes this split without a switch:
`FocalLargeMatrix()` scores peak-level **gains** on total copy number alone — the
allele-aware variant is the one commented out in the published source — while its
**losses**, and the gene-level deletions in `MutationsMatrixGenerator()`, are
allele-aware. `REQUIRE.MINOR.ALLELE` governs trisomy calls only.

## Residual disagreement

2,586 cells (1.22%) with the defaults, or 2,053 (0.97%) with `"non-HRD"`; the
difference is the 535 trisomy calls above. Setting that aside, the rest is not
spread evenly:

- **988 samples have none at all** (1,045 with `"non-HRD"`). The top 100 samples
  hold roughly half of it.
- **96% is one-directional**: the published matrix scores one hit more than the
  pipeline (1,126 cells where it has 1 and we have 0, 401 where it has 2 and we
  have 1).
- Most of it sits in the tumour-suppressor block, spread evenly across all 32
  genes (median 47.5 cells per gene, range 24–93). It is not concentrated on the
  four genes that also have a GISTIC peak proxy (16% of the disagreement across
  4 of 32 features, against 12% expected if uniform), so it is not a peak-vs-gene
  artefact.
- The rate tracks allele-specific data availability, exactly as the trisomy
  finding does:

  | Segments with `minor = NA` | Samples | TSG disagreements per sample |
  |---|---|---|
  | none | 924 | 0.42 |
  | under 25% | 556 | 0.94 |
  | 25–50% | 53 | 2.28 |
  | all | 159 | 3.45 |

So the remaining gap is the same phenomenon one level down: the gene-level
copy-number call in the published matrix is more permissive about lost minor
alleles than the allele-aware rule in the published code. The gene-level call is
left as published, since unlike the trisomy case there is no single alternative
rule that closes the gap cleanly.

## What `checkCNV()` says about this cohort

The input QC, run on the same 1,693 copy-number profiles, explains most of the
residual disagreement before the pipeline is even run:

| Finding | Samples |
|---|---|
| clean | 918 |
| partial allele-specific data (median 10.6% of segments, worst 52.9%) | 541 |
| no allele-specific data at all | 159 |
| overlapping segments | 46 |
| implausible ploidy (outside 1.2–5) | 20 |
| malformed segments (zero width) | 9 |

Three things worth noting, because they are properties of the data rather than of
the code:

- **All 159 Moffitt samples have no minor allele on any segment**, and every one of
  them also has overlapping segments — autosomal coverage up to 2.5× the genome.
  MPG, MSK and UAMS are free of overlaps.
- **20 samples have an implausible ploidy.** `UAMS_35830` sits at total copy number
  23–36 across whole chromosomes and `UAMS_26497` at 0 — calls that were not
  normalised for purity and ploidy. Both are scored genome-wide by absolute
  thresholds, so they come out gained (or lost) everywhere.
- The QC recommends `REQUIRE.MINOR.ALLELE = "non-HRD"` for this cohort, because of
  the 159. For a clean allele-specific WGS cohort it recommends `TRUE`, the default.

## Coverage

| | Target | Pipeline |
|---|---|---|
| Samples | 1,933 | 1,739 |
| with copy number | 1,693 | 1,693 ✓ |
| with mutations | 1,727 | 1,727 ✓ |

Copy-number availability matches the published matrix exactly — 1,693 samples,
zero mismatches in either direction. The 194 samples the pipeline does not emit
have neither copy number nor mutations in the supplied inputs; the published
matrix carries them as all-`NA` rows.

Of the 10,438 raw cell differences, **8,385 (80%) are `NA`-versus-value**, not
disagreeing values. Most of that is one thing the inputs do not describe: **MSK
and UAMS are panel-sequenced**. 108 MSK samples have exactly 29 `NA` mutation
columns and 93 UAMS samples exactly 33 — a per-cohort gene mask marking genes that
were never on the panel. No file in the inputs encodes it, so the pipeline
correctly reports `0` (sequenced, not mutated) where the published matrix reports
`NA` (not assayed). If you are assembling a matrix from panel data, apply that
mask yourself after `FinalMatrixGenerator()`.

## Reproducing this

The comparison needs the cohort's own inputs, which are not redistributed here.
Given a reference matrix and your own run, `scripts/compare_matrices.R` performs
the same cell-by-cell comparison:

```bash
Rscript scripts/compare_matrices.R reference_matrix.txt FinalMatrix_MyCohort.txt
```

Two things the comparison relies on, which are worth repeating for any such check:

- **Sample identifiers.** The inputs used three different id spaces (canonical,
  `cnv.sample`, `snv.sample`), the mapping is many-to-many, and a named-vector
  lookup silently drops samples — it has to be a join. Mutation ids also needed a
  `_1_BM` suffix and a ` (SAMEA…)` parenthetical stripped.
- **Genome build.** The cohort is hg19; segment ends exceed hg38 chromosome
  lengths on 14 of 23 chromosomes. Running it as hg38 would produce plausible
  nonsense.

## One regression this exercise caught

Running on real data found a bug introduced by the refactor, now fixed. The gain
filter was written as

```r
cnv.n[!is.na(cnv.n$chrom) & cnv.n$total >= 2.5 & cnv.n$minor >= 0.5, ]
```

but a logical index carrying `NA` does not skip those rows — it inserts an all-`NA`
row for each one. With 22,820 segments having an unknown minor allele, that
injected 3,767 phantom rows and `GRanges()` then refused the `NA` coordinates. The
published code hit the same thing and swept the phantom rows up afterwards with
`complete.cases()`; the fix selects through `which()`, which never creates them.
Every logical row-subset in the pipeline was audited for the same pattern.
Synthetic test data had no missing minor alleles, so `scripts/selftest.R` could not
have caught this — real data did.
