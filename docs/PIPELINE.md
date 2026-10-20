# What each step computes

Read [INPUT_FORMATS.md](INPUT_FORMATS.md) first for the inputs. This is what
happens to them.

Throughout, `DATABASE` is your cohort label and `OUT.DIR` the working folder.
Steps find each other's output by those two, so keep them constant across a run.

---

## 1. `inputCNV()` — normalise segments

Selects and renames your columns to `sample chrom start end total minor`, recodes
chromosomes to integers, drops non-primary contigs and sorts. No thresholds
applied yet.

## 1b. `checkCNV()` — is the input fit for this pipeline?

Run before anything else. Every threshold downstream is an **absolute copy
number**, so the pipeline assumes calls already corrected for purity and ploidy,
with the minor allele present. `checkCNV()` reports, per sample: allele-specific
completeness, overlapping segments, autosomal coverage, length-weighted ploidy,
and malformed segments — then gives one verdict per sample and recommends the right
`REQUIRE.MINOR.ALLELE` setting.

Three failures it exists to catch, all otherwise silent:

| Problem | Consequence |
|---|---|
| `minor = NA` on a segment | can never satisfy `minor >= 0.5`, so it never counts towards a gain; the sample comes back non-hyperdiploid rather than unassessable |
| overlapping segments | gained bases are summed per arm, so overlap double-counts them — a chromosome can pass the 60% threshold on far less than 60% of real coverage |
| uncorrected ploidy | a sample at total copy number 20 is "gained" everywhere and is hyperdiploid automatically |

It reports `minor > major` counts too, but does **not** treat them as errors:
with fractional copy number `total = 1.5, minor = 1` is a subclonal loss, and
nothing in the pipeline uses the major allele.

## 2. `FocalLargeMatrix()` — recurrent peaks, focal vs large

Overlaps every segment with the recurrent GISTIC peaks in
`reference/<build>/Genomic_Coordinates_GISTICpeaks_<build>.txt`, then scores each
peak per sample. A segment is **focal** at or below `FOCAL.LENGTH` (5 Mb by
default) and **large** above it — the same absolute copy-number change means more
when it is focal, because a focal event implicates the peak's target gene rather
than the whole arm.

| Gains (`lesion == "Amp"`) | Score | | Losses (`lesion == "Del"`) | Score |
|---|---|---|---|---|
| `total >= 3.5`, focal | 2 | | `total < 0.5`, focal | 2 |
| `total >= 3.5`, large | 1.5 | | `total < 0.5`, large | 1.5 |
| `2.5 <= total < 3.5`, focal | 1 | | one copy lost, focal | 1 |
| `2.5 <= total < 3.5`, large | 0.5 | | one copy lost, large | 0.5 |
| otherwise | 0 | | otherwise | 0 |

Where several segments hit one peak, the highest score wins. A peak no sample
touched still gets a column, as `NA`.

→ `CNV_<DATABASE>matrix_focal_large_5Mb.txt`

## 3. `HYidentification()` — hyperdiploidy

Consecutive segments with equal rounded copy number are collapsed
(`SegmentCollapsing()`), then segments at `total >= 2.5` with `minor >= 0.5` —
true gains, not LOH — are summed per arm and per chromosome and divided by the arm
or chromosome length.

An arm or chromosome is **gained** above `PERCENTAGE` (0.6). A sample is
**hyperdiploid** when at least two of chromosomes 3, 5, 7, 9, 11, 15, 19, 21 are
gained. chr18 is reported but excluded from the call. Above
`DOUBLE.GENOME.ARMS` (20) gained arms the sample is labelled `Double_Genomes`
instead — whole-genome doubling, not hyperdiploidy.

Requiring 60% of the chromosome is what separates a trisomy from a large
segmental gain.

→ `60%_ARMS_<DATABASE>hyperdiploid.txt`, `...no.hyperdiploid.txt`

## 4. `CytogeneticMatrixGenerator()` — three feature tables

Collapses the focal/large scale to severity — `{0.5, 1} → 1`, `{1.5, 2} → 2` —
and splits the peaks by role:

- `HRD_chromosomes<DATABASE>.txt` — per-chromosome trisomy flags. Samples absent
  from the hyperdiploid table are 0.
- `gistic.other_All_peaks<DATABASE>.txt` — the recurrent gains and losses.
- `x1q_amp_gain<DATABASE>.txt` — the three 1q peaks, kept separate because 1q is
  collapsed into a single feature later.

## 5. `MutationsMatrixGenerator()` — per-gene drivers

Keeps genes called by more than one method or flagged as known myeloma drivers,
splits them by role, and builds:

- **Oncogenes** — mutation only, 0/1.
- **Tumour suppressors** — mutation *and* deletion. A gene is 1 for one hit and 2
  for two, so biallelic loss and deletion-plus-mutation score alike. This is the
  two-hit logic, and it is why `minor` copy number matters upstream.
- **Unknown role** — mutation only.
- **X-linked** (`ATRX`, `KDM5C`, `KDM6A`, `RPL10`) — mutation only, because a
  single-copy loss on X is not interpretable the way it is on an autosome.

Copy number over a gene is taken from segments overlapping its transcript
± 50 kb; the worst call wins. `0` means retained, `1` one copy lost, `2` both.

`NA` versus `0` comes from the metadata: a sample without mutation calls is `NA`
across mutation features, and one without copy number is `NA` across deletion
features. Where only one of the two is missing, the other one's call stands; only
when both are missing is the combined feature `NA`.

→ `oncogene_`, `tsg_mutated_genes_`, `tsg_deleted_genes_`, `tsg_`,
`unk_mutated_genes_`, `chrx_tsg_`

## 6. `APOBECanalyzer()` — mutational signatures

Fits SBS1, SBS2, SBS5, SBS8, SBS9, SBS13 and SBS18 with mmsig, taking
trinucleotide contexts from the build set by `MMGC_setup()`. APOBEC is the
combined relative contribution of **SBS2 + SBS13**.

The result is a continuous value here; `FinalMatrixGenerator()` bins it.

→ `APOBEC_relative_contribution_matrix<DATABASE>.txt`

## 7. `FinalMatrixGenerator()` — the feature matrix

Merges everything from steps 4–6 and reduces it to the published 126 features.

Three things happen here beyond merging:

- **1q collapses** to one feature: gain in any of `Amp_1q21.1`, `Amp_1q21.3`,
  `Amp_1q44` → 1, amplification in any → 2.
- **8p collapses** to `CNV_chr8p.loss` from `Del_8p21.3` and `Del_8p23.1`.
- **APOBEC is binned** on the 6th and 10th deciles of the 1,933-sample cohort the
  classification was derived on — `< 0.0289` → 0, `< 0.1120` → 1, else 2. These
  are fixed, not recomputed from your cohort, so a new cohort is scored on the
  published scale. A small or unusual cohort will not have 60% of samples at 0.

Peaks whose signal a gene-level feature already carries are dropped here
(`Del_1p32.3`, `Del_1p12`, `Del_1p22.1`, `Del_13q14.2`, `Del_16p13.3`,
`Del_16q12.1`, `Del_17p13.1`), as are the immunoglobulin loci
(`Amp_14q32.33`, `Del_14q32.33`, `Del_22q11.22`, `Amp_2p11.2`, `Del_2p11.2`),
where copy-number change is physiological rearrangement rather than a driver.

Feature **order** is part of the contract: `MMGenoCla()` derives the deletion
burden from the loss and tumour-suppressor block, so those stay contiguous.

→ your `OUTPUT` path

## 8. `dis3_hs()` then `MMGenoCla()` — the classification

`dis3_hs()` splits `SNV_DIS3` into `SNV.HS_DIS3` and `SNV.noHS_DIS3`. The three
RNB-domain hotspots — **D479, D488, R780** — sit in the catalytic site and
abolish exonuclease activity; other DIS3 mutations do not carry the same weight,
and the hotspot state is used to exclude samples from the MAF/APOBEC groups.

`MMGenoCla()` then derives four composite features:

| Feature | Definition |
|---|---|
| `HRD` | ≥ 2 of the eight trisomy chromosomes gained |
| `add_gains` | ≥ 3 of the secondary gain peaks |
| `deletions` | ≥ 5 hits across the 50 loss and tumour-suppressor features (10%) |
| `x1q` | the collapsed 1q feature |

and applies the group rules **in order** — later rules overwrite earlier ones,
which is how a sample matching several definitions lands in the most specific
group. Broadly: hyperdiploid groups first, then MAF, then CCND1, then
NSD2/1q/del(13q), then hyper-APOBEC reassignments, then the remainder split by
deletion burden, then a few late corrections.

Samples with **any** missing feature are dropped before classification
(`complete.cases`), so gaps in the translocation, chromothripsis or APOBEC calls
cost whole samples. The count is reported.

Returns `$genomic` (matrix), `$annotation` (group per sample), `$colors`
(palette), `$clusters` (internal rule labels), `$missing.features` and
`$high.missing.features`.

---

## Runtime

Dominated by two things: `APOBECanalyzer()` scales with `ITER` × samples
(`ITER = 1000` on a few hundred samples takes hours), and the overlap steps scale
with the number of segments. Everything else is minutes. Copy number and
signatures are independent, so they can run in parallel.
