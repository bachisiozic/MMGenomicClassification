# Known issues and differences from the published scripts

`legacy/` holds the scripts exactly as published. The implementation in `R/` is
the same analysis reorganised, but it is **not bug-for-bug identical**. This file
lists every behavioural difference, what it affects, and how to switch back.

Three of these change results. Each defaults to the behaviour that matches the
published feature matrix and each has an opt-out, so output from the original
scripts remains reproducible.

The pipeline has been run against the published matrix for the paper's own cohort
and agrees on 99.03% of comparable cells — see [VALIDATION.md](VALIDATION.md),
which is also where items 1 and 2 below were measured.

---

## 1. Hyperdiploidy was called by comparing text, not numbers

**Changes results. Opt out with `HYidentification(LEGACY.STRING.COMPARISON = TRUE)`.**

In the published `HYidentification()`, the per-arm and per-chromosome gain
fractions were accumulated with `rbind()` on `c(...)` vectors. `c()` coerces mixed
types to character, so `hyper.new.1` was a table of *strings*, and

```r
hyper.new.2$arm.gain = ifelse(hyper.new.2$seg.gain.prop > PERCENTAGE, 1, 0)
hyper.new.2$chr.gain = ifelse(hyper.new.2$perc.gain.chr > PERCENTAGE, 1, 0)
```

compared a character vector to `0.6` — a lexicographic comparison, not a numeric
one. The same applies to the `hyp.chrN` flags computed from `hyper1`.

Lexicographic order happens to agree with numeric order for values R prints in
full, so `"0.05" < "0.6"` is correct. It disagrees for values small enough that R
switches to scientific notation:

```r
as.character(2.5e-05)          # "2.5e-05"
"2.5e-05" > "0.6"              # TRUE   <- treated as a gain
as.numeric("2.5e-05") > 0.6    # FALSE
```

The crossover is at about `1e-4`. So a chromosome carrying a gain of under
roughly 10–25 kb — an ordinary focal event — was scored as gained along its
entire length. Two such chromosomes among 3, 5, 7, 9, 11, 15, 19, 21 were enough
to call a sample hyperdiploid.

**Effect:** false-positive HRD calls, in samples whose only gain on an HRD
chromosome is very small. How many samples this touches depends entirely on the
segmentation: callers that emit small gained segments are affected more.
Downstream, an HRD call routes a sample into `HRD_RAS`, `HRD_Gains` or
`HRD_Complex_Cytogenetic` instead of the group its other features imply.

**Now:** `.mmgc_gt()` compares numerically.

---

## 2. The additional-gains burden was always zero

**Changes results. Opt out with `MMGenoCla(LEGACY.ADD.GAINS = TRUE)`.**

`MMGenoCla()` computed the secondary-gain burden as

```r
all_gen$add_gains[all_gen$CNV_Amp_14q11.2 + all_gen$CNV_Amp_17q22 + ... >= 3] <- 1
```

over eight features. But `FinalMatrixGenerator()` deliberately excludes
`Amp_14q11.2` from the matrix ("Removing the GISTIC peaks Amp_14q11.2 and
Del_21p11.2"), so `all_gen$CNV_Amp_14q11.2` is `NULL`. In R, `NULL + numeric`
yields `numeric(0)`, the comparison yields `logical(0)`, and the subset assignment
does nothing at all — silently:

```r
d <- data.frame(a = 1:3); d$x <- 0
d$x[d$a + d$missing >= 3] <- 1
d$x        # 0 0 0
```

**Effect:** `add_gains` is 0 for every sample, so the rule

```r
clusters[HRD == 1 & deletions == 0 & add_gains == 1 & ...] <- "B_HRD_gains"
```

never fires and the **`HRD_Gains` group is empty**. Those samples fall through to
the later catch-all and are labelled `HRD_RAS` instead.

**Important qualification.** This bites only when the matrix comes from the
`clinic10` column selection that is *active* in the published script. The
published feature matrix itself **does** carry `CNV_Amp_14q11.2` (and
`CNV_Del_21p11.2`): it was built with the alternative selection that sits
commented out just above, which keeps both peaks. So `add_gains` did work in the
published analysis, and `HRD_Gains` was populated there. The bug affects anyone
running `FinalMatrixGenerator()` as shipped and then classifying — which is the
normal path, and why it is fixed here.

**Now:** the burden is summed over whichever of the eight features are present —
seven, in practice — keeping the `>= 3` threshold, and a message reports which
were used. The threshold was calibrated on eight features, so with seven it is
marginally easier to reach; that is a deliberate trade against the feature being
inert.

---

## 3. Hyperdiploidy required the minor allele; the published matrix did not

**Changes trisomy calls. Opt out with `HYidentification(REQUIRE.MINOR.ALLELE = TRUE)`.**


The published function selects gained segments with

```r
cnv.gain = cnv.n[cnv.n$total >= 2.5 & cnv.n$minor >= 0.5, ]
```

requiring the second allele to be retained, so a chromosome at three copies with
the minor allele lost (`3+0`) is not a trisomy.

Checked against the published matrix for the paper's own cohort, that is not how
the matrix was made:

| Gain definition | Trisomy calls | Agreement with the matrix |
|---|---|---|
| `total >= 2.5 & minor >= 0.5` (as published) | 4,878 | 96.49% — misses 535 |
| `total >= 2.5` (total copy number only) | 5,413 | **99.99%** |

The matrix holds exactly 5,413 positive trisomy calls. Every one of the 535 the
allele-aware rule misses is either a `3+0` gain or a sample with no
allele-specific data — `MMRF_1719`, for instance, is `total = 3, minor = 0` across
all of chromosome 9.

It matters beyond reproduction: allele-specific copy number is absent for whole
cohorts. All 159 Moffitt samples in that dataset have `minor = NA` on every
segment, and the published filter silently reports every one of them as diploid.

**Now:** `REQUIRE.MINOR.ALLELE` makes the choice explicit, and **defaults to
`TRUE`, the published behaviour** — so this function reproduces the published code
unless you ask otherwise. It warns when samples carry no minor-allele values at
all, which is the case that would otherwise pass silently (159 of 1,693 samples on
the published cohort).

Set `REQUIRE.MINOR.ALLELE = "non-HRD"` to reproduce the published *matrix*: total
copy number alone decides a trisomy, while a gain anywhere else must still retain
its minor allele. That also keeps the whole-genome-doubling call strict, since
doubling means gaining both alleles. `FALSE` drops the requirement everywhere and
gives an identical feature matrix to `"non-HRD"` (222,592 of 222,592 cells), so
`"non-HRD"` is the better of the two.

The equivalent gene-level rule inside `MutationsMatrixGenerator()` is left as
published — see [VALIDATION.md](VALIDATION.md#residual-disagreement).

## 4. hg38 mutational signatures were fitted against hg19 contexts

**Changes results for hg38 and mm10. No opt-out — the old path cannot run.**

`mmsig::mm_fit_signatures()` hardcodes `BSgenome.Hsapiens.UCSC.hg19` for the
trinucleotide context lookup, and `mmsig:::getStrandBias()` hardcodes the hg19
TxDb. The published hg38 script worked around this by sourcing a privately
patched copy:

```r
source("~/Google Drive/.../mmsig-master/R/mmsig_functions_hg38.R")
```

which fails on any other machine and leaves no record of what was patched.

Fitting hg38 coordinates against hg19 contexts is wrong: the flanking bases at a
given coordinate differ between builds. Measured on 368 substitutions on chr1,
**82 of the 96 mutation classes come out with different counts, and about 28% of
contexts are misassigned.**

**Now:** `R/09-APOBECanalyzer.R` reimplements only the two genome-dependent steps
with the genome as an argument, calling mmsig's own internals for the fitting
itself. hg19, hg38 and mm10 all work with no patched copy.

Note this means hg38 APOBEC values from this repository will differ from any
produced with unpatched mmsig — correctly so.

---

## 5. Genes with a blank role were dropped from every matrix

**Affected mm10 only.**

Driver role was tested as

```r
ukn = drivers2[drivers2$Role == "Unknown",]; ukn1 = drivers2[is.na(drivers2$Role),]
```

`Driver_summary_hg19.txt` and `Driver_summary_hg38.txt` write unassigned roles as
`NA`, but `Driver_summary_mm10.txt` writes an empty string, which matches neither
test. Nine mouse mutation features (`DDX3X`, `INO80E`, `MRGPRX4`, `NDUFC2`,
`PANK3`, `SLC35G5`, `SYT17`, `TNFRSF11B`, `USP7`) were silently absent.

**Now:** blank, `NA` and the literal `"NA"` are all treated as unknown role.

---

## 6. Smaller fixes

None of these changed published results — they are latent failures that a
different cohort would hit.

**`X` and `Y` chromosome lengths were hardcoded to hg19 values.**
`HYidentification()` appended `c(23, 155270560)` and `c(24, 59373566)` to whatever
size file was passed, including the hg38 one, even though both files already carry
`X` and `Y` rows. Only per-chromosome gain fractions for X/Y were affected and
those feed no HRD call. Now read from the file.

**A logical index carrying `NA` inserted phantom all-`NA` rows.** Introduced by
this refactor and caught by running on real data: `x[cond, ]` where `cond` has
`NA` does not skip those rows, it adds an all-`NA` row for each. With 22,820
segments having an unknown minor allele, the gain filter produced 3,767 phantom
rows and `GRanges()` then rejected the `NA` coordinates. The published code hit
the same thing and removed the phantom rows afterwards with `complete.cases()`.
Every logical row-subset in the pipeline now selects through `which()`.

**`-grep()` returned a zero-column table when the pattern matched nothing.**
`FinalMatrixGenerator()` used `clinic7[, -grep("IGLL5", colnames(clinic7))]`. With
no match, `grep()` returns `integer(0)`, and `x[, -integer(0)]` is *every column
dropped*, not none — so a cohort that never called IGLL5 lost its whole matrix.
Now guarded by `.mmgc_drop_cols()`.

**Tumour-suppressor combining assumed exactly 36 genes.** The mutation and
deletion matrices were merged by hardcoded column positions (`1:36`, `37:72`,
`73:108`). Any driver list with a different number of tumour suppressors produced
silently wrong output. Now paired by gene symbol.

**The final matrix was assembled by hardcoded column indices.**
`clinic7[, c(1, 2:45, 68:99, 46:67, 100:132, 135:138, 134)]` assumed exactly 138
columns. Its only effect was to drop the raw `apobec` column and reorder, and the
next step re-selected by name, so it is now a named drop.

**Missing peaks aborted with an unhelpful error.** Peak columns were selected by
name directly, so one absent peak gave "undefined columns selected" with no
indication which. Now reported by name and carried as `NA`.

**Only one missing peak was ever added back.** `FocalLargeMatrix()` had
`for (jjj in length(mis.peak))` instead of `seq_along(...)`, so when several peaks
were untouched by a cohort only the last got a column.

**Synonymous DIS3 mutations counted as hotspots.** `dis3_hs()` labelled
synonymous rows `"Synonymous"`, then set `cristal_code <- 1` for everything not
equal to `0` — including them. Harmless as called (the input was already filtered
to non-synonymous) but wrong if handed a full `annotmuts`. Now filtered
explicitly.

**A gene covered only by segments of unknown copy number became `-Inf`.** The
per-gene call took `max()` over the covering segments without `na.rm`, so one
missing value poisoned the gene; the reverse case now returns `NA`.

**`OUT.DIR` had to end in `/`.** Some paths were pasted with a separator and
others without, so the pipeline only found its own intermediate files when
`OUT.DIR` carried a trailing slash. Now joined with `file.path()` throughout, so
either form works.

**`SegmentCollapsing()` paired run starts with run ends positionally.** It split
each run's first and last row into two tables and zipped them by row number,
which is only correct while the input is sorted. Now grouped by run id.

---

## 7. Things left as they are

**`CYTOBAND.FILE` and `GISTIC.ARMS` do nothing.** `HYidentification()` reads both
into variables nothing else touches, and `FocalLargeMatrix()` reads `CYTOBAND.FILE`
and `METADATA` without using either. They remain in the signatures so existing
call sites keep working, and are documented as ignored.

**The `HH_CCND1` rule is unreachable.** It tests `clusters == 0` one line after
every remaining `0` has been set to `RR_Simple`. It is not in the published group
list either. Omitted.

**`MMGenoCla()` depends on feature order for the deletion burden.** The published
code took the span from the first column matching `Del`/`loss`/`CNV.SNV_` to the
last and summed everything between, which is the intended 50 features only
because `FinalMatrixGenerator()` emits them contiguously. Now selected by name —
same features, no longer order-dependent.

**The interactive prompts are off by default.** `MMGenoCla()` blocked on two
`readline()` calls, which is why the published run script has bare `y` lines after
the call. Pass `CONFIRM = TRUE` to restore them.
