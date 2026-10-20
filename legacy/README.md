# Original scripts, as published

These three files are the code behind
[doi:10.1200/JCO.23.01277](https://doi.org/10.1200/JCO.23.01277), **unmodified**.
They are kept so that published results stay reproducible and so the reorganised
version in `../R/` can be diffed against its origin.

| File | What it is |
|---|---|
| `JCO_GenomicClassification.R` | the hg19 function file |
| `JCO_GenomicClassification_hg38.R` | the hg38 function file |
| `RUN_JCO_Genomic_Classification.R` | the driver script, with the authors' own paths |

## Do not run these directly

They will not work outside the machine they were written on:

- `JCO_GenomicClassification_hg38.R:1551` sources a privately patched copy of
  mmsig from an absolute Google Drive path.
- `RUN_JCO_Genomic_Classification.R` hardcodes absolute paths to reference files,
  input data and output folders throughout.
- Both open with 25 `library()` calls, ten of them for packages no function uses.

Use `../R/` instead, via `../MMGenomicClassification.R`.

## They also contain two defects that changed results

Both are fixed in `../R/`, with an opt-out to reproduce the behaviour here:

- Hyperdiploidy was called by comparing gain fractions as text, so a gain under
  about 1e-4 of a chromosome scored as a whole-chromosome gain.
- The additional-gains burden summed a feature the pipeline had already dropped,
  so it was always 0 and the `HRD_Gains` group could never be assigned.

`../docs/KNOWN_ISSUES.md` has the detail, including how to switch each one back.
