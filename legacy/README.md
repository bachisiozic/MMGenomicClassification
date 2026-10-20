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
