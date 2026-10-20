# mm10 references for the MM genomic classification (vkMYC mouse model)

These files are the **mouse (mm10 / GRCm38)** equivalents of the human (hg38)
references in `../references/`. They are produced by
[`../script/build_mm10_references.R`](../script/build_mm10_references.R).

## Method — gene-based, NOT coordinate liftover

Synteny is broken between human and mouse, so GISTIC peaks are **relocated by
orthologous gene**, not by lifting coordinates. Each human recurrent lesion is
anchored on its target gene's **mouse ortholog** (via the `babelgene` package),
and the region is placed at that ortholog's mm10 position (± 1 Mb flank). The
human peak name / `code` is kept unchanged, so downstream classification
categories stay identical and comparable (e.g. `Amp_1q21.1` still means "MCL1
gain" — but is now detected on **mouse chr3**, where *Mcl1* lives).

Anchoring priority per peak (`anchor_level` column):
1. `curated`        — from the editable `peak_target_override.txt` (well-known MM
                      targets; the GISTIC *peaks* are narrow and often do not
                      physically contain their target, e.g. MCL1 falls between
                      the two 1q21 peaks, TP53 sits just outside Del_17p13.1).
2. `MM_driver_near` / `driver_near` — nearest MM driver gene within 2 Mb.
3. `all_genes`      — fallback: the peak's genes, majority mouse chromosome.

## Files
| File | Purpose |
|---|---|
| `Genomic_Coordinates_GISTICpeaks_mm10.txt` | GISTIC peaks relocated to mm10 (+ anchor annotation columns) |
| `peak_target_override.txt` | **editable** curated peak→target-gene map; edit & re-run |
| `Driver_summary_mm10.txt`  | MM drivers mapped to mouse orthologs (`gene_name`=mouse, `human_gene_name` kept) |
| `UCSC_mm10_ref_genes.dms`  | mm10 gene models (UCSC refGene) |
| `chrom_sizes_mm10.txt`     | mm10 chromosome sizes |
| `cytoBand_mm10.txt`        | mm10 cytobands |
| `QC_unmapped_peaks.txt`    | human peaks with no confident mouse anchor (mostly IG loci) |
| `QC_split_peaks.txt`       | peaks whose genes scatter across ≥2 mouse chromosomes |
| `QC_unmapped_drivers.txt`  | driver genes with no mouse ortholog |
| `raw/`                     | downloaded UCSC source tables |

## To review / refine
- **Check the `all_genes` peaks** in the GISTIC file (`anchor_level` column) —
  their anchor is automatic and may not be the true CNV target. To fix one, add
  a row to `peak_target_override.txt` (`code <TAB> human_gene_symbol`) and re-run
  `build_mm10_references.R` (your edits are preserved).
- **Unmapped peaks** are mostly immunoglobulin loci (IGH 14q32.33, IGK 2p11.2,
  15q11.1) with no meaningful mouse focal target — usually fine to leave out.

## Regenerate
```bash
Rscript ../script/build_mm10_references.R
```
Needs: R with `data.table`, `babelgene`, `GenomicRanges`, and internet (to pull
UCSC mm10 refGene / cytoBand / chrom.sizes into `raw/`).
