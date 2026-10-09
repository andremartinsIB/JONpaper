# Tomato root metabolome under *Metarhizium anisopliae* and *Meloidogyne incognita*: factorial reanalysis

This repository contains the data and the R script used to analyse the tomato root metabolome in a 2 × 2 factorial experiment with the entomopathogenic fungus *Metarhizium anisopliae* (MA) and the root-knot nematode *Meloidogyne incognita* (MI). The code is the version run for the associated manuscript.

**Associated manuscript:** Metabolic responses of tomato roots to Meloidogyne incognita infection and Metarhizium anisopliae application.

## Experimental design

- Factors: MA (absent/present) × MI (absent/present), completely randomised design.
- Six plants per treatment, 24 plants in total.
- Treatments and sample numbers in the data matrix:

| Treatment | Code in data file | Samples |
|---|---|---|
| Control | `Teste` | 1–6 |
| *M. anisopliae* | `MA` | 7–12 |
| *M. incognita* | `MI` | 13–18 |
| *M. anisopliae* + *M. incognita* | `MA_MI` | 19–24 |

## Repository structure

```
├── CompleteAnalysis.R           complete analysis (statistics, tables and figures)
├── MetabolomicData.csv          GC–MS peak table (relative abundance)
├── Metadata.csv                 sample sheet with the MA and MI factors
└── resultados_FDR010/
    └── sessionInfo.txt          R session information from the run of the script
```

### Input format

`MetabolomicData.csv` is a semicolon-separated text file with a decimal point:

- row 1: sample number (1–24);
- row 2: treatment code (`Teste`, `MA`, `MI`, `MA_MI`);
- rows 3–87: one metabolite per row (85 metabolites), relative abundance per sample;
- `NA`: value not reported in the peak table.

The script reads the treatment of each sample from row 2 of this file. `Metadata.csv` is provided as the sample sheet of the experiment and is not read by the script.

## Analyses

All parameters are set in section 0 of `CompleteAnalysis.R`.

1. Metabolite selection: of the 85 metabolites, 49 have values in at least 3 of 6 plants in every treatment; 3 of these are on the script's list of probable non-biological compounds and are excluded, leaving 46 metabolites tested.
2. Two-way ANOVA per metabolite on log2-transformed values (`lm(y ~ MA * MI)` with sum-to-zero contrasts, type III sums of squares with `car::Anova`), missing values omitted, with a Shapiro–Wilk test on the residuals and a Brown–Forsythe test for homogeneity of variance.
3. limma with a cell-means model, factorial contrasts and empirical Bayes moderation; simple-effect contrasts are also reported.
4. Two-way PERMANOVA on the global profile (Euclidean distance; log2 transformation, mean centring and Pareto scaling; missing values replaced by half of the metabolite minimum in this step only), checked against `vegan::adonis2`, and PERMDISP.
5. Exact pairwise PERMANOVA for the four simple-effect comparisons.
6. Distances between treatment centroids with bootstrap intervals.
7. PERMANOVA by chemical class, with classes defined in the script before testing.
8. Sensitivity analysis without MA plants 10 and 12.

Benjamini–Hochberg false discovery rate (FDR) correction is applied within each model term or contrast. Results are reported at FDR ≤ 0.10 (`FDR_LIMIAR`), and results with FDR ≤ 0.05 (`FDR_ESTRITO`) are flagged separately in tables and figures.

## How to run

1. Line 1 of the script (`setwd(...)`) sets the working directory on the authors' computer. Replace that path with the folder that contains this repository, or remove the line and set the working directory to that folder.
2. Run the whole script, for example with `Rscript CompleteAnalysis.R` from that folder, or with Source in RStudio.

The script creates the folder `resultados_FDR010/` with a spreadsheet of all results (`resultados_completos.xlsx`), figures in vector PDF and 600 dpi PNG (`figuras/`) and `sessionInfo.txt`. The run takes a few minutes because of the permutation and bootstrap steps.

## Software

The analysis was run with R 4.5.3 (conda installation) on Ubuntu 24.04.5 LTS, using:

| Package | Version | Source |
|---|---|---|
| car | 3.1-5 | CRAN |
| carData | 3.0-6 | CRAN |
| vegan | 2.7-6 | CRAN |
| permute | 0.9-10 | CRAN |
| ggplot2 | 4.0.3 | CRAN |
| ggrepel | 0.9.8 | CRAN |
| scales | 1.4.0 | CRAN |
| writexl | 2.0.1 | CRAN |
| limma | 3.66.0 | Bioconductor |

The full session information is in `resultados_FDR010/sessionInfo.txt`.

## Reproducibility

- Random seed: 20260930.
- Permutations: 9,999 (PERMANOVA, PERMDISP); bootstrap resamples: 9,999 (centroid distances).
- Exact pairwise PERMANOVA P values enumerate all possible labellings and do not depend on the seed.
- Section 9 of the script compares the deterministic results with reference values from an independent Python implementation of the same analysis (not included in this repository). The `VALIDACAO` sheet of the output spreadsheet should show `TRUE` in the `confere` column for every row.

Code comments and figure labels are in English; variable names, console messages, spreadsheet sheet and column names, and output folder and file names are in Portuguese BR.

## Licence

No licence has been assigned yet. Until a licence is added, all rights are reserved by the authors.
