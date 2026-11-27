# Eilat paper package

Code, data and outputs for the Eilat coral depth-gradient paper, collected in one folder. Everything here is a copy; no other file in the repository was changed, moved or deleted.

## Contents
- `Scripts/`
  - `Bentic_script.R` and `FIRe_boxplots_PCA.R`: the repository scripts, unchanged.
  - `Older_benthic_script_Shifroni-Taylor_2024-04-04_paths_removed.R`: the older script that produced the reported PERMANOVA (lines 779-787). Only the working-directory lines were changed.
- `Data/`: data files that the scripts read or write through explicit file calls in the code (all inputs found in those calls are included; dynamically built file names were not checked).
- `Data/provenance_uncertain/`: data files whose names were not found in literal references in the three scripts (this is a text match, not proof that no script uses them) (`FIRe_all.csv`, `FIRe_averaged.csv`, `formatted_fire_data.csv`, `HC_LDM_res_t_with_RA.csv`, `LDM_HC_results.csv`, `Hard_corals_totals_prepare_2plot.csv`). Their role in the paper is not yet confirmed by the authors, so they are kept, not judged.
- `Output/figures/`, `Output/stats_tables/`: outputs whose file names appear literally in the scripts.
- `Output/provenance_uncertain/`: outputs whose names were not found literally in the scripts (an earlier ART-ANOVA workbook, a combined plot, an `Alga cover` figure, an LDM workbook, a draft PCA slide file). Kept for the authors to decide.
- `Reproducibility/`: re-runs of the PERMANOVA lines and of `Bentic_script.R`, with logs and notes (`README_reproducibility.md`).

## Using the scripts
- `Bentic_script.R` and `FIRe_boxplots_PCA.R` contain absolute paths to the author's computer (`setwd(...)` and output folders). Before running, change these to the `Data/` and `Output/` folders here. The scripts were not edited.
- The older script needs `setwd` set to `Data/`.
- Every input file found in the scripts' explicit read calls is present in `Data/`; this is not a proof of all dynamic dependencies.

## Still open (authors)
- Which FIRe file feeds which analysis: `FIRe_averaged_122025.csv` (statistics) and `FIRe_averaged_122025_1.csv` (PCA) are the files the script reads. The PCA file has two colonies fewer (Pocillopora 93, Porites 117); the reason is not yet known.
- Which LDM version corresponds to Table S7.
- The Fig. 1 map has no source file here.
- Only the PERMANOVA lines and the repository benthic script were re-run; other analyses were not independently reproduced.
