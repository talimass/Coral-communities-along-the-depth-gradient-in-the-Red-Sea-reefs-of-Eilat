# Source of the reported PERMANOVA values 

Author line in the script header: Tamar Shifroni-Taylor, 04/04/2024. The file here is a copy with local folder paths removed (only the working-directory lines changed); the unchanged original is kept by the authors.

## What the script does (lines 779-787)
adonis2 on the Hellinger-transformed hard-coral counts, formula `DEPTH + SITE + DEPTH * SITE`, 1,000 permutations, with DEPTH as a factor. `adonis2` has no `distance` argument, so `distance = "jaccard"` in the call is ignored and the analysis uses Bray-Curtis dissimilarities (sums of squares identical to `method = "bray"`: 32.80974 in vegan 2.6-4, 57.13557 in vegan 2.7-5).

## Result of the re-run (464 images; log files in this folder)
R 4.3.3, vegan 2.6-4 (the versions named in the manuscript):

| Term | Df | SumOfSqs | R2 | F | Pr(>F) |
|---|---|---|---|---|---|
| DEPTH | 4 | 32.810 | 0.25835 | 52.716 | 0.000999 |
| SITE | 2 | 11.005 | 0.08665 | 35.363 | 0.000999 |
| DEPTH:SITE | 8 | 13.321 | 0.10489 | 10.702 | 0.000999 |
| Residual | 449 | 69.863 | 0.55011 | | |
| Total | 463 | 126.998 | 1.00000 | | |

These match the manuscript values (R2 = 0.26 depth, 0.09 site, total 0.45). The same R2 values were obtained in R 4.5.3 with vegan 2.7-5 (`R4.5.3_vegan2.7-5.log`).

## Repository script (`Bentic_script.R`), run log `repository_benthic_script_R_run.log`
`adonis2(bray ~ DEPTH)` on Bray-Curtis dissimilarities of raw hard-coral counts, 999 permutations: R2 = 0.2265, F = 33.60. The dispersion test by depth in that script (betadisper, same data and distance): F = 3.52, p = 0.0076. This dispersion result belongs to that depth-only model only. No dispersion test was run for the analysis above.

## Not tested
Pairwise PERMANOVA, the NMDS distance setting, and which LDM version produced Table S7.
