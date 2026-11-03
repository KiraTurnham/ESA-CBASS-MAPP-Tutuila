# CBASS and MAPP of American Samoa Corals

## Overview

This repository contains data and analysis scripts for investigating thermal susceptibility and photophysiological variation in four common forereef corals from Tutuila, American Samoa: *Acropora globiceps*, *Isopora crateriformis*, *Acropora abrotanoides*, and *Acropora hyacinthus*.

The study combines acute thermal-stress experiments using the Coral Bleaching Automated Stress System (CBASS) with Multispectral Aquatic Phenotyping Platform (MAPP) measurements. It asks:

1. How is baseline photophysiological trait variation partitioned among coral holobionts?
2. Do host species and/or MAPP-derived photophysiological phenotypes differ in thermal tolerance, measured using ED50 and color-slope paling responses?
3. Which photophysiological metrics are associated with different phenotypes?

The analysis estimates thermal-performance metrics and color-loss slopes, compiles and filters MAPP measurements, explores multivariate patterns using ordination and PERMANOVA, and identifies functional clusters using model-based clustering and AIC/BIC comparisons. Mixed-effects models then compare thermal-performance measures across species and functional clusters. Trait vectors are projected onto baseline ordinations to identify metrics associated with species and cluster patterns.

## Repository contents

- [Data/](Data/) — CBASS/PAM measurements, MAPP phenotyping data and metadata, lookup tables, and other study inputs.
- [Scripts/](Scripts/) — R scripts and R Markdown analyses for thermal-response calculations, data compilation, multivariate analyses, clustering, model comparisons, and figures.
- [Outputs/](Outputs/) — processed datasets, model results, and analysis products used by downstream scripts.

Analysis scripts write intermediate data to `Data/` and `Outputs/`, and figures to `Plots/`. These directories are included in the repository.

## Prerequisites

- R and Pandoc (for knitting R Markdown documents).
- Open the repository root as the RStudio project/working directory. Scripts use `here()` for project-relative paths.
- Install the packages loaded by the script you plan to run. The workflows use packages including `tidyverse`, `here`, `vegan`, `CBASSED50`, `lme4`, `lmerTest`, `emmeans`, `pvclust`, `pheatmap`, `caret`, `pairwiseAdonis`, `ggtree`, and `ggnewscale`.
- `CBASSED50` and `ggtree` may require installation from their project/Bioconductor sources. There is not yet a package lockfile, so package versions are not pinned.

## Getting started

1. Run [Scripts/00_Temp_dataset_comparison.R](Scripts/00_Temp_dataset_comparison.R) to fit ED values from the eight named raw CBASS/PAM runs. It writes `Data/EDsdf_all_ramps.csv` and refreshes `Data/EDsdf-all.csv`, the ED table used by the MAPP compiler.
2. Run [Scripts/01b_color_calc.R](Scripts/01b_color_calc.R) to join color-score and PAM data and write `Outputs/color_slope.csv`.
3. Knit [Scripts/03_MAPP_DataCompiler.Rmd](Scripts/03_MAPP_DataCompiler.Rmd). It requires the ED and color-slope files from steps 1–2, collection/MAPP metadata, the metric lookup table, and processed CSVs under `Data/phenotyping_data/normalized/`. It writes `Outputs/photophys_preprocessed_data2.RData`, including the filtered metric matrix, metadata, clusters, clustering fit, and metric-family annotations.
4. Run [Scripts/04_MAPP_Multivariate_NMDS.R](Scripts/04_MAPP_Multivariate_NMDS.R) to run dispersion/PERMANOVA tests and nMDS/envfit. It creates `Outputs/ef_filtered2.rds`, `Outputs/top_vectors2.rds`, and the top-vector table.
5. Run [Scripts/05_MAPP_Photophys_traces_script_filtered_2_outliers.R](Scripts/05_MAPP_Photophys_traces_script_filtered_2_outliers.R), [Scripts/PhenotypeHeatmap_FILTERED.Rmd](Scripts/PhenotypeHeatmap_FILTERED.Rmd), or [Scripts/06_Dendrogram.Rmd](Scripts/06_Dendrogram.Rmd) for downstream MAPP plots. The trace script needs the top-vector file from step 4.

### CBASS response figures

- [Scripts/02_Taxa_differences_in_EDs.R](Scripts/02_Taxa_differences_in_EDs.R) uses `Data/EDsdf_all_ramps.csv` from step 1 and writes `Outputs/plot_ed50.rds`.
- [Scripts/06a_Plot_ED_emmeans.R](Scripts/06a_Plot_ED_emmeans.R) and [Scripts/06b_Plot_Color_emmeans.R](Scripts/06b_Plot_Color_emmeans.R) use the compiled data from step 3. They create the ED50 and color-slope cluster/species plot objects.
- After those scripts and step 2, run [Scripts/01c_combined_cbass_response_plot.R](Scripts/01c_combined_cbass_response_plot.R) to combine the saved plots.

### Other entry points

- [Scripts/01_ED50_calc.R](Scripts/01_ED50_calc.R) is an interactive per-run tool: select a raw CBASS file, and it writes a run-specific ED table and summary beside that file. It is not required for the all-ramp workflow above.
- [Scripts/06_model_comparisons.R](Scripts/06_model_comparisons.R) runs additional model comparisons and produces supplementary tables/figures from the compiled data.
- [Scripts/Fig1_TUT_map.R](Scripts/Fig1_TUT_map.R), [Scripts/S01a_Rarefaction_SE.R](Scripts/S01a_Rarefaction_SE.R), and [Scripts/SFig_dose_response_curves.R](Scripts/SFig_dose_response_curves.R) are supplementary analyses; see their inputs and outputs in the scripts.
- Script names with `2` or `filtered` identify the selected MAPP data/filtering iteration. They are retained to avoid breaking saved output references.

`Scripts/00_visualization_prep.R` defines shared palettes and plotting helpers only; it does not load analysis data. Scripts load their input `.RData` explicitly. `prepare_site_metadata()` applies the site labels and coordinates where needed.

## Reproducibility notes

- The main script order above follows actual file and object dependencies; standalone supplementary scripts are not one linear pipeline.
- Seeds are set for randomized analyses and clustering is run sequentially in the MAPP compiler. Model and permutation settings remain in the scripts.
- `Data/EDsdf-all.csv` is regenerated from the all-ramp fit rather than treated as a manually maintained prerequisite.
- `Outputs/photophys_preprocessed_data.RData` is a legacy saved snapshot; current downstream scripts use `photophys_preprocessed_data2.RData` from the compiler.
- Use `sessionInfo()` to record the R/package versions used for a run. Cross-machine reproducibility is limited until a package lockfile is added.
- Do not overwrite or redistribute source data without confirming the applicable data-use terms.

## Resources

- CBASS: Coral Bleaching Automated Stress System.
- MAPP: Multispectral Aquatic Phenotyping Platform.
- Analysis entry points and data-specific details are documented in the scripts under [Scripts/](Scripts/).

## Version Control Platform
- Git

## License

## Disclaimer
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
