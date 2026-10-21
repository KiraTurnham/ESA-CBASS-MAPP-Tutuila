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

Some scripts write figures to a `Plots/` directory. Create that directory if it is not present before running a script that saves plots there.

## Prerequisites

- R (the repository does not currently specify a tested R version).
- RStudio is recommended; some scripts use RStudio file-selection and viewing functions.
- R packages used across the analyses include `tidyverse`, `here`, `vegan`, `CBASSED50`, `lme4`, `lmerTest`, `emmeans`, `pvclust`, `pheatmap`, `caret`, and `pairwiseAdonis`, among others. Individual scripts load their own dependencies.

There is no package lockfile or centralized dependency-installation script in this repository. Install the required packages in your R environment before running each analysis. `CBASSED50` is an external package and may require installation from its maintained source rather than the standard CRAN repository; follow that package's installation instructions.

## Getting started

1. Clone or download this repository and open the repository root as your RStudio project or working directory. Scripts use `here()` to locate files relative to the project root.
2. Install the R packages required by the script(s) you intend to run.
3. Run the analyses in dependency order, checking each script's input paths and comments before execution. A typical workflow is:
	- [Scripts/01_ED50_calc.R](Scripts/01_ED50_calc.R) — calculate CBASS ED values from a selected run file. This script opens a file chooser and is intended for interactive use.
	- [Scripts/01b_color_calc.R](Scripts/01b_color_calc.R) — combine color-score and PAM data and estimate color-slope responses.
	- [Scripts/03_DataCompiler_Outliers.Rmd](Scripts/03_DataCompiler_Outliers.Rmd) — compile MAPP phenotype data, apply the documented filtering and preprocessing, and assemble analysis metadata.
	- [Scripts/04_MAPP_Multivariate_filtered_2_outliers.R](Scripts/04_MAPP_Multivariate_filtered_2_outliers.R) — run multivariate analyses, including PERMANOVA, dispersion tests, nMDS, and trait-vector summaries.
	- [Scripts/05_model_comparisons.R](Scripts/05_model_comparisons.R) — compare species- and cluster-based mixed-effects models for thermal and other response metrics.
	- Other scripts in [Scripts/](Scripts/) generate figures, evaluate alternative CBASS temperature ramps, or perform supplementary analyses.

Several downstream scripts load precomputed files from [Outputs/](Outputs/). Run the relevant upstream analyses first, or confirm that the expected input files are already present. File names and filtering choices are encoded in individual scripts and may need to be updated when using a different dataset or study iteration.

## Reproducibility notes

- Results depend on the input data, filtering rules, model settings, and package versions. Check script comments and recorded output files when reproducing a particular analysis.
- Some analyses use random procedures; where set, the random seed is specified in the script.
- The scripts are research analyses rather than a single automated pipeline; they may require interactive steps and manual review.
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
