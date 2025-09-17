# sim_ancestry_GWAS
Code for the following paper:

[Yang, Y. & Lawson, D.J. From individuals to ancestries: towards attributing trait variation to haplotypes. PLoS Genet (2025) to appear. Available from medRxiv (2025). doi:10.1101/2024.03.13.24304206](https://www.medrxiv.org/content/10.1101/2025.03.13.25323895v1)

The file `sim_Tractor.R` is for simulating ancestry-specific GWAS.

The file `PCA_UKB_analysis_xgboost.R` performs xgboost on birth location within the UK.

The file `visxgboost.R` visualises these results.

Note that the UKB analysis is replicable if you have access to the UK Biobank dataset and generate your own individual-level PCs and HCs. Because consent can be withdrawn at any stage, it is likely that you will obtain a different set of training/test individuals and therefore some variation for outlying regions is expected.


