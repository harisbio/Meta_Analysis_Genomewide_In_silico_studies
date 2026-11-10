# Data

The datasets are not distributed with this repository. They are available from the corresponding author on
reasonable request.

## `Data_cleaned.csv`: the 310 included studies

Used by `figures/main/` and `analysis/manuscript_values.py`. UTF-8 CSV, one row per study, with these columns:

`Paper ID`, `Country of Origin`, `Title`, `Authors`, `Year`, `Journal`, `DOI`, `Keywords`, `Primary Organism`,
`Secondary Organism`, `Gene Family`, `Genome Source`, `HMMER (with Pfam domains ID)`, `BLAST (E-value cut-off used)`,
`SMART`, `Other (Tools for Domain Analysis)`, `CDD`, `MEGA`, `IQ-TREE`, `MEME`, `Others (Tools for Motif Analysis)`,
`MCScanX`, `Other (Tools for Synteny)`, `TBTools`, `Other (Tools for Visualization)`, `Experimental Validation (qPCR)`,
`Other Tools Used In Methodology`, `Number of Genes Identified`, `Phylogenetic Structure`,
`Gene Duplication Events (TANDEM/SEGMENTAL)`, `Expression Analysis (RNA-seq)`, `Expression (microarray/EST)`,
`Genes per species`.

Field definitions and coding rules are in the codebook of the extraction program (gene-family-extractor).

## `abstract_dataset.csv`: the global abstract dataset

Used by `figures/supplementary/`. One row per paper from the broader search described in the Supplementary File;
only rows with `in_analysis == True` are analysed. The scripts use the bibliographic fields (`year`, `journal`,
`doi`), the abstract-extracted fields (`species_norm`, `n_species_norm`, `gene_family`, `family_class`, `analyses`,
`expression`, `validation`, `conditions`, `genes_total`, `groups`, `duplication_norm`, `selection`) and the OpenAlex
fields (`corr_country`, `corr_country_code`, `all_countries`, `international`, `cited_by`, `fwci`, `top10pct`).
