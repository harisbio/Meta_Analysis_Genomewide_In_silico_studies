# Meta-Analysis of Genome-Wide In Silico Studies

Code for the article *Twenty-Five Years of Genome-Wide Gene Family Identification Analysis in Plants: Trends, Tools,
and Methodological Shifts – A Systematic Review and Meta-Analysis* (Haris, Azeem, Golokhvast and Nawaz).

This repository contains the scripts that produce every figure and every value reported in the article and its
Supplementary File. It contains code only: the datasets and the manuscript are not included.

The program used to extract data from the included articles is published separately:
**gene-family-extractor** ([extraction repository URL]).

## Contents

```
├── config.py                         input and output locations
├── run_all.py                        regenerate all figures and values
├── analysis/
│   └── manuscript_values.py          every number reported in the article, recalculated from the dataset
├── figures/
│   ├── main/                         Figures 2-11 of the article
│   │   ├── fig2_country.py           country distribution (world map and ranked bars)
│   │   ├── fig3_trends_network.py    papers per year with fitted trend (negative binomial), cumulative papers
│   │   ├── fig4_team_size.py         authors per paper and trend (Poisson, robust standard errors)
│   │   ├── fig5_radial.py            journals and study design
│   │   ├── fig6_organisms.py         network of study and reference organisms
│   │   ├── fig7_families.py          method use and growth of the most studied families
│   │   ├── fig8_genes.py             genes identified per paper over time
│   │   ├── fig9_pipeline.py          genome sources, study design and tool use
│   │   ├── fig10_phylo_dup.py        phylogenetic grouping terms and duplication types
│   │   └── fig11_validation.py       RNA-seq and qPCR validation over time
│   └── supplementary/                Figures S1-S26 and the tables behind them
│       ├── global_analysis.py
│       ├── extra_figures.py
│       └── family_names.py           gene family name harmonization
├── data/README.md                    which data files the scripts expect
└── legacy/genomewide_pdf_local_extractor/
                                      the earlier rule-based extraction pipeline (see below)
```

## Data access

The datasets are not public. The 310-study dataset is available from the corresponding author on reasonable request.
The included studies are listed in Supplementary Table S7 of the article. See [data/README.md](data/README.md) for the
files the scripts expect and where to put them.

## Installation

Python 3.10 or later.

```bash
python -m pip install -r requirements.txt
```

The figures use the Arial font. Figure 2 and Figure S2 are drawn with Plotly and exported with Kaleido, which needs a
local Chrome or Chromium installation (`kaleido_get_chrome` installs one if needed).

## Usage

Place the data files in `data/` (or set the environment variables below), then:

```bash
# all main figures (2-11) and the reported values
python run_all.py

# also the Supplementary File figures and tables
python run_all.py --global

# a single figure
python figures/main/fig3_trends_network.py
```

Outputs are written to `output/`: one folder per figure with PNG (300 dpi), PDF and SVG files and the CSV tables
behind each figure, and `output/Supplementary/` for the supplementary figures and tables.

Paths can be changed without editing the code:

| Variable | Default | Used for |
|---|---|---|
| `META_ANALYSIS_DATA_FILE` | `data/Data_cleaned.csv` | main figures and reported values |
| `META_ANALYSIS_GLOBAL_DATA_FILE` | `data/abstract_dataset.csv` | Supplementary File figures |
| `META_ANALYSIS_OUTPUT_DIR` | `output/` | all outputs |

## Statistical models

- **Annual papers (Figure 3):** negative binomial regression of papers per year on year, with the dispersion
  estimated from a Poisson fit; growth is reported as the percentage change per year with a 95% confidence interval.
- **Team size (Figure 4):** Poisson regression of authors per paper on year, with robust (HC1) standard errors.
- **Genes per paper (Figure 8):** ordinary least squares regression of log10 genes per paper on year.

All other results are counts, proportions and medians. Unless stated otherwise, percentages use all 310 studies as
the denominator; where a variable was not reported by all studies, the denominator is the studies that reported it.

## Earlier extraction pipeline

`legacy/genomewide_pdf_local_extractor/` is the rule-based pipeline, with optional local language-model repair through
Ollama, that produced the first version of the dataset. It is kept here for transparency. The data reported in the
article were extracted with gene-family-extractor ([extraction repository URL]) and compared with this earlier
extraction field by field.

## Citation

Please cite the article. See [CITATION.cff](CITATION.cff).

## License

The code is released under the MIT License ([LICENSE](LICENSE)). The licence covers code only, not the datasets,
the manuscript or third-party content.
