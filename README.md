# covid-brain-aging-thesis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Data: GSE188847](https://img.shields.io/badge/Data-GSE188847-informational.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188847)
[![R: ≥4.3](https://img.shields.io/badge/R-%E2%89%A54.3-blue.svg)](https://cran.r-project.org/)

Analytical pipeline for the manuscript:

> **Cellular composition shifts dominate the COVID-19 brain transcriptomic
> response: NF-κB, IFN/JAK-STAT and REST master regulators converge on a
> five-cell-type neurovascular signature in postmortem frontal cortex.**
> Elechi KW,. *Manuscript in preparation*, 2026.

This repository contains the complete R pipeline that reproduces the 13 figures
and headline statistics reported in the manuscript. The analysis re-examines
the public GSE188847 frontal cortex bulk RNA-seq cohort (Mavrikaki et al.,
*Nature Aging* 2022) using master regulator inference (VIPER + DoRothEA),
cell-type deconvolution (BRETIGEA), and composition-corrected differential
expression.

---

## What's in this repo

```
covid-brain-aging-thesis/
├── README.md                       (this file)
├── LICENSE                         (MIT)
├── CITATION.cff                    (GitHub "Cite this repo" metadata)
├── config.R                        (path configuration; edit before running)
├── scripts/                        (the 9 analysis scripts)
│   ├── step1_import_qc_pca_n54_v2.R
│   ├── step2_deseq2_n45.R
│   ├── step3_gsea_n45.R
│   ├── step4_viper_n45.R
│   ├── step5_bretigea_composition_n45.R
│   ├── step6_heatmap_regulons.R
│   ├── step7_clean_figures.R
│   ├── step7b_regulons_and_composition.R
│   └── step7d_volcano_labels_v2.R
├── data/
│   ├── README.md                   (data acquisition instructions)
│   ├── GSE188847_metadata_clean.csv
│   └── tx2gene.csv                 (Ensembl v104 tx → gene mapping)
└── docs/
    └── pipeline_overview.md        (one-page walkthrough of each script)
```

---

## Quick-start reproduction

### 1. Clone the repository

```bash
git clone https://github.com/kelyto30/covid-brain-aging-thesis.git
cd covid-brain-aging-thesis
```

### 2. Download GSE188847 quants

The 54 Salmon `.gz` quantification files are too large to commit to Git
and must be downloaded separately. See [`data/README.md`](data/README.md)
for instructions. Approximately 1–2 GB of disk space required.

### 3. Edit `config.R`

Open `config.R` at the repo root and set two paths:

```r
PROJECT_ROOT <- "/path/to/your/clone/of/covid-brain-aging-thesis"
QUANTS_DIR   <- "/path/to/folder/containing/GSE188847/.gz/files"
```

The other paths (metadata CSV, tx2gene CSV) are derived from `PROJECT_ROOT`
and require no editing.

### 4. Install R dependencies

The pipeline requires R ≥ 4.3 and a number of CRAN and Bioconductor
packages. From an R session, run:

```r
# CRAN
install.packages(c("ggplot2", "ggrepel", "dplyr", "tidyr",
                   "pheatmap", "tibble"))

# Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("tximport", "DESeq2", "apeglm",
                       "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
                       "viper", "dorothea", "BRETIGEA"))
```

### 5. Run the pipeline

From the repo root, run scripts in order:

```bash
Rscript scripts/step1_import_qc_pca_n54_v2.R
Rscript scripts/step2_deseq2_n45.R
Rscript scripts/step3_gsea_n45.R
Rscript scripts/step4_viper_n45.R
Rscript scripts/step5_bretigea_composition_n45.R
Rscript scripts/step6_heatmap_regulons.R
Rscript scripts/step7_clean_figures.R
Rscript scripts/step7b_regulons_and_composition.R
Rscript scripts/step7d_volcano_labels_v2.R
```

Or interactively from an R session, set the working directory to the repo
root and `source()` each script.

All outputs are written to `outliers_included/{figures,tables,rdata}/`
inside `PROJECT_ROOT`. End-to-end run time is approximately 30–45 minutes
on a modern laptop. The slowest steps are GSEA (~5 min) and VIPER (~3 min);
everything else is seconds to a minute or two.

---

## Headline statistics this pipeline reproduces

| Quantity | Value |
|---|---|
| Donors analyzed | 54 (22 COVID-19, 23 Control, 9 ICU/VENT) |
| Age range | 23 – 85 years |
| DEGs (COVID-19 vs Control, padj < 0.05) | **7,079** (4,086 up / 2,993 down) |
| Significantly enriched GO:BP terms | **1,708** (1,169 up / 539 down) |
| Significant TFs (VIPER, padj < 0.05) | **31** (25 activated, 6 suppressed) |
| Top activated TFs | RELA (+6.91), NFKB1 (+6.83), IRF2 (+5.40), STAT3 (+5.29), STAT2 (+5.17), STAT1 (+4.94) |
| Top suppressed TFs | HMBOX1 (-3.89), NR2F1 (-2.93) |
| Aging regulator REST | NES = +3.29, padj = 0.016, 8/11 targets sig |
| Cell types altered (BRETIGEA, padj < 0.05) | 5 of 6 (neurons ↓; astrocytes, microglia, OPCs, endothelial ↑) |
| DEGs cell-intrinsic (composition-corrected) | **110 (1.7%)** |
| DEGs composition-driven | **6,493 (98.3%)** |

---

## Pipeline structure at a glance

```
                ┌──────────┐
                │  Step 1  │  Salmon import + QC + PCA (n=54)
                └────┬─────┘
                     │ dds_n54_all_groups.rds, vsd_n54.rds
                     ▼
                ┌──────────┐
                │  Step 2  │  DESeq2 differential expression (n=45)
                └────┬─────┘
                     │ res_apeglm_n45.rds, res_raw_df_n45.rds
                     ▼
       ┌─────────────┼─────────────┬──────────────┐
       ▼             ▼             ▼              ▼
   ┌───────┐    ┌────────┐    ┌────────┐     ┌────────┐
   │Step 3 │    │ Step 4 │    │ Step 5 │     │ Step 6 │
   │ GSEA  │    │ VIPER  │    │BRETIGEA│     │Heatmap │
   └───┬───┘    └───┬────┘    └────┬───┘     │+regulon│
       │            │              │         └────┬───┘
       └────────────┴──────┬───────┴──────────────┘
                           ▼
                   ┌────────────────┐
                   │  Step 7 / 7b / 7d │  Final clean figure pass
                   └────────────────┘
```

See [`docs/pipeline_overview.md`](docs/pipeline_overview.md) for a
script-by-script walkthrough including inputs, outputs, and
reported statistics.

---

## Methodological notes

### Outlier handling

This is the **outliers-included** sensitivity analysis. The original
dissertation analysis used n=42 COVID-vs-Control samples after removing
three samples (COVID11, COVID14, COVID15) flagged as PCA outliers. The
analysis here re-runs the full pipeline on n=45 (i.e. all 22 COVID-19 and
23 Control samples) to test the robustness of the headline findings to
outlier inclusion.

Each script's console output prints a side-by-side comparison to the
original n=42 numbers (e.g. `step2` shows 5,637 DEGs originally vs 7,079
in this re-run). The findings of the dissertation are robust: every
headline qualitative result (NF-κB activation, REST engagement,
composition-driven dominance) is preserved.

### Composition correction

The composition-corrected model uses BRETIGEA neuronal and astrocyte
scores as continuous covariates added to the DESeq2 design:

```
~ age + sex + neu_score + ast_score + group
```

Including additional cell-type covariates (microglia, OPC, endothelial,
oligodendrocyte) produced model singularity at the n=45 sample size.
Neurons and astrocytes were chosen as the two largest cell-type effects
in the unadjusted Wilcoxon comparison.

### Data citation

This work re-analyzes data originally generated by:

> Mavrikaki M, Lee JD, Solomon IH, Slack FJ. **Severe COVID-19 is
> associated with molecular signatures of aging in the human brain.**
> *Nature Aging* 2:1130–1137 (2022).
> doi: [10.1038/s43587-022-00321-w](https://doi.org/10.1038/s43587-022-00321-w)

The original investigators retain all rights to the underlying tissue
collection and raw data. Original IRB approval was obtained at Beth
Israel Deaconess Medical Center.

### Reproducibility caveats

- Results are deterministic up to floating-point precision when the
  random seed is preserved (`set.seed(42)` is set in the GSEA, VIPER,
  and BRETIGEA steps).
- Package version differences may produce small numerical shifts. This
  pipeline was developed with the package versions listed in the
  `sessionInfo()` block at the bottom of this README.
- DEG counts may differ marginally from the Mavrikaki et al. report
  (7,079 here vs 6,993 originally) due to apeglm shrinkage and DESeq2
  version updates; the qualitative biology is identical.

---

## License

This code is released under the [MIT License](LICENSE).

You are free to use, modify, and redistribute this code with attribution.
The underlying transcriptomic data is governed by the original
investigators' deposition terms at GEO accession GSE188847.

---

## Citation

If you use this code in your own work, please cite both the manuscript
and this repository. See [`CITATION.cff`](CITATION.cff) for the GitHub-
formatted citation, or use:

```
Elechi KW. Cellular composition shifts dominate the COVID-19 brain
transcriptomic response: NF-κB, IFN/JAK-STAT and REST master regulators
converge on a five-cell-type neurovascular signature in postmortem
frontal cortex. Manuscript in preparation, 2026.

Code repository: https://github.com/kelyto30/covid-brain-aging-thesis
```

---

## Contact

**Kelechi Wisdom Elechi**
PhD Candidate, Integrated Biomedical Sciences
Sam and Ann Barshop Institute for Longevity and Aging Studies
University of Texas Health Science Center at San Antonio

Email: kelechiwisdom20@gmail.com
GitHub: [@kelyto30](https://github.com/kelyto30)

---

## sessionInfo placeholder

Run `sessionInfo()` after sourcing all scripts to capture exact package
versions for reproducibility. Save the output as `sessionInfo.txt` and
commit alongside the code if reproducing in a new environment.

```r
# Example
cat(capture.output(sessionInfo()), file = "sessionInfo.txt", sep = "\n")
```
