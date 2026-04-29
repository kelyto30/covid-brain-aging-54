# Pipeline overview

This document walks through the nine analysis scripts in `/scripts/`,
in the order they should be executed. Each script reads inputs from
the previous step's saved RDS objects (in `outliers_included/rdata/`)
and writes new objects, tables, and figures to that same parent folder.

A complete run from raw quants to all 13 manuscript figures takes
approximately 30–45 minutes on a modern laptop. The slowest steps are
GSEA (step 3, ~5 minutes) and VIPER (step 4, ~3 minutes); everything
else is seconds to a minute or two.

---

## Step 1 — `step1_import_qc_pca_n54_v2.R`

Imports Salmon transcript-level quantifications, summarizes to
gene-level counts via `tximport`, runs library-size QC, builds the
54-sample DESeq2 object, and produces a diagnostic PCA.

**Inputs**: 54 Salmon `.gz` quants from GEO; metadata CSV; tx2gene CSV.
**Outputs**:
- `dds_n54_all_groups.rds` (DESeq2 object, all 54 samples)
- `vsd_n54.rds` (variance-stabilized expression matrix)
- `meta_n54.rds` (54-sample metadata)
- `pca_data_n54.rds` (PCA coordinates + 3-SD outlier flag)
- `Figure1_PCA_n54.png`, `Figure2_library_sizes_n54.png`
- `library_sizes_n54.csv`, `pca_coordinates_n54.csv`

Reports the 17%/7% PC1/PC2 variance and confirms n=22/23/9 group sizes.
Originally-flagged outliers (COVID11/14/15) are shown as triangles in
the PCA plot but retained for analysis (this is the
"outliers-included" sensitivity analysis).

---

## Step 2 — `step2_deseq2_n45.R`

Subsets to the 45-sample COVID-19 vs Control comparison (drops
ICU/VENT) and runs DESeq2 differential expression with the design
`~ age + sex + group`. Extracts both the unshrunken Wald-statistic
result (for downstream GSEA ranking) and the apeglm-shrunken result
(for effect-size visualization).

**Inputs**: `dds_n54_all_groups.rds`, `meta_n54.rds`, tx2gene CSV.
**Outputs**:
- `dds_cc_n45.rds`, `res_apeglm_n45.rds`, `res_raw_n45.rds`,
  `res_raw_df_n45.rds`
- `DEG_results_n45_apeglm.csv`, `DEG_results_n45_raw_Wald.csv`
- `Figure3_volcano_n45.png` (an early version; final labeled
  volcano is regenerated in step 7d)

Reports the 7,079 DEG headline (4,086 up / 2,993 down) and a side-
by-side comparison to the original n=42 analysis (5,637 / 2,869 / 2,768).

---

## Step 3 — `step3_gsea_n45.R`

GSEA against GO Biological Process using `clusterProfiler::gseGO`.
Ranking is the unshrunken DESeq2 Wald statistic, deduplicated by
gene symbol (keeping the row with largest |stat| per symbol).

**Inputs**: `res_raw_df_n45.rds`.
**Outputs**:
- `gsea_bp_n45.rds`
- `GSEA_GOBP_n45_outliers_included.csv` (all significant terms)
- `GSEA_aging_relevant_terms_n45.csv` (immune/synaptic/neuronal subset)
- `GSEA_dissertation_term_spotcheck_n45.csv`
- `Figure4_GSEA_dotplot_n45.png` (top 15 up + top 15 down terms)

Reports the 1,708 significant GO:BP terms headline (1,169 up / 539 down).

---

## Step 4 — `step4_viper_n45.R`

VIPER master regulator analysis using DoRothEA confidence A/B/C
regulons. Uses the gene-symbol-keyed Wald-statistic vector as the
gene expression signature. Builds the regulon list in VIPER's
expected format (`tfmode` and `likelihood` slots per TF) and runs
`viper::msviper`.

**Inputs**: `res_raw_df_n45.rds`, DoRothEA package data.
**Outputs**:
- `msviper_n45.rds`, `viper_df_n45.rds`
- `VIPER_results_n45_outliers_included.csv`
- `VIPER_headline_TFs_comparison_n45.csv`
- `Figure5_VIPER_volcano_n45.png` (TF-level volcano)
- `Figure6_aging_TF_barplot_n45.png` (curated aging-relevant TFs)

Reports 252 TFs, 12,088 TF-target interactions, and 31 significant
TFs at padj < 0.05 (25 activated, 6 suppressed). Headline TFs:
NFKB1 (+6.83), RELA (+6.91), STAT1/2/3 and IRF2 (NES +4.94 to +5.40),
REST (+3.29, padj 0.016), HMBOX1 (-3.89), NR2F1 (-2.93).

---

## Step 5 — `step5_bretigea_composition_n45.R`

Three components:

**5A.** BRETIGEA cell-type deconvolution (`brainCells()` with 50
markers per cell type) on the variance-stabilized expression matrix
of all 54 samples. Wilcoxon rank-sum tests compare COVID-19 vs
Control scores per cell type, with BH correction across the 6 tests.

**5B.** Composition-corrected DESeq2 with the design
`~ age + sex + neu_score + ast_score + group` on the 45-sample
subset. Apeglm shrinkage applied to the disease coefficient.

**5C.** Comparison plots: scatter of uncorrected vs corrected log2
fold change, and a summary bar of cell-intrinsic vs composition-driven
gene counts.

**Inputs**: `vsd_n54.rds`, `dds_cc_n45.rds`, `res_raw_df_n45.rds`,
tx2gene CSV.
**Outputs**:
- `bret_n54.rds`, `dds_corr_n45.rds`, `res_corr_apeglm_n45.rds`
- `BRETIGEA_scores_n54.csv`, `BRETIGEA_wilcoxon_n45.csv`
- `DEG_corrected_vs_uncorrected_n45.csv`
- `Figure7_celltype_boxplots_n54.png`
- `Figure11_LFC_uncorr_vs_corr_n45.png`
- `Figure12_composition_summary_n45.png`

Reports the headline composition-correction result: of the
symbol-mapped DEGs in the n=45 analysis, only ~110 (~1.7%) survive
composition correction, while ~6,493 (~98.3%) lose significance —
indicating that the bulk transcriptomic signal is overwhelmingly
composition-driven rather than cell-intrinsic.

---

## Step 6 — `step6_heatmap_regulons.R`

Builds three figure deliverables that depend only on previously
saved RDS objects (no new statistical computation):

1. Sample-sample Pearson correlation heatmap on the VST matrix
   (n = 54), with row/column annotation tracks for Group, Sex, and
   Original QC status.
2. Three VIPER regulon plots (RELA, NR2F1, REST) each showing
   target-level log2 fold change stratified by regulatory
   relationship to the TF.

**Inputs**: `vsd_n54.rds`, `res_apeglm_n45.rds`, `viper_df_n45.rds`,
`meta_n54.rds`, tx2gene CSV.
**Outputs**:
- `Figure8_sample_correlation_heatmap_n54.png`
- `Figure9a_VIPER_RELA_regulon_n45.png`
- `Figure9b_VIPER_NR2F1_regulon_n45.png`
- `Figure9c_VIPER_REST_regulon_n45.png`

Regulon target counts: RELA 121 (72 sig, 60%); NR2F1 64 (30 sig, 47%);
REST 11 (8 sig, 73%).

---

## Step 7 — `step7_clean_figures.R`

Re-plots all 13 manuscript figures with a neutral aesthetic suitable
for publication: drops the "outliers retained" subtitles, removes the
triangle markers for COVID11/14/15, and removes the OrigQC legend.
Reads raw analysis objects from steps 1–6 and rebuilds plot dataframes
inline.

**Outputs**: clean overwrites of all 13 PNGs in `outliers_included/figures/`.

---

## Step 7b — `step7b_regulons_and_composition.R`

Patch script that handles the Ensembl-vs-symbol rownames issue in
some of the regulon and composition figures. Maps Ensembl IDs to
gene symbols via `org.Hs.eg.db`, then regenerates Figures 9a/9b/9c,
11, and 12 cleanly.

---

## Step 7d — `step7d_volcano_labels_v2.R`

Final labeled volcano (Figure 4 of the manuscript). Keeps the full
Ensembl-keyed DESeq2 result for the plot cloud and DEG counts (so
totals match step 2's 7,079 / 4,086 / 2,993), but uses a
symbol-to-Ensembl reverse lookup to label the 9 dissertation genes:

- Up: MT2A, SOCS3, CXCL10, MMP9, BCL2A1
- Down: GRIN3A, KIF21B, BDNF, LINGO1

Replaces an earlier version (`step7c`) that had a symbol-mapping
filter applied to the entire dataframe, which dropped DEG counts
to 6,603 / 3,715 / 2,888.

**Output**: `Figure3_volcano_n45.png` (final labeled version).

---

## Run order summary

```
step1  → step2  → step3  → step4  → step5  → step6
                                                  ↓
                                step7  →  step7b  →  step7d
```

All step 7 scripts depend on outputs from step 1 through step 6.
Step 7b depends on the same RDS pool as step 7. Step 7d only needs
`res_apeglm_n45.rds` from step 2.
