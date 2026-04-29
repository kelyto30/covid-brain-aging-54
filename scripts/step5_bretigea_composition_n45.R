# ================================================================
# Step 5 - Aim 2: BRETIGEA deconvolution + composition-corrected DE
# OUTLIERS-INCLUDED sensitivity analysis
# Dissertation: COVID-19 Brain Aging (GSE188847)
# Kelechi Wisdom Elechi | UT Health San Antonio
# ----------------------------------------------------------------
# Three components:
#   5A  BRETIGEA brainCells() deconvolution on all 54 VST-transformed
#       samples (6 cell types), then Wilcoxon per cell type between
#       COVID-19 and Control (n = 45). Matches dissertation Table 5.
#   5B  Composition-corrected DESeq2:
#       design = ~ age + sex + neu_score + ast_score + group
#       on the 45-sample COVID-vs-Control subset.
#       Reports DEG retention (cell-intrinsic vs composition-driven).
#   5C  Comparison figures: scatter (uncorrected vs corrected LFC)
#       and summary bar (13.3 %/86.7 % equivalent).
# ================================================================


# ---- 0. PATHS (loaded from config.R at repo root) -------------
if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists("../config.R")) {
  source("../config.R")
} else {
  stop("Cannot find config.R. Set your working directory to the ",
       "repository root before running this script.")
}


# ---- 1. LIBRARIES ----------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(BRETIGEA)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(dplyr)
})


# ---- 2. LOAD STEP 1 / STEP 2 OBJECTS --------------------------
vsd_n54    <- readRDS(file.path(OUT_RDS, "vsd_n54.rds"))
dds_cc     <- readRDS(file.path(OUT_RDS, "dds_cc_n45.rds"))
res_raw_df <- readRDS(file.path(OUT_RDS, "res_raw_df_n45.rds"))
tx2gene    <- read.csv(TX2GENE_PATH, stringsAsFactors = FALSE)

cat("VST (n=54):", nrow(vsd_n54), "genes x", ncol(vsd_n54), "samples\n")
cat("dds_cc (n=45):", nrow(dds_cc), "genes x", ncol(dds_cc), "samples\n")


# ---- 3. CONVERT VST MATRIX TO GENE SYMBOLS (for BRETIGEA) -----
vst_mat <- assay(vsd_n54)
rownames(vst_mat) <- sub("\\.[0-9]+$", "", rownames(vst_mat))

tx2gene_u <- unique(tx2gene[, c("GENEID", "GENENAME")])
tx2gene_u$GENEID_clean <- sub("\\.[0-9]+$", "", tx2gene_u$GENEID)

symbol_map <- tx2gene_u$GENENAME[match(rownames(vst_mat), tx2gene_u$GENEID_clean)]

# Drop rows without a symbol
keep_idx <- !is.na(symbol_map) & nzchar(symbol_map)
vst_mat_sym <- vst_mat[keep_idx, ]
rownames(vst_mat_sym) <- symbol_map[keep_idx]

# Collapse duplicate symbols by keeping the row with max mean expression
if (any(duplicated(rownames(vst_mat_sym)))) {
  mean_expr <- rowMeans(vst_mat_sym)
  ord       <- order(-mean_expr)
  vst_mat_sym <- vst_mat_sym[ord, ]
  vst_mat_sym <- vst_mat_sym[!duplicated(rownames(vst_mat_sym)), ]
}
cat("Symbol-keyed VST matrix:", nrow(vst_mat_sym), "symbols x",
    ncol(vst_mat_sym), "samples\n")


# ---- 4. RUN BRETIGEA ------------------------------------------
cat("\nRunning BRETIGEA::brainCells() on all 54 samples ...\n")
set.seed(42)
bret <- brainCells(vst_mat_sym, nMarker = 50, species = "human")
# bret is a samples x 6 matrix (ast, end, mic, neu, oli, opc)
cat("BRETIGEA output:", nrow(bret), "samples x", ncol(bret), "cell types\n")
print(head(round(bret, 3)))

# Save
write.csv(data.frame(sample = rownames(bret), bret),
          file.path(OUT_TBL, "BRETIGEA_scores_n54.csv"),
          row.names = FALSE)
saveRDS(bret, file.path(OUT_RDS, "bret_n54.rds"))


# ---- 5. WILCOXON: COVID-19 vs CONTROL (n = 45) ----------------
meta_n54 <- as.data.frame(colData(vsd_n54))
sample_ids <- rownames(bret)
meta_n54 <- meta_n54[match(sample_ids, meta_n54$sample_id), ]

covid_idx <- meta_n54$group == "COVID-19"
ctrl_idx  <- meta_n54$group == "Control"

cell_labels <- c(neu = "Neurons", ast = "Astrocytes", oli = "Oligodendrocytes",
                 mic = "Microglia", end = "Endothelial", opc = "OPCs")

ct_res <- do.call(rbind, lapply(colnames(bret), function(ct) {
  w    <- wilcox.test(bret[covid_idx, ct], bret[ctrl_idx, ct])
  mC   <- mean(bret[covid_idx, ct])
  mK   <- mean(bret[ctrl_idx, ct])
  data.frame(
    cell_code  = ct,
    cell_type  = cell_labels[[ct]],
    mean_ctrl  = mK,
    mean_covid = mC,
    delta      = mC - mK,
    pvalue     = w$p.value
  )
}))
ct_res$padj <- p.adjust(ct_res$pvalue, method = "BH")
ct_res$sig_label <- dplyr::case_when(
  ct_res$padj <  0.001 ~ "***",
  ct_res$padj <  0.01  ~ "**",
  ct_res$padj <  0.05  ~ "*",
  ct_res$padj <  0.10  ~ "\u2020",
  TRUE                 ~ "NS"
)
ct_res <- ct_res[order(-abs(ct_res$delta)), ]

write.csv(ct_res, file.path(OUT_TBL, "BRETIGEA_wilcoxon_n45.csv"),
          row.names = FALSE)

cat("\n--- BRETIGEA cell-type Wilcoxon (COVID-19 vs Control, n = 22 vs 23) ---\n")
print(ct_res, row.names = FALSE)


# ---- 6. FIGURE 7: CELL-TYPE SCORES BY GROUP -------------------
bret_long <- data.frame(sample = rownames(bret), bret, check.names = FALSE) %>%
  pivot_longer(cols = -sample, names_to = "cell_code", values_to = "score")
bret_long$group    <- meta_n54$group[match(bret_long$sample, meta_n54$sample_id)]
bret_long$cell_type <- cell_labels[bret_long$cell_code]
bret_long$cell_type <- factor(bret_long$cell_type,
                              levels = c("Neurons", "Astrocytes", "Oligodendrocytes",
                                         "Microglia", "Endothelial", "OPCs"))

p_ct <- ggplot(bret_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.alpha = 0.4, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Control" = "#3498DB",
                               "COVID-19" = "#E74C3C",
                               "ICU/VENT" = "#27AE60")) +
  labs(x = NULL, y = "BRETIGEA cell-type enrichment score",
       title = "Cell-type composition by group (n = 54, outliers retained)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(OUT_FIG, "Figure7_celltype_boxplots_n54.png"),
       p_ct, width = 9, height = 6, dpi = 300)


# ---- 7. COMPOSITION-CORRECTED DESeq2 --------------------------
# Add neu and ast BRETIGEA scores as continuous covariates
colData(dds_cc)$neu_score <- bret[match(colnames(dds_cc), rownames(bret)), "neu"]
colData(dds_cc)$ast_score <- bret[match(colnames(dds_cc), rownames(bret)), "ast"]

stopifnot(!any(is.na(dds_cc$neu_score)))
stopifnot(!any(is.na(dds_cc$ast_score)))

design(dds_cc) <- ~ age + sex + neu_score + ast_score + group

cat("\nRunning composition-corrected DESeq2 ...\n")
dds_corr <- DESeq(dds_cc)

coef_name <- grep("COVID.*Control", resultsNames(dds_corr), value = TRUE)
cat("Using coefficient:", coef_name, "\n")

res_corr_raw <- results(dds_corr, name = coef_name, alpha = 0.05)
res_corr     <- lfcShrink(dds_corr, coef = coef_name, type = "apeglm")

# Save
saveRDS(dds_corr, file.path(OUT_RDS, "dds_corr_n45.rds"))
saveRDS(res_corr, file.path(OUT_RDS, "res_corr_apeglm_n45.rds"))


# ---- 8. MERGE CORRECTED AND UNCORRECTED RESULTS --------------
res_corr_df <- as.data.frame(res_corr)
res_corr_df$gene_id       <- rownames(res_corr_df)
res_corr_df$gene_id_clean <- sub("\\.[0-9]+$", "", res_corr_df$gene_id)
res_corr_df$gene_symbol   <- tx2gene_u$GENENAME[match(res_corr_df$gene_id_clean,
                                                      tx2gene_u$GENEID_clean)]
names(res_corr_df)[names(res_corr_df) == "log2FoldChange"] <- "LFC_corr"
names(res_corr_df)[names(res_corr_df) == "padj"]           <- "padj_corr"
names(res_corr_df)[names(res_corr_df) == "pvalue"]         <- "pvalue_corr"

# Load the uncorrected (Step 2) results for comparison
res_uncorr_df <- readRDS(file.path(OUT_RDS, "res_raw_df_n45.rds"))
names(res_uncorr_df)[names(res_uncorr_df) == "log2FoldChange"] <- "LFC_uncorr"
names(res_uncorr_df)[names(res_uncorr_df) == "padj"]           <- "padj_uncorr"
names(res_uncorr_df)[names(res_uncorr_df) == "pvalue"]         <- "pvalue_uncorr"
names(res_uncorr_df)[names(res_uncorr_df) == "stat"]           <- "stat_uncorr"

merged <- merge(
  res_uncorr_df[, c("gene_id", "gene_symbol", "LFC_uncorr", "padj_uncorr")],
  res_corr_df[, c("gene_id", "LFC_corr", "padj_corr")],
  by = "gene_id")

# Flag status
merged$sig_uncorr <- !is.na(merged$padj_uncorr) & merged$padj_uncorr < 0.05
merged$sig_corr   <- !is.na(merged$padj_corr)   & merged$padj_corr   < 0.05

merged$composition_status <- with(merged, dplyr::case_when(
  sig_uncorr &  sig_corr ~ "Cell-intrinsic",
  sig_uncorr & !sig_corr ~ "Composition-driven",
  !sig_uncorr &  sig_corr ~ "Emerged after correction",
  TRUE                    ~ "NS in both"
))

write.csv(merged,
          file.path(OUT_TBL, "DEG_corrected_vs_uncorrected_n45.csv"),
          row.names = FALSE)


# ---- 9. RETENTION SUMMARY -------------------------------------
n_uncorr  <- sum(merged$sig_uncorr)
n_corr    <- sum(merged$sig_corr)
n_intrinsic  <- sum(merged$composition_status == "Cell-intrinsic")
n_comp_drv   <- sum(merged$composition_status == "Composition-driven")
n_emerged    <- sum(merged$composition_status == "Emerged after correction")

pct_intrinsic <- 100 * n_intrinsic / n_uncorr
pct_comp_drv  <- 100 * n_comp_drv  / n_uncorr

cat("\n================= STEP 5 SUMMARY =================\n")
cat(sprintf("%-36s %14s %14s\n", "Metric", "Original n=42", "New n=45"))
cat(sprintf("%-36s %14d %14d\n",
            "DEGs uncorrected",            5637, n_uncorr))
cat(sprintf("%-36s %14d %14d\n",
            "DEGs composition-corrected",   751, n_corr))
cat(sprintf("%-36s %14s %14d\n",
            "  Cell-intrinsic (both sig)", "(= 751)", n_intrinsic))
cat(sprintf("%-36s %14s %14d\n",
            "  Composition-driven (lost)",
            sprintf("%d (86.7%%)", 4886), n_comp_drv))
cat(sprintf("%-36s %14s %14.1f%%\n",
            "  Retention rate (cell-intrinsic)",
            "13.3%", pct_intrinsic))
cat(sprintf("%-36s %14s %14.1f%%\n",
            "  Composition-driven fraction",
            "86.7%", pct_comp_drv))
cat(sprintf("%-36s %14s %14d\n",
            "  Emerged after correction",   "(not reported)", n_emerged))
cat("===================================================\n")


# ---- 10. FIGURE 11 - SCATTER UNCORR vs CORR LFC ---------------
scat <- merged[merged$sig_uncorr, ]
scat$color_cat <- with(scat, dplyr::case_when(
  sig_corr &  LFC_uncorr > 0 ~ "Cell-intrinsic (up)",
  sig_corr &  LFC_uncorr < 0 ~ "Cell-intrinsic (down)",
  !sig_corr & LFC_uncorr > 0 ~ "Composition-driven (up)",
  !sig_corr & LFC_uncorr < 0 ~ "Composition-driven (down)"
))

# Label key genes
key_genes <- c("MT2A", "GRIN3A", "STAT1", "KIF21B", "SOCS3", "CXCL10",
               "BDNF", "BEGAIN", "MMP9", "BCL2A1", "LINGO1", "CPNE7",
               "MT1E", "HMBOX1")
scat$label <- ifelse(scat$gene_symbol %in% key_genes & !is.na(scat$gene_symbol),
                     scat$gene_symbol, "")

p_scat <- ggplot(scat, aes(x = LFC_uncorr, y = LFC_corr, color = color_cat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_point(alpha = 0.55, size = 0.8) +
  geom_text_repel(aes(label = label), size = 3, color = "black",
                  max.overlaps = 30, box.padding = 0.35,
                  min.segment.length = 0) +
  scale_color_manual(values = c(
    "Cell-intrinsic (up)"     = "#8B0000",
    "Cell-intrinsic (down)"   = "#00008B",
    "Composition-driven (up)" = "#F5B7B1",
    "Composition-driven (down)" = "#AED6F1"
  ), name = NULL) +
  labs(x = "LFC uncorrected (~ age + sex + group)",
       y = "LFC composition-corrected (+ neu_score + ast_score)",
       title = "DEG effect sizes: uncorrected vs composition-corrected",
       subtitle = sprintf("n = 45, outliers retained. %d cell-intrinsic, %d composition-driven.",
                          n_intrinsic, n_comp_drv)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(OUT_FIG, "Figure11_LFC_uncorr_vs_corr_n45.png"),
       p_scat, width = 8.5, height = 6.5, dpi = 300)


# ---- 11. FIGURE 12 - COMPOSITION SUMMARY BAR ------------------
summary_bar <- data.frame(
  category = c("Cell-intrinsic", "Composition-driven"),
  n        = c(n_intrinsic, n_comp_drv),
  pct      = c(pct_intrinsic, pct_comp_drv)
)

p_bar <- ggplot(summary_bar,
                aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            vjust = -0.3, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Cell-intrinsic"     = "#2ECC71",
                               "Composition-driven" = "#E74C3C")) +
  labs(x = NULL, y = "Number of DEGs",
       title = "COVID-19 brain transcriptomic signal: cell-intrinsic vs composition-driven",
       subtitle = sprintf("n = %d total DEGs (uncorrected) partitioned by composition correction",
                          n_uncorr)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, max(summary_bar$n) * 1.15))

ggsave(file.path(OUT_FIG, "Figure12_composition_summary_n45.png"),
       p_bar, width = 7, height = 5.5, dpi = 300)


# ---- 12. KEY GENE COMPOSITION STATUS (Table 6 equiv) ----------
cat("\n--- Key gene composition status (Table 6 spot-check) ---\n")
key_spot <- merged[merged$gene_symbol %in% key_genes &
                   !is.na(merged$gene_symbol), ]
key_spot <- key_spot[order(key_spot$padj_corr), ]
print(key_spot[, c("gene_symbol", "LFC_uncorr", "LFC_corr",
                   "padj_uncorr", "padj_corr", "composition_status")],
      row.names = FALSE)

write.csv(key_spot,
          file.path(OUT_TBL, "key_gene_composition_status_n45.csv"),
          row.names = FALSE)

cat("\nOutputs written to:", OUT, "\n")
