## ============================================================
##  step7_clean_figures.R   (v2 — reads actual RDS filenames)
##
##  Re-plot all 13 dissertation figures with NEUTRAL subtitles —
##  no "outliers retained", no triangle markers for COVID11/14/15,
##  no OrigQC legend.
##
##  Reads the raw analysis objects produced by Steps 1–6 and
##  derives plot dataframes inline.
##
##  Run on local Mac. Outputs overwrite PNGs in
##  outliers_included/figures/.
## ============================================================

## ---- 0. Paths & libraries ----------------------------------
# Load paths from config.R at repo root
if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists("../config.R")) {
  source("../config.R")
} else {
  stop("Cannot find config.R. Set your working directory to the ",
       "repository root before running this script.")
}
proj   <- PROJECT_ROOT
rd_dir <- OUT_RDS
fg_dir <- OUT_FIG
stopifnot(dir.exists(rd_dir), dir.exists(fg_dir))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(DESeq2)
  library(SummarizedExperiment)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
})

## Common neutral theme used throughout
neutral_theme <- function(base = 13) {
  theme_bw(base_size = base) +
    theme(
      plot.title       = element_text(face = "bold", size = base + 2),
      plot.subtitle    = element_text(size = base - 1, color = "grey30"),
      panel.grid.minor = element_blank()
    )
}

## ---- 1. Load core objects ----------------------------------
cat("Loading RDS objects...\n")
meta_n54        <- readRDS(file.path(rd_dir, "meta_n54.rds"))
vsd_n54         <- readRDS(file.path(rd_dir, "vsd_n54.rds"))
dds_n54         <- readRDS(file.path(rd_dir, "dds_n54_all_groups.rds"))
bret_n54        <- readRDS(file.path(rd_dir, "bret_n54.rds"))

dds_n45         <- readRDS(file.path(rd_dir, "dds_cc_n45.rds"))
dds_corr        <- readRDS(file.path(rd_dir, "dds_corr_n45.rds"))
res_apeglm      <- readRDS(file.path(rd_dir, "res_apeglm_n45.rds"))
res_raw         <- readRDS(file.path(rd_dir, "res_raw_n45.rds"))
res_corr_apeglm <- readRDS(file.path(rd_dir, "res_corr_apeglm_n45.rds"))

gsea_bp         <- readRDS(file.path(rd_dir, "gsea_bp_n45.rds"))
msviper_obj     <- readRDS(file.path(rd_dir, "msviper_n45.rds"))
viper_df        <- readRDS(file.path(rd_dir, "viper_df_n45.rds"))

## ---- normalize meta column names ---------------------------
mcols <- colnames(meta_n54)
cat("meta_n54 columns:", paste(mcols, collapse = ", "), "\n")

sample_col <- intersect(c("sample","sample_id","geo_accession","Sample"), mcols)[1]
group_col  <- intersect(c("group","Group","condition"), mcols)[1]
sex_col    <- intersect(c("sex","Sex"), mcols)[1]
age_col    <- intersect(c("age","Age"), mcols)[1]

stopifnot(!is.na(sample_col), !is.na(group_col), !is.na(sex_col))

meta <- meta_n54
meta$sample <- meta[[sample_col]]
meta$group  <- meta[[group_col]]
meta$sex    <- meta[[sex_col]]
if (!is.na(age_col)) meta$age <- meta[[age_col]]

# Normalize group labels
meta$group <- as.character(meta$group)
meta$group[grepl("control",   meta$group, ignore.case = TRUE)] <- "Control"
meta$group[grepl("covid",     meta$group, ignore.case = TRUE) &
           !grepl("icu|vent", meta$group, ignore.case = TRUE)] <- "COVID-19"
meta$group[grepl("icu|vent",  meta$group, ignore.case = TRUE)] <- "ICU/VENT"
meta$group <- factor(meta$group, levels = c("Control","COVID-19","ICU/VENT"))

# Normalize sex labels
meta$sex <- as.character(meta$sex)
meta$sex[grepl("^m", meta$sex, ignore.case = TRUE)] <- "M"
meta$sex[grepl("^f", meta$sex, ignore.case = TRUE)] <- "F"
meta$sex <- factor(meta$sex, levels = c("M","F"))

GROUP_COLORS <- c(Control = "#2E86C1",
                  `COVID-19` = "#C0392B",
                  `ICU/VENT` = "#27AE60")

## Helpers
res_to_df <- function(r) {
  if (inherits(r, "DESeqResults")) {
    df <- as.data.frame(r)
    df$gene <- rownames(df)
    return(df)
  }
  df <- as.data.frame(r)
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  df
}

## ============================================================
## FIGURE 1 — PCA of all 54 samples
## ============================================================
cat("\n[1/13] Figure 1: PCA n=54 ...\n")

pca_data <- tryCatch(readRDS(file.path(rd_dir, "pca_data_n54.rds")),
                     error = function(e) NULL)

build_pca <- function() {
  vst_mat <- assay(vsd_n54)
  rv  <- matrixStats::rowVars(vst_mat)
  top <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(vst_mat[top, ]))
  pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  df  <- data.frame(
    sample = rownames(pca$x),
    PC1    = pca$x[, 1],
    PC2    = pca$x[, 2]
  )
  df <- merge(df, meta[, c("sample","group")], by = "sample", all.x = TRUE)
  list(df = df, pct = pct[1:2])
}

if (!is.null(pca_data) && all(c("PC1","PC2") %in% colnames(pca_data))) {
  pca_df <- as.data.frame(pca_data)
  if (!"sample" %in% colnames(pca_df)) pca_df$sample <- rownames(pca_df)
  if (!"group"  %in% colnames(pca_df)) {
    pca_df <- merge(pca_df, meta[, c("sample","group")], by = "sample", all.x = TRUE)
  } else {
    pca_df$group <- as.character(pca_df$group)
    pca_df$group[grepl("control",  pca_df$group, ignore.case = TRUE)] <- "Control"
    pca_df$group[grepl("covid",    pca_df$group, ignore.case = TRUE) &
                 !grepl("icu|vent",pca_df$group, ignore.case = TRUE)] <- "COVID-19"
    pca_df$group[grepl("icu|vent", pca_df$group, ignore.case = TRUE)] <- "ICU/VENT"
    pca_df$group <- factor(pca_df$group, levels = c("Control","COVID-19","ICU/VENT"))
  }
  pcv_attr <- attr(pca_data, "percentVar")
  pc_var   <- if (!is.null(pcv_attr) && length(pcv_attr) >= 2) pcv_attr[1:2] else build_pca()$pct
} else {
  built  <- build_pca()
  pca_df <- built$df
  pc_var <- built$pct
}

g1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3.2, alpha = 0.85) +
  scale_color_manual(values = GROUP_COLORS) +
  labs(
    title    = "Diagnostic PCA — frontal cortex (n = 54)",
    subtitle = sprintf("PC1: %.0f%% variance   |   PC2: %.0f%% variance",
                       pc_var[1] * 100, pc_var[2] * 100),
    x = sprintf("PC1 (%.0f%%)", pc_var[1] * 100),
    y = sprintf("PC2 (%.0f%%)", pc_var[2] * 100),
    color = "Group"
  ) +
  neutral_theme() + theme(legend.position = "right")

ggsave(file.path(fg_dir, "Figure1_PCA_n54.png"),
       plot = g1, width = 8, height = 5.5, dpi = 300, bg = "white")

## ============================================================
## FIGURE 2 — Library size bar chart (no per-sample labels)
## ============================================================
cat("[2/13] Figure 2: library sizes n=54 ...\n")

lib_M <- colSums(counts(dds_n54)) / 1e6
lib_df <- data.frame(sample = names(lib_M), lib_M = as.numeric(lib_M))
lib_df <- merge(lib_df, meta[, c("sample","group")], by = "sample", all.x = TRUE)
lib_df <- lib_df %>% arrange(group, lib_M) %>%
  mutate(sample = factor(sample, levels = sample))

g2 <- ggplot(lib_df, aes(x = sample, y = lib_M, fill = group)) +
  geom_col(width = 0.85) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(
    title    = "Library size per sample (n = 54)",
    subtitle = "Gene-level mapped read counts after tximport summarization",
    x = NULL,
    y = "Mapped reads (millions)",
    fill = "Group"
  ) +
  neutral_theme() +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.position = "right")

ggsave(file.path(fg_dir, "Figure2_library_sizes_n54.png"),
       plot = g2, width = 8, height = 8, dpi = 300, bg = "white")

## ============================================================
## FIGURE 3 — Volcano plot of DEGs
## ============================================================
cat("[3/13] Figure 3: volcano n=45 ...\n")

res_df <- res_to_df(res_apeglm)
if (!"gene" %in% colnames(res_df) && "symbol" %in% colnames(res_df)) {
  res_df$gene <- res_df$symbol
}

res_df <- res_df %>%
  mutate(status = case_when(
    !is.na(padj) & padj < 0.05 & log2FoldChange >  0 ~ "Up",
    !is.na(padj) & padj < 0.05 & log2FoldChange <  0 ~ "Down",
    TRUE                                              ~ "NS"
  ))

n_up   <- sum(res_df$status == "Up",   na.rm = TRUE)
n_down <- sum(res_df$status == "Down", na.rm = TRUE)

label_genes <- c("MT2A","MT1E","SOCS3","CXCL10","MMP9","BCL2A1",
                 "GRIN3A","KIF21B","BDNF","LINGO1","CPNE7","BEGAIN")
lab_df <- res_df %>% filter(gene %in% label_genes, !is.na(padj))

g3 <- ggplot(res_df %>% filter(!is.na(padj)),
             aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.55, size = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text_repel(data = lab_df, aes(label = gene),
                  size = 3.4, color = "black", max.overlaps = Inf,
                  segment.color = "grey40", min.segment.length = 0) +
  scale_color_manual(values = c(Up = "#C0392B", Down = "#2E86C1", NS = "grey70"),
                     labels = c(Up   = sprintf("Up (n=%d)",   n_up),
                                Down = sprintf("Down (n=%d)", n_down),
                                NS   = "NS")) +
  labs(
    title    = "Differential expression: COVID-19 vs Control",
    subtitle = sprintf("%d DEGs at padj < 0.05  |  design: ~ age + sex + group  (n = 45)",
                       n_up + n_down),
    x = "Log2 fold change (COVID-19 / Control) — apeglm shrunken",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  neutral_theme() + theme(legend.position = "right")

ggsave(file.path(fg_dir, "Figure3_volcano_n45.png"),
       plot = g3, width = 8, height = 5.5, dpi = 300, bg = "white")

## ============================================================
## FIGURE 4 — GSEA dot plot
## ============================================================
cat("[4/13] Figure 4: GSEA dotplot n=45 ...\n")

gsea_to_df <- function(g) {
  if (inherits(g, c("gseaResult","enrichResult"))) return(as.data.frame(g))
  as.data.frame(g)
}
gsea_df_full <- gsea_to_df(gsea_bp)

# Harmonize column names
if (!"p.adjust" %in% colnames(gsea_df_full)) {
  alt <- intersect(c("padj","qvalue","FDR","adj.P.Val"), colnames(gsea_df_full))[1]
  if (!is.na(alt)) gsea_df_full$p.adjust <- gsea_df_full[[alt]]
}
stopifnot(all(c("Description","NES","p.adjust","setSize") %in% colnames(gsea_df_full)))

gsea_sig <- gsea_df_full %>% filter(p.adjust < 0.05)
top_up   <- gsea_sig %>% filter(NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = 15)
top_down <- gsea_sig %>% filter(NES < 0) %>% arrange(NES)        %>% slice_head(n = 15)

gsea_top <- bind_rows(top_up, top_down) %>%
  arrange(NES) %>%
  mutate(Description = factor(Description, levels = Description),
         neglogp     = -log10(p.adjust))

g4 <- ggplot(gsea_top,
             aes(x = NES, y = Description,
                 size = setSize, color = neglogp)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#9EB3C2", high = "#C0392B", name = "-log10(padj)") +
  scale_size_continuous(range = c(2.5, 7), name = "Set size") +
  labs(
    title    = "Top 30 GO:BP terms — GSEA",
    subtitle = "15 most upregulated and 15 most downregulated terms (n = 45)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  neutral_theme(base = 11) +
  theme(axis.text.y = element_text(size = 9), legend.position = "right")

ggsave(file.path(fg_dir, "Figure4_GSEA_dotplot_n45.png"),
       plot = g4, width = 9, height = 7.5, dpi = 300, bg = "white")

## ============================================================
## FIGURE 5 — VIPER volcano
## ============================================================
cat("[5/13] Figure 5: VIPER volcano n=45 ...\n")

vdf <- as.data.frame(viper_df)

# Harmonize tf / nes / padj columns
if (!"tf" %in% colnames(vdf)) {
  alt <- intersect(c("Regulon","TF","gene","Gene","name"), colnames(vdf))[1]
  if (!is.na(alt)) vdf$tf <- vdf[[alt]] else vdf$tf <- rownames(vdf)
}
if (!"nes" %in% colnames(vdf)) {
  alt <- intersect(c("NES","es","ES","nes_value"), colnames(vdf))[1]
  if (!is.na(alt)) vdf$nes <- vdf[[alt]]
}
if (!"padj" %in% colnames(vdf)) {
  alt <- intersect(c("FDR","qvalue","p.adjust","adj.P.Val","fdr"), colnames(vdf))[1]
  if (!is.na(alt)) vdf$padj <- vdf[[alt]]
}
stopifnot(all(c("tf","nes","padj") %in% colnames(vdf)))

vdf <- vdf %>%
  mutate(status = case_when(
    !is.na(padj) & padj < 0.05 & nes > 0 ~ "Activated",
    !is.na(padj) & padj < 0.05 & nes < 0 ~ "Suppressed",
    TRUE                                  ~ "NS"
  ))

n_act <- sum(vdf$status == "Activated",  na.rm = TRUE)
n_sup <- sum(vdf$status == "Suppressed", na.rm = TRUE)

label_tfs <- c("RELA","NFKB1","STAT1","STAT2","STAT3","IRF2",
               "REST","NR2F1","HMBOX1","TP53","FOXO3","NFE2L2","RELB")
lab_tf <- vdf %>% filter(tf %in% label_tfs)

g5 <- ggplot(vdf %>% filter(!is.na(padj)),
             aes(x = nes, y = -log10(padj), color = status)) +
  geom_point(size = 2.3, alpha = 0.75) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text_repel(data = lab_tf, aes(label = tf),
                  color = "black", size = 3.4,
                  max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = c(Activated = "#C0392B",
                                Suppressed = "#2E86C1",
                                NS = "grey70"),
                     labels = c(Activated  = sprintf("Activated (n=%d)",  n_act),
                                Suppressed = sprintf("Suppressed (n=%d)", n_sup),
                                NS         = "NS")) +
  labs(
    title    = "VIPER master regulator analysis",
    subtitle = sprintf("DoRothEA A/B/C regulons  |  %d significant TFs at padj < 0.05  (n = 45)",
                       n_act + n_sup),
    x = "Normalized Enrichment Score (NES)",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  neutral_theme() + theme(legend.position = "right")

ggsave(file.path(fg_dir, "Figure5_VIPER_volcano_n45.png"),
       plot = g5, width = 8, height = 5.5, dpi = 300, bg = "white")

## ============================================================
## FIGURE 6 — Aging-relevant TF barplot
## ============================================================
cat("[6/13] Figure 6: aging TF barplot n=45 ...\n")

aging_panel <- tribble(
  ~tf,       ~pathway,
  "RELA",    "NF-κB / TNF",
  "NFKB1",   "NF-κB / TNF",
  "RELB",    "NF-κB / TNF",
  "STAT1",   "IFN/JAK-STAT",
  "STAT2",   "IFN/JAK-STAT",
  "STAT3",   "IFN/JAK-STAT",
  "IRF2",    "IFN/JAK-STAT",
  "REST",    "Neuronal aging",
  "NR2F1",   "Synaptic regulator",
  "HMBOX1",  "Synaptic regulator",
  "TP53",    "DNA damage / stress",
  "FOXO3",   "DNA damage / stress",
  "NFE2L2",  "DNA damage / stress",
  "MYC",     "Other",
  "E2F1",    "Other",
  "SP1",     "Other",
  "CEBPA",   "Other",
  "SPI1",    "Other"
)

aging_tf <- aging_panel %>%
  inner_join(vdf %>% select(tf, nes, padj), by = "tf") %>%
  mutate(
    sig_label = case_when(
      is.na(padj)  ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      padj < 0.10  ~ "†",
      TRUE         ~ ""
    ),
    alpha_val = ifelse(!is.na(padj) & padj < 0.10, 1, 0.35)
  ) %>%
  arrange(pathway, nes) %>%
  mutate(tf = factor(tf, levels = tf))

g6 <- ggplot(aging_tf, aes(x = nes, y = tf, fill = pathway, alpha = alpha_val)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sig_label,
                x = nes + ifelse(nes > 0, 0.25, -0.25)),
            size = 3.5, color = "black", show.legend = FALSE) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_alpha_identity() +
  scale_fill_manual(values = c(
    "IFN/JAK-STAT"        = "#1C7293",
    "NF-κB / TNF"         = "#C0392B",
    "Synaptic regulator"  = "#27AE60",
    "Neuronal aging"      = "#8E44AD",
    "DNA damage / stress" = "#B8860B",
    "Other"               = "grey50"
  )) +
  labs(
    title    = "Aging-relevant transcription factor activity",
    subtitle = "Curated panel of 18 TFs grouped by pathway (n = 45)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    fill = "Pathway"
  ) +
  neutral_theme() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 11))

ggsave(file.path(fg_dir, "Figure6_aging_TF_barplot_n45.png"),
       plot = g6, width = 9, height = 6.5, dpi = 300, bg = "white")

## ============================================================
## FIGURE 7 — BRETIGEA cell-type boxplots
## ============================================================
cat("[7/13] Figure 7: cell-type boxplots n=54 ...\n")

bret_mat <- as.data.frame(bret_n54)
if (!"sample" %in% colnames(bret_mat)) bret_mat$sample <- rownames(bret_mat)

ct_cols <- setdiff(colnames(bret_mat), c("sample","group","sex","age"))

ct_rename <- c(
  neu = "Neurons", ast = "Astrocytes", mic = "Microglia",
  oli = "Oligodendrocytes", opc = "OPCs", end = "Endothelial",
  Neurons = "Neurons", Astrocytes = "Astrocytes", Microglia = "Microglia",
  Oligodendrocytes = "Oligodendrocytes", OPCs = "OPCs", Endothelial = "Endothelial"
)

ct_long <- bret_mat %>%
  select(sample, all_of(ct_cols)) %>%
  pivot_longer(-sample, names_to = "cell_type_raw", values_to = "score") %>%
  mutate(cell_type = ifelse(cell_type_raw %in% names(ct_rename),
                            ct_rename[cell_type_raw], cell_type_raw)) %>%
  inner_join(meta[, c("sample","group")], by = "sample")

ct_long$cell_type <- factor(
  ct_long$cell_type,
  levels = intersect(c("Neurons","Astrocytes","Microglia","Oligodendrocytes","OPCs","Endothelial"),
                     unique(ct_long$cell_type))
)

g7 <- ggplot(ct_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.7) +
  geom_jitter(width = 0.18, size = 1.4, alpha = 0.65, color = "black") +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(
    title    = "BRETIGEA cell-type enrichment scores by group (n = 54)",
    subtitle = "Six brain cell types  |  Wilcoxon rank-sum tests reported in summary table",
    x = NULL,
    y = "Cell-type enrichment score"
  ) +
  neutral_theme() +
  theme(legend.position = "none",
        strip.text     = element_text(face = "bold"))

ggsave(file.path(fg_dir, "Figure7_celltype_boxplots_n54.png"),
       plot = g7, width = 9, height = 6, dpi = 300, bg = "white")

## ============================================================
## FIGURE 8 — Sample-sample correlation heatmap (NO OrigQC)
## ============================================================
cat("[8/13] Figure 8: sample correlation heatmap n=54 ...\n")

vst_mat    <- assay(vsd_n54)
sample_cor <- cor(vst_mat, method = "pearson")

ann <- data.frame(
  Group = meta$group[match(colnames(sample_cor), meta$sample)],
  Sex   = meta$sex[  match(colnames(sample_cor), meta$sample)],
  row.names = colnames(sample_cor)
)

ann_colors <- list(
  Group = c(Control = "#2E86C1", `COVID-19` = "#C0392B", `ICU/VENT` = "#27AE60"),
  Sex   = c(M = "#444444", F = "#BBBBBB")
)

png(file.path(fg_dir, "Figure8_sample_correlation_heatmap_n54.png"),
    width = 11, height = 10, units = "in", res = 300)
pheatmap(
  sample_cor,
  annotation_col    = ann,
  annotation_row    = ann,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 6,
  fontsize_col      = 6,
  color             = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
  main              = "Sample-sample Pearson correlation (VST, n = 54)",
  border_color      = "grey90",
  treeheight_row    = 30,
  treeheight_col    = 30
)
dev.off()

## ============================================================
## FIGURES 9a / 9b / 9c — VIPER regulon plots
## ============================================================
cat("[9/13]  Figure 9a: RELA regulon ...\n")
cat("[10/13] Figure 9b: NR2F1 regulon ...\n")
cat("[11/13] Figure 9c: REST regulon ...\n")

# msviper_obj contains $regulon — use it
regulon_list <- if (!is.null(msviper_obj$regulon)) {
  msviper_obj$regulon
} else if (!is.null(attr(msviper_obj, "regulon"))) {
  attr(msviper_obj, "regulon")
} else {
  NULL
}

stopifnot(!is.null(regulon_list))

res_apeglm_df <- res_to_df(res_apeglm)
if (!"gene" %in% colnames(res_apeglm_df)) res_apeglm_df$gene <- rownames(res_apeglm)

build_regulon_df <- function(tf_name) {
  if (is.null(regulon_list[[tf_name]])) {
    stop(sprintf("Regulon for %s not found in msviper$regulon", tf_name))
  }
  tfmode <- regulon_list[[tf_name]]$tfmode
  df <- data.frame(target = names(tfmode), mor = as.numeric(tfmode))
  df <- df %>%
    inner_join(res_apeglm_df %>% select(gene, log2FoldChange, padj),
               by = c("target" = "gene")) %>%
    rename(lfc = log2FoldChange) %>%
    filter(!is.na(padj))
  df
}

get_tf_stats <- function(tf_name) {
  row <- vdf %>% filter(tf == tf_name)
  if (nrow(row) == 0) {
    list(nes = NA_real_, padj = NA_real_)
  } else {
    list(nes = as.numeric(row$nes[1]), padj = as.numeric(row$padj[1]))
  }
}

plot_regulon <- function(tf_name, out_file) {
  df    <- build_regulon_df(tf_name)
  stats <- get_tf_stats(tf_name)

  df <- df %>%
    mutate(
      direction  = ifelse(mor > 0, "Activated by TF", "Repressed by TF"),
      status_lab = case_when(
        padj < 0.05 & lfc >  0 ~ "Up in COVID",
        padj < 0.05 & lfc <  0 ~ "Down in COVID",
        TRUE                   ~ "NS"
      )
    )

  n_total <- nrow(df)
  n_sig   <- sum(df$status_lab != "NS")

  top_lab <- df %>%
    filter(status_lab != "NS") %>%
    group_by(direction) %>%
    arrange(desc(abs(lfc))) %>%
    slice_head(n = 6) %>%
    ungroup()

  g <- ggplot(df, aes(x = direction, y = lfc, color = status_lab)) +
    geom_boxplot(aes(group = direction), outlier.shape = NA,
                 fill = "white", color = "grey30", width = 0.55) +
    geom_jitter(aes(size = ifelse(status_lab != "NS", 2.5, 1.0)),
                width = 0.18, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_lab, aes(label = target),
                    size = 3.3, max.overlaps = Inf,
                    min.segment.length = 0, segment.color = "grey50",
                    show.legend = FALSE) +
    scale_color_manual(values = c(`Up in COVID`   = "#C0392B",
                                  `Down in COVID` = "#2E86C1",
                                  NS              = "grey70"),
                       breaks = c("Up in COVID","Down in COVID","NS")) +
    scale_size_identity() +
    labs(
      title    = sprintf("%s regulon activity in COVID-19 brains", tf_name),
      subtitle = sprintf("NES = %+.2f  |  padj = %s  |  %d targets in dataset  |  %d significant DEGs  (n = 45)",
                         stats$nes,
                         format(stats$padj, digits = 2, scientific = TRUE),
                         n_total, n_sig),
      x = sprintf("Regulatory relationship to %s", tf_name),
      y = "Log2 Fold Change (COVID vs Control)",
      color = "DEG status"
    ) +
    neutral_theme() + theme(legend.position = "right")

  ggsave(out_file, plot = g, width = 8.5, height = 5.5, dpi = 300, bg = "white")
}

plot_regulon("RELA",  file.path(fg_dir, "Figure9a_VIPER_RELA_regulon_n45.png"))
plot_regulon("NR2F1", file.path(fg_dir, "Figure9b_VIPER_NR2F1_regulon_n45.png"))
plot_regulon("REST",  file.path(fg_dir, "Figure9c_VIPER_REST_regulon_n45.png"))

## ============================================================
## FIGURE 11 — LFC scatter (uncorrected vs corrected)
## ============================================================
cat("[12/13] Figure 11: LFC scatter uncorr vs corr n=45 ...\n")

uncorr_df <- res_to_df(res_apeglm) %>%
  select(gene, lfc_uncorr = log2FoldChange, padj_uncorr = padj)
corr_df <- res_to_df(res_corr_apeglm) %>%
  select(gene, lfc_corr = log2FoldChange, padj_corr = padj)

lfc_df <- inner_join(uncorr_df, corr_df, by = "gene") %>%
  filter(!is.na(padj_uncorr) & padj_uncorr < 0.05) %>%
  mutate(
    status = case_when(
      !is.na(padj_corr) & padj_corr < 0.05 ~ "cell_intrinsic",
      TRUE                                  ~ "composition_driven"
    )
  )

label_lfc <- c("MT2A","MT1E","BCL2A1","CXCL10","STAT1",
               "GRIN3A","KIF21B","BDNF","SOCS3","MMP9")
lab_lfc_df <- lfc_df %>% filter(gene %in% label_lfc)

g11 <- ggplot(lfc_df, aes(x = lfc_uncorr, y = lfc_corr, color = status)) +
  geom_point(alpha = 0.5, size = 1.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_text_repel(data = lab_lfc_df, aes(label = gene),
                  size = 3.3, color = "black",
                  max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = c(cell_intrinsic     = "#C0392B",
                                composition_driven = "#9EB3C2"),
                     labels = c(cell_intrinsic     = "Cell-intrinsic",
                                composition_driven = "Composition-driven")) +
  labs(
    title    = "DEG effect sizes: uncorrected vs composition-corrected models",
    subtitle = "Diagonal = preserved after correction (cell-intrinsic);  collapsed = composition-driven  (n = 45)",
    x = "Log2 fold change — uncorrected model",
    y = "Log2 fold change — composition-corrected model",
    color = NULL
  ) +
  neutral_theme() + theme(legend.position = "right")

ggsave(file.path(fg_dir, "Figure11_LFC_uncorr_vs_corr_n45.png"),
       plot = g11, width = 8, height = 5.5, dpi = 300, bg = "white")

## ============================================================
## FIGURE 12 — Composition vs cell-intrinsic summary
## ============================================================
cat("[13/13] Figure 12: composition summary n=45 ...\n")

n_total_uncorr <- nrow(lfc_df)
n_intrinsic    <- sum(lfc_df$status == "cell_intrinsic")
n_compdriven   <- n_total_uncorr - n_intrinsic

comp_summary <- data.frame(
  category = factor(c("Cell-intrinsic","Composition-driven"),
                    levels = c("Composition-driven","Cell-intrinsic")),
  n        = c(n_intrinsic, n_compdriven)
) %>%
  mutate(
    pct   = 100 * n / sum(n),
    label = sprintf("%s\n%d genes  (%.1f%%)", category, n, pct)
  )

g12 <- ggplot(comp_summary, aes(x = "", y = pct, fill = category)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c(`Cell-intrinsic`     = "#C0392B",
                               `Composition-driven` = "#1C7293")) +
  labs(
    title    = "Composition-driven vs cell-intrinsic DEGs",
    subtitle = sprintf("%d total DEGs at padj < 0.05  (n = 45)", n_total_uncorr),
    x = NULL,
    y = "Percent of DEGs"
  ) +
  neutral_theme() +
  theme(legend.position = "right",
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file.path(fg_dir, "Figure12_composition_summary_n45.png"),
       plot = g12, width = 7, height = 6, dpi = 300, bg = "white")

cat("\n================================================\n")
cat("STEP 7 COMPLETE — clean figures regenerated.\n")
cat("Files overwritten in:\n  ", fg_dir, "\n", sep = "")
cat("================================================\n")
