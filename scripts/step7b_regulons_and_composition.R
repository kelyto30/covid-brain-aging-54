## ============================================================
##  step7b_regulons_and_composition.R
##  Patch for step7: rownames are Ensembl IDs, not symbols.
##  Maps to gene symbols via org.Hs.eg.db, then produces:
##    - Figure 9a / 9b / 9c (RELA, NR2F1, REST regulons)
##    - Figure 11 (LFC scatter uncorr vs corr)
##    - Figure 12 (composition summary)
## ============================================================

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

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

neutral_theme <- function(base = 13) {
  theme_bw(base_size = base) +
    theme(
      plot.title       = element_text(face = "bold", size = base + 2),
      plot.subtitle    = element_text(size = base - 1, color = "grey30"),
      panel.grid.minor = element_blank()
    )
}

## ---- Load objects ----
cat("Loading RDS objects...\n")
res_apeglm      <- readRDS(file.path(rd_dir, "res_apeglm_n45.rds"))
res_corr_apeglm <- readRDS(file.path(rd_dir, "res_corr_apeglm_n45.rds"))
msviper_obj     <- readRDS(file.path(rd_dir, "msviper_n45.rds"))
viper_df        <- readRDS(file.path(rd_dir, "viper_df_n45.rds"))

## ---- Helper: convert DESeqResults to data.frame keyed by SYMBOL ----
ens_to_symbol <- function(ens_ids) {
  # Strip Ensembl version suffix if present (ENSG00000000003.14 -> ENSG00000000003)
  ens_clean <- sub("\\..*$", "", ens_ids)
  syms <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = ens_clean,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  unname(syms)
}

res_to_symbol_df <- function(r) {
  df <- as.data.frame(r)
  df$ensembl <- rownames(df)
  df$gene    <- ens_to_symbol(df$ensembl)
  df <- df %>% filter(!is.na(gene), !duplicated(gene))
  df
}

cat("Mapping Ensembl IDs to gene symbols...\n")
res_apeglm_df      <- res_to_symbol_df(res_apeglm)
res_corr_apeglm_df <- res_to_symbol_df(res_corr_apeglm)

cat(sprintf("  res_apeglm:      %d rows mapped to symbols\n",      nrow(res_apeglm_df)))
cat(sprintf("  res_corr_apeglm: %d rows mapped to symbols\n", nrow(res_corr_apeglm_df)))

# Detach AnnotationDbi/org.Hs.eg.db so dplyr verbs (select, rename, filter)
# are no longer masked for the rest of the script
suppressWarnings({
  try(detach("package:org.Hs.eg.db",  unload = TRUE, character.only = TRUE), silent = TRUE)
  try(detach("package:AnnotationDbi", unload = TRUE, character.only = TRUE), silent = TRUE)
})
# Re-attach dplyr to ensure its verbs are at top of search path
library(dplyr)

## ---- Harmonize viper_df ----
vdf <- as.data.frame(viper_df)
if (!"tf"   %in% colnames(vdf)) vdf$tf   <- if ("Regulon" %in% colnames(vdf)) vdf$Regulon else rownames(vdf)
if (!"nes"  %in% colnames(vdf)) vdf$nes  <- if ("NES"     %in% colnames(vdf)) vdf$NES     else NA
if (!"padj" %in% colnames(vdf)) {
  alt <- intersect(c("FDR","qvalue","p.adjust","adj.P.Val"), colnames(vdf))[1]
  if (!is.na(alt)) vdf$padj <- vdf[[alt]]
}

## ============================================================
## FIGURES 9a / 9b / 9c — VIPER regulon plots
## ============================================================

regulon_list <- if (!is.null(msviper_obj$regulon)) {
  msviper_obj$regulon
} else if (!is.null(attr(msviper_obj, "regulon"))) {
  attr(msviper_obj, "regulon")
} else {
  NULL
}
stopifnot(!is.null(regulon_list))

# Detect whether DoRothEA regulons are keyed by symbol or Ensembl
regulon_keys_sample <- head(names(regulon_list), 5)
cat("Regulon list sample keys:", paste(regulon_keys_sample, collapse = ", "), "\n")

# Detect target ID type from one regulon
target_sample <- head(names(regulon_list[[regulon_keys_sample[1]]]$tfmode), 5)
cat("Target sample (first regulon):", paste(target_sample, collapse = ", "), "\n")

target_is_ensembl <- any(grepl("^ENSG", target_sample))
cat(sprintf("Targets are Ensembl IDs: %s\n", target_is_ensembl))

build_regulon_df <- function(tf_name) {
  if (is.null(regulon_list[[tf_name]])) {
    stop(sprintf("Regulon for %s not found", tf_name))
  }
  tfmode <- regulon_list[[tf_name]]$tfmode
  df <- data.frame(
    target_raw = names(tfmode),
    mor        = as.numeric(tfmode),
    stringsAsFactors = FALSE
  )

  if (target_is_ensembl) {
    df$target <- ens_to_symbol(df$target_raw)
    df <- df %>% filter(!is.na(target))
    df <- df %>%
      inner_join(res_apeglm_df %>% dplyr::select(gene, log2FoldChange, padj),
                 by = c("target" = "gene"))
  } else {
    df$target <- df$target_raw
    df <- df %>%
      inner_join(res_apeglm_df %>% dplyr::select(gene, log2FoldChange, padj),
                 by = c("target" = "gene"))
  }

  df %>%
    dplyr::rename(lfc = log2FoldChange) %>%
    dplyr::filter(!is.na(padj))
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

  cat(sprintf("  %s: %d targets joined, %d significant\n",
              tf_name, nrow(df), sum(df$padj < 0.05, na.rm = TRUE)))

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

cat("\n[1/5] Figure 9a: RELA regulon ...\n")
plot_regulon("RELA",  file.path(fg_dir, "Figure9a_VIPER_RELA_regulon_n45.png"))

cat("[2/5] Figure 9b: NR2F1 regulon ...\n")
plot_regulon("NR2F1", file.path(fg_dir, "Figure9b_VIPER_NR2F1_regulon_n45.png"))

cat("[3/5] Figure 9c: REST regulon ...\n")
plot_regulon("REST",  file.path(fg_dir, "Figure9c_VIPER_REST_regulon_n45.png"))

## ============================================================
## FIGURE 11 — LFC scatter (uncorrected vs corrected)
## ============================================================
cat("\n[4/5] Figure 11: LFC scatter uncorr vs corr ...\n")

uncorr_df <- res_apeglm_df %>%
  dplyr::select(gene, lfc_uncorr = log2FoldChange, padj_uncorr = padj)
corr_df <- res_corr_apeglm_df %>%
  dplyr::select(gene, lfc_corr = log2FoldChange, padj_corr = padj)

lfc_df <- inner_join(uncorr_df, corr_df, by = "gene") %>%
  filter(!is.na(padj_uncorr) & padj_uncorr < 0.05) %>%
  mutate(
    status = case_when(
      !is.na(padj_corr) & padj_corr < 0.05 ~ "cell_intrinsic",
      TRUE                                  ~ "composition_driven"
    )
  )

cat(sprintf("  joined %d genes, %d significant uncorr, %d cell-intrinsic\n",
            nrow(inner_join(uncorr_df, corr_df, by = "gene")),
            nrow(lfc_df),
            sum(lfc_df$status == "cell_intrinsic")))

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
cat("[5/5] Figure 12: composition summary ...\n")

n_total      <- nrow(lfc_df)
n_intrinsic  <- sum(lfc_df$status == "cell_intrinsic")
n_compdriven <- n_total - n_intrinsic

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
    subtitle = sprintf("%d total DEGs at padj < 0.05  (n = 45)", n_total),
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
cat("STEP 7b COMPLETE — regulons + composition figures done.\n")
cat("Files in: ", fg_dir, "\n", sep = "")
cat("================================================\n")
