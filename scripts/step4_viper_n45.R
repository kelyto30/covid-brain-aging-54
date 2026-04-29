# ================================================================
# Step 4 - VIPER Master Regulator Analysis
# OUTLIERS-INCLUDED sensitivity analysis (n = 45)
# Dissertation: COVID-19 Brain Aging (GSE188847)
# Kelechi Wisdom Elechi | UT Health San Antonio
# ----------------------------------------------------------------
# Mirrors the dissertation Aim 1 method:
#   signature : DESeq2 Wald statistic (unshrunken, from Step 2)
#   regulons  : DoRothEA A/B/C-confidence (human)
#   engine    : viper::msviper
#   minsize   : 10 target genes per TF
#   threshold : BH padj < 0.05
#
# Input : res_raw_df_n45.rds from Step 2
# Output: TF activity table (Table 4 equiv), VIPER volcano (Fig 5),
#         aging-relevant TF barplot (Fig 6), RDS for downstream
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
  library(viper)
  library(dorothea)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})


# ---- 2. LOAD WALD-STAT SIGNATURE -------------------------------
res_raw <- readRDS(file.path(OUT_RDS, "res_raw_df_n45.rds"))

# Build gene-symbol-keyed signature (deduplicate by max |stat|)
sig_df <- res_raw[!is.na(res_raw$stat) & !is.na(res_raw$gene_symbol) &
                  nzchar(res_raw$gene_symbol), ]
sig_df <- sig_df[order(-abs(sig_df$stat)), ]
sig_df <- sig_df[!duplicated(sig_df$gene_symbol), ]

signature <- setNames(sig_df$stat, sig_df$gene_symbol)
cat("Input signature:", length(signature), "unique gene symbols\n")


# ---- 3. LOAD DoRothEA REGULONS (A/B/C) -------------------------
data("dorothea_hs", package = "dorothea")

regulons_abc <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

n_tfs <- length(unique(regulons_abc$tf))
n_int <- nrow(regulons_abc)
cat("DoRothEA A/B/C regulons:", n_tfs, "TFs,",
    format(n_int, big.mark = ","), "interactions\n")


# ---- 4. CONVERT TO VIPER REGULON FORMAT ------------------------
# VIPER expects a list: one entry per TF, each with $tfmode and $likelihood.
df2regulon <- function(df) {
  regulon_list <- split(df, df$tf)
  lapply(regulon_list, function(r) {
    tfmode <- setNames(r$mor, r$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
}

regulon_viper <- df2regulon(regulons_abc)
cat("VIPER regulon object built:", length(regulon_viper), "TFs\n")


# ---- 5. RUN msviper --------------------------------------------
set.seed(42)
cat("\nRunning viper::msviper ...\n")
msv <- msviper(ges     = signature,
               regulon = regulon_viper,
               minsize = 10,
               verbose = FALSE)

# Assemble results data frame
viper_df <- data.frame(
  TF     = names(msv$es$nes),
  NES    = as.numeric(msv$es$nes),
  pvalue = as.numeric(msv$es$p.value),
  size   = msv$es$size[names(msv$es$nes)],
  stringsAsFactors = FALSE
)
viper_df$padj <- p.adjust(viper_df$pvalue, method = "BH")
viper_df$direction <- ifelse(viper_df$NES > 0, "Activated", "Suppressed")
viper_df$significant <- viper_df$padj < 0.05
viper_df <- viper_df[order(-abs(viper_df$NES)), ]


# ---- 6. SAVE RESULTS ------------------------------------------
write.csv(viper_df,
          file.path(OUT_TBL, "VIPER_results_n45_outliers_included.csv"),
          row.names = FALSE)
saveRDS(msv,      file.path(OUT_RDS, "msviper_n45.rds"))
saveRDS(viper_df, file.path(OUT_RDS, "viper_df_n45.rds"))


# ---- 7. SUMMARY AND DISSERTATION SPOT-CHECK --------------------
n_sig   <- sum(viper_df$significant)
n_up    <- sum(viper_df$significant & viper_df$NES > 0)
n_down  <- sum(viper_df$significant & viper_df$NES < 0)

ORIG_SIG  <- 47; ORIG_UP <- 33; ORIG_DOWN <- 14

cat("\n================= STEP 4 SUMMARY =================\n")
cat(sprintf("%-28s %14s %12s %+10s\n",
            "Metric", "Original n=42", "New n=45", "Change"))
cat(sprintf("%-28s %14d %12d %+10d\n", "Significant TFs (padj<0.05)",
            ORIG_SIG, n_sig, n_sig - ORIG_SIG))
cat(sprintf("%-28s %14d %12d %+10d\n", "  Activated",
            ORIG_UP,   n_up,   n_up   - ORIG_UP))
cat(sprintf("%-28s %14d %12d %+10d\n", "  Suppressed",
            ORIG_DOWN, n_down, n_down - ORIG_DOWN))
cat("===================================================\n")

cat("\n--- Dissertation headline TFs (Table 4 spot-check) ---\n")
headline_tfs <- c("STAT1", "RELA", "NFKB1", "STAT3", "NR2F1",
                  "HMBOX1", "REST", "TP53", "FOXO3", "NFE2L2")
orig_nes <- c(STAT1 = 7.23, RELA = 5.44, NFKB1 = 4.79, STAT3 = 4.51,
              NR2F1 = -4.09, HMBOX1 = NA, REST = NA,
              TP53 = NA, FOXO3 = NA, NFE2L2 = NA)

spot <- viper_df[match(headline_tfs, viper_df$TF), ]
spot$Original_NES <- orig_nes[spot$TF]
spot <- spot[, c("TF", "Original_NES", "NES", "padj", "size", "direction", "significant")]
names(spot)[names(spot) == "NES"] <- "New_NES"
print(spot, row.names = FALSE)
write.csv(spot, file.path(OUT_TBL, "VIPER_headline_TFs_comparison_n45.csv"),
          row.names = FALSE)


# ---- 8. FIGURE 5 - VIPER VOLCANO ------------------------------
vdf <- viper_df
vdf$sig_cat <- ifelse(!vdf$significant, "NS",
                      ifelse(vdf$NES > 0, "Activated", "Suppressed"))
vdf$neglog10p <- -log10(pmax(vdf$padj, 1e-300))

# Label: dissertation headline TFs + top 10 by |NES|
top_by_nes <- vdf$TF[order(-abs(vdf$NES))][1:10]
to_label   <- unique(c(headline_tfs, top_by_nes))
vdf$label  <- ifelse(vdf$TF %in% to_label, vdf$TF, "")

# REST trend annotation (original dissertation reported padj = 0.057)
vdf$rest_trend <- vdf$TF == "REST" & !vdf$significant & vdf$padj < 0.10

p_vol <- ggplot(vdf, aes(x = NES, y = neglog10p, color = sig_cat)) +
  geom_point(aes(shape = rest_trend), size = 2.3, alpha = 0.75) +
  geom_text_repel(aes(label = label), size = 3.1, color = "black",
                  max.overlaps = 30, box.padding = 0.4,
                  min.segment.length = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0,            linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Activated"  = "#E74C3C",
                                "Suppressed" = "#3498DB",
                                "NS"         = "grey70"),
                     name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     guide = "none") +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "-log10(adj. p-value)",
       title = "VIPER master regulators (n = 45, outliers retained)",
       subtitle = sprintf("%d significant TFs (%d activated, %d suppressed) at padj < 0.05",
                          n_sig, n_up, n_down)) +
  theme_classic(base_size = 11)

ggsave(file.path(OUT_FIG, "Figure5_VIPER_volcano_n45.png"),
       p_vol, width = 8, height = 6.5, dpi = 300)


# ---- 9. FIGURE 6 - AGING-RELEVANT TF BARPLOT ------------------
# 21 TFs matching the dissertation's curated aging-relevant set
aging_set <- list(
  "IFN/JAK-STAT"      = c("STAT1", "STAT2", "STAT3", "IRF1", "IRF7", "IRF9"),
  "NF-kB/TNF"         = c("RELA", "NFKB1", "RELB", "NFKB2"),
  "Synaptic"          = c("NR2F1", "NR2F2", "HMBOX1", "FOXG1"),
  "Neuronal aging"    = c("REST"),
  "DNA damage/stress" = c("TP53", "FOXO3", "NFE2L2", "E2F1", "MYC")
)

tf_cat <- do.call(rbind, lapply(names(aging_set), function(cat) {
  data.frame(TF = aging_set[[cat]], category = cat, stringsAsFactors = FALSE)
}))

bar_df <- merge(tf_cat, viper_df, by = "TF", all.x = TRUE)
bar_df <- bar_df[!is.na(bar_df$NES), ]   # drop TFs that weren't tested
bar_df$sig_label <- dplyr::case_when(
  bar_df$padj <  0.001 ~ "***",
  bar_df$padj <  0.01  ~ "**",
  bar_df$padj <  0.05  ~ "*",
  bar_df$padj <  0.10  ~ "\u2020",   # dagger for trend
  TRUE                 ~ ""
)
bar_df$alpha_level <- ifelse(bar_df$padj < 0.10, 1.0, 0.5)
bar_df <- bar_df[order(bar_df$NES), ]
bar_df$TF <- factor(bar_df$TF, levels = bar_df$TF)

cat_colors <- c("IFN/JAK-STAT"      = "#E74C3C",
                "NF-kB/TNF"         = "#E67E22",
                "Synaptic"          = "#3498DB",
                "Neuronal aging"    = "#8E44AD",
                "DNA damage/stress" = "#8B4513")

p_bar <- ggplot(bar_df, aes(x = NES, y = TF, fill = category, alpha = alpha_level)) +
  geom_col() +
  geom_text(aes(label = sig_label,
                x = NES + ifelse(NES > 0, 0.25, -0.25)),
            size = 5, fontface = "bold") +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = cat_colors, name = "Pathway") +
  scale_alpha_identity() +
  labs(x = "Normalized Enrichment Score (NES)",
       y = NULL,
       title = "Aging-relevant transcription factors (n = 45, outliers retained)",
       subtitle = "* padj<0.05  ** padj<0.01  *** padj<0.001  \u2020 trend padj<0.10") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(OUT_FIG, "Figure6_aging_TF_barplot_n45.png"),
       p_bar, width = 8, height = 7, dpi = 300)


# ---- 10. FINAL CONSOLE TABLE ----------------------------------
cat("\n--- Aging-relevant TFs (Table 4 + barplot) ---\n")
print(bar_df[order(-bar_df$NES),
             c("TF", "category", "NES", "padj", "size", "sig_label")],
      row.names = FALSE)

cat("\nOutputs written to:", OUT, "\n")
