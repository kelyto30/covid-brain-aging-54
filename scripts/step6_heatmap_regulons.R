# ================================================================
# Step 6 - Sample correlation heatmap + VIPER regulon plots
# OUTLIERS-INCLUDED sensitivity analysis
# Dissertation: COVID-19 Brain Aging (GSE188847)
# Kelechi Wisdom Elechi | UT Health San Antonio
# ----------------------------------------------------------------
# Four figure deliverables, no new statistical computation:
#   Figure 8  Sample-sample Pearson correlation heatmap (n = 54)
#             Counterpart to QC_SampleCorrelation_GSE188847.pdf
#   Figure 9a VIPER RELA regulon plot  (n = 45)
#   Figure 9b VIPER NR2F1 regulon plot (n = 45)
#   Figure 9c VIPER REST regulon plot  (n = 45)
#             Counterparts to VIPER_{RELA,NR2F1,REST}_regulon.png
#
# All stats come from previously-saved RDS objects; this script
# only reshapes and plots.
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
  library(pheatmap)
  library(dorothea)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})


# ---- 2. LOAD OBJECTS ------------------------------------------
vsd_n54    <- readRDS(file.path(OUT_RDS, "vsd_n54.rds"))
res_apeglm <- readRDS(file.path(OUT_RDS, "res_apeglm_n45.rds"))
viper_df   <- readRDS(file.path(OUT_RDS, "viper_df_n45.rds"))
meta_n54   <- readRDS(file.path(OUT_RDS, "meta_n54.rds"))
tx2gene    <- read.csv(TX2GENE_PATH, stringsAsFactors = FALSE)

cat("Loaded:\n")
cat("  vsd_n54:   ", nrow(vsd_n54), "genes x", ncol(vsd_n54), "samples\n")
cat("  res_apeglm:", nrow(res_apeglm), "genes (n=45 COVID-19 vs Control)\n")
cat("  viper_df:  ", nrow(viper_df), "TFs tested\n")


# ---- 3. SAMPLE-SAMPLE CORRELATION HEATMAP (Figure 8) ----------
vst_mat <- assay(vsd_n54)
cor_mat <- cor(vst_mat, method = "pearson")

# Annotation: Group, Sex, Original QC status
annot <- data.frame(
  Group  = meta_n54$group,
  Sex    = meta_n54$sex,
  OrigQC = ifelse(meta_n54$sample_id %in% c("COVID11", "COVID14", "COVID15"),
                  "Removed orig.", "Retained orig.")
)
rownames(annot) <- meta_n54$sample_id

annot_colors <- list(
  Group  = c("Control" = "#3498DB", "COVID-19" = "#E74C3C", "ICU/VENT" = "#27AE60"),
  Sex    = c("M" = "#444444", "F" = "#BBBBBB"),
  OrigQC = c("Retained orig." = "#FFFFFF", "Removed orig." = "#FFD700")
)

cat("\nBuilding sample-sample correlation heatmap (n = 54)...\n")
png(file.path(OUT_FIG, "Figure8_sample_correlation_heatmap_n54.png"),
    width = 12, height = 11, units = "in", res = 300)
pheatmap(cor_mat,
         annotation_col    = annot,
         annotation_row    = annot,
         annotation_colors = annot_colors,
         color             = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
         fontsize          = 7,
         fontsize_row      = 6,
         fontsize_col      = 6,
         main              = "Sample-sample Pearson correlation (VST)  |  n = 54 (outliers retained)")
dev.off()
cat("  wrote Figure8_sample_correlation_heatmap_n54.png\n")


# ---- 4. PREP SHARED DE TABLE FOR REGULON PLOTS ----------------
data("dorothea_hs", package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

# Build DE table keyed by gene symbol (shrunken apeglm LFC)
de_df <- as.data.frame(res_apeglm)
de_df$gene_id       <- rownames(de_df)
de_df$gene_id_clean <- sub("\\.[0-9]+$", "", de_df$gene_id)

tx2gene_u <- unique(tx2gene[, c("GENEID", "GENENAME")])
tx2gene_u$GENEID_clean <- sub("\\.[0-9]+$", "", tx2gene_u$GENEID)
de_df$gene_symbol <- tx2gene_u$GENENAME[match(de_df$gene_id_clean, tx2gene_u$GENEID_clean)]

# Keep one row per symbol (the one with the smallest padj)
de_df <- de_df[!is.na(de_df$gene_symbol) & nzchar(de_df$gene_symbol), ]
de_df <- de_df[order(de_df$padj, na.last = NA), ]
de_df <- de_df[!duplicated(de_df$gene_symbol), ]

# DEG status
de_df$deg_status <- dplyr::case_when(
  !is.na(de_df$padj) & de_df$padj < 0.05 & de_df$log2FoldChange > 0 ~ "Up in COVID",
  !is.na(de_df$padj) & de_df$padj < 0.05 & de_df$log2FoldChange < 0 ~ "Down in COVID",
  TRUE ~ "NS"
)


# ---- 5. REGULON-PLOT FUNCTION ---------------------------------
fmt_p <- function(p) {
  if (is.na(p))        return("NA")
  if (p < 0.001)       return(sprintf("%.2e", p))
  if (p < 0.01)        return(sprintf("%.4f", p))
  sprintf("%.3f", p)
}

make_regulon_plot <- function(tf_name, regulons, de_df, viper_df,
                              outpath, n_label = 6) {
  reg <- regulons %>% filter(tf == tf_name)
  reg$relationship <- ifelse(reg$mor > 0, "Activated by TF", "Repressed by TF")

  reg_de <- merge(reg, de_df, by.x = "target", by.y = "gene_symbol")

  tf_stats <- viper_df[viper_df$TF == tf_name, ]
  nes  <- tf_stats$NES
  padj <- tf_stats$padj

  n_targets  <- nrow(reg_de)
  n_sig_degs <- sum(reg_de$deg_status != "NS")

  # Label the top |LFC| significant genes per regulatory relationship
  label_genes <- character(0)
  for (rel in unique(reg_de$relationship)) {
    sub <- reg_de[reg_de$relationship == rel & reg_de$deg_status != "NS", ]
    sub <- sub[order(-abs(sub$log2FoldChange)), ]
    label_genes <- c(label_genes, head(sub$target, n_label))
  }
  reg_de$label <- ifelse(reg_de$target %in% label_genes, reg_de$target, "")

  # Order factor levels so Activated is on the left
  reg_de$relationship <- factor(reg_de$relationship,
                                levels = c("Activated by TF", "Repressed by TF"))

  p <- ggplot(reg_de, aes(x = relationship, y = log2FoldChange)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white") +
    geom_jitter(aes(color = deg_status, size = deg_status),
                width = 0.18, alpha = 0.8) +
    geom_text_repel(aes(label = label, color = deg_status),
                    size = 3, max.overlaps = 25,
                    box.padding = 0.4, min.segment.length = 0,
                    fontface = "italic", show.legend = FALSE) +
    scale_color_manual(values = c("Up in COVID"   = "#E74C3C",
                                  "Down in COVID" = "#3498DB",
                                  "NS"            = "grey70"),
                       breaks = c("Up in COVID", "Down in COVID", "NS"),
                       name   = "DEG status") +
    scale_size_manual(values = c("Up in COVID"   = 2.8,
                                 "Down in COVID" = 2.8,
                                 "NS"            = 1.1),
                      guide  = "none") +
    labs(x = paste("Regulatory relationship to", tf_name),
         y = "Log2 Fold Change (COVID vs Control)",
         title = sprintf("%s Regulon Activity in COVID-19 Brains  |  n = 45 (outliers retained)",
                         tf_name),
         subtitle = sprintf("NES = %.2f  |  padj = %s  |  %d targets in dataset  |  %d significant DEGs",
                            nes, fmt_p(padj), n_targets, n_sig_degs)) +
    theme_classic(base_size = 11) +
    theme(plot.subtitle = element_text(color = "grey30"),
          legend.position = "right")

  ggsave(outpath, p, width = 7.8, height = 5.6, dpi = 300)

  cat(sprintf("  %-6s  NES = %+.2f  padj = %s  -  %d / %d targets significant\n",
              tf_name, nes, fmt_p(padj), n_sig_degs, n_targets))
}


# ---- 6. GENERATE THREE REGULON PLOTS (Figs 9a, 9b, 9c) -------
cat("\nBuilding VIPER regulon plots (n = 45):\n")
make_regulon_plot("RELA",  regulons, de_df, viper_df,
                  file.path(OUT_FIG, "Figure9a_VIPER_RELA_regulon_n45.png"))
make_regulon_plot("NR2F1", regulons, de_df, viper_df,
                  file.path(OUT_FIG, "Figure9b_VIPER_NR2F1_regulon_n45.png"))
make_regulon_plot("REST",  regulons, de_df, viper_df,
                  file.path(OUT_FIG, "Figure9c_VIPER_REST_regulon_n45.png"))


# ---- 7. DONE ---------------------------------------------------
cat("\n================ STEP 6 COMPLETE ================\n")
cat("Four new figures written to:\n")
cat("  ", OUT_FIG, "\n")
cat("Files:\n")
cat("   Figure8_sample_correlation_heatmap_n54.png\n")
cat("   Figure9a_VIPER_RELA_regulon_n45.png\n")
cat("   Figure9b_VIPER_NR2F1_regulon_n45.png\n")
cat("   Figure9c_VIPER_REST_regulon_n45.png\n")
cat("================================================\n")
