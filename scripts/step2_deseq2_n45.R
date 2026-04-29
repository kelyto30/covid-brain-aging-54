# ================================================================
# Step 2 - Differential expression (DESeq2): COVID-19 vs Control
# OUTLIERS-INCLUDED sensitivity analysis (n = 45: 22 COVID + 23 CTRL)
# Dissertation: COVID-19 Brain Aging (GSE188847)
# Kelechi Wisdom Elechi | UT Health San Antonio
# ----------------------------------------------------------------
# Mirrors the original dissertation model:
#   design = ~ age + sex + group
#   shrinkage = apeglm
#   significance = BH padj < 0.05
# Input: dds_n54_all_groups.rds from Step 1
# Output: DEG table, volcano plot (Figure 3), objects for Step 3 (GSEA)
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
  library(ggplot2)
  library(ggrepel)
})


# ---- 2. LOAD STEP 1 OBJECTS ------------------------------------
dds     <- readRDS(file.path(OUT_RDS, "dds_n54_all_groups.rds"))
meta    <- readRDS(file.path(OUT_RDS, "meta_n54.rds"))
tx2gene <- read.csv(TX2GENE_PATH, stringsAsFactors = FALSE)

cat("Loaded dds_n54:", nrow(dds), "genes x", ncol(dds), "samples\n")


# ---- 3. SUBSET TO COVID-19 vs CONTROL (drop ICU/VENT) ---------
keep <- dds$group %in% c("Control", "COVID-19")
dds_cc <- dds[, keep]
dds_cc$group <- droplevels(dds_cc$group)

cat("\nAnalytic subset for primary DE comparison:\n")
print(table(dds_cc$group))
stopifnot(ncol(dds_cc) == 45)


# ---- 4. RUN DESeq2 ---------------------------------------------
cat("\nRunning DESeq() ...\n")
dds_cc <- DESeq(dds_cc)

cat("\nModel coefficients:\n")
print(resultsNames(dds_cc))

# Identify the COVID-vs-Control coefficient (DESeq2 rewrites the hyphen in "COVID-19")
coef_name <- grep("COVID.*Control", resultsNames(dds_cc), value = TRUE)
if (length(coef_name) != 1) {
  stop("Could not uniquely identify the COVID-vs-Control coefficient. Got: ",
       paste(coef_name, collapse = ", "))
}
cat("Using coefficient:", coef_name, "\n")


# ---- 5. EXTRACT RESULTS: RAW + SHRUNKEN ------------------------
# Raw Wald stats (for GSEA ranking in Step 3)
res_raw <- results(dds_cc, name = coef_name, alpha = 0.05)
cat("\n--- Unshrunken results summary ---\n"); summary(res_raw)

# Apeglm-shrunken (for effect sizes and volcano)
res_shr <- lfcShrink(dds_cc, coef = coef_name, type = "apeglm")
cat("\n--- Apeglm-shrunken results summary ---\n"); summary(res_shr)


# ---- 6. ANNOTATE WITH GENE SYMBOLS -----------------------------
res_df <- as.data.frame(res_shr)
res_df$gene_id <- rownames(res_df)
res_df$gene_id_clean <- sub("\\.[0-9]+$", "", res_df$gene_id)

tx2gene_u <- unique(tx2gene[, c("GENEID", "GENENAME")])
tx2gene_u$GENEID_clean <- sub("\\.[0-9]+$", "", tx2gene_u$GENEID)
res_df$gene_symbol <- tx2gene_u$GENENAME[match(res_df$gene_id_clean, tx2gene_u$GENEID_clean)]

# Also annotate the raw (Wald stat) results for GSEA
res_raw_df <- as.data.frame(res_raw)
res_raw_df$gene_id <- rownames(res_raw_df)
res_raw_df$gene_id_clean <- sub("\\.[0-9]+$", "", res_raw_df$gene_id)
res_raw_df$gene_symbol <- tx2gene_u$GENENAME[match(res_raw_df$gene_id_clean,
                                                   tx2gene_u$GENEID_clean)]

# Reorder and sort
res_df <- res_df[, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange",
                     "lfcSE", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj), ]

res_raw_df <- res_raw_df[, c("gene_id", "gene_symbol", "baseMean",
                             "log2FoldChange", "lfcSE", "stat",
                             "pvalue", "padj")]


# ---- 7. SAVE TABLES AND OBJECTS --------------------------------
write.csv(res_df,
          file.path(OUT_TBL, "DEG_results_n45_apeglm.csv"),
          row.names = FALSE)
write.csv(res_raw_df,
          file.path(OUT_TBL, "DEG_results_n45_raw_Wald.csv"),
          row.names = FALSE)

saveRDS(dds_cc,     file.path(OUT_RDS, "dds_cc_n45.rds"))
saveRDS(res_shr,    file.path(OUT_RDS, "res_apeglm_n45.rds"))
saveRDS(res_raw,    file.path(OUT_RDS, "res_raw_n45.rds"))
saveRDS(res_raw_df, file.path(OUT_RDS, "res_raw_df_n45.rds"))


# ---- 8. COUNT DEGs AND COMPARE TO ORIGINAL ---------------------
sig      <- !is.na(res_df$padj) & res_df$padj < 0.05
n_total  <- sum(sig)
n_up     <- sum(sig & res_df$log2FoldChange > 0)
n_down   <- sum(sig & res_df$log2FoldChange < 0)

# Original dissertation numbers (from Table 3)
ORIG_TOTAL <- 5637; ORIG_UP <- 2869; ORIG_DOWN <- 2768

cat("\n================= STEP 2 SUMMARY =================\n")
cat(sprintf("%-26s %12s %12s %12s\n", "Metric", "Original n=42", "New n=45", "Change"))
cat(sprintf("%-26s %12d %12d %+12d\n", "Total DEGs (padj<0.05)",
            ORIG_TOTAL, n_total, n_total - ORIG_TOTAL))
cat(sprintf("%-26s %12d %12d %+12d\n", "Upregulated",
            ORIG_UP,    n_up,    n_up    - ORIG_UP))
cat(sprintf("%-26s %12d %12d %+12d\n", "Downregulated",
            ORIG_DOWN,  n_down,  n_down  - ORIG_DOWN))
cat("===================================================\n\n")


# ---- 9. VOLCANO PLOT (Figure 3 equivalent) --------------------
vdf <- res_df
vdf$sig <- "NS"
vdf$sig[sig & vdf$log2FoldChange > 0] <- "Up"
vdf$sig[sig & vdf$log2FoldChange < 0] <- "Down"
vdf$neglog10p <- -log10(vdf$padj)

# Label the same key genes used in the original dissertation volcano
key_genes <- c("MT2A", "MT1E", "STAT1", "CXCL10", "GRIN3A", "SOCS3",
               "KIF21B", "BDNF", "BEGAIN", "MMP9", "BCL2A1", "CPNE7", "LINGO1")
vdf$label <- ifelse(vdf$gene_symbol %in% key_genes & sig, vdf$gene_symbol, "")

p_vol <- ggplot(vdf, aes(x = log2FoldChange, y = neglog10p, color = sig)) +
  geom_point(alpha = 0.55, size = 0.8) +
  geom_text_repel(aes(label = label), size = 3, color = "black",
                  max.overlaps = 30, box.padding = 0.35,
                  min.segment.length = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("NS"   = "grey75",
                                "Up"   = "#E74C3C",
                                "Down" = "#3498DB"),
                     breaks = c("Up", "Down", "NS"),
                     labels = c(paste0("Up (n=",   n_up, ")"),
                                paste0("Down (n=", n_down, ")"),
                                "NS")) +
  labs(x = "Log2 fold change (COVID-19 / Control) - apeglm shrunken",
       y = "-log10(adj. p-value)",
       title = "Volcano plot - DESeq2 (n = 45, outliers retained)",
       subtitle = sprintf("%d DEGs at padj < 0.05  |  age + sex + group",
                          n_total),
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(OUT_FIG, "Figure3_volcano_n45.png"),
       p_vol, width = 8, height = 6.5, dpi = 300)


# ---- 10. KEY GENE SPOT-CHECK ----------------------------------
cat("\n--- Key gene effect sizes (for comparison to dissertation) ---\n")
spot <- res_df[res_df$gene_symbol %in% key_genes, ]
spot <- spot[order(spot$padj), ]
print(spot, row.names = FALSE)

cat("\nOutputs written to:", OUT, "\n")
