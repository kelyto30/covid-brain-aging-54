# ================================================================
# Step 3 - Gene Set Enrichment Analysis (GO:BP)
# OUTLIERS-INCLUDED sensitivity analysis (n = 45)
# Dissertation: COVID-19 Brain Aging (GSE188847)
# Kelechi Wisdom Elechi | UT Health San Antonio
# ----------------------------------------------------------------
# Mirrors the dissertation GSEA:
#   ranking  : DESeq2 Wald statistic (unshrunken, from Step 2)
#   gene sets: GO Biological Process (org.Hs.eg.db)
#   engine   : clusterProfiler::gseGO (an fgsea wrapper)
#   perms    : 1000
#   minGSSize: 10
#   threshold: BH padj < 0.05
# Input : res_raw_df_n45.rds from Step 2
# Output: GSEA results CSV, dot plot (Figure 4), key-term comparison
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
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})


# ---- 2. LOAD RAW WALD RESULTS ----------------------------------
res_raw <- readRDS(file.path(OUT_RDS, "res_raw_df_n45.rds"))
cat("Raw DE results loaded:", nrow(res_raw), "genes\n")


# ---- 3. BUILD RANKED LIST (Wald stat, by symbol) --------------
r <- res_raw[!is.na(res_raw$stat) & !is.na(res_raw$gene_symbol) &
             nzchar(res_raw$gene_symbol), ]

# Deduplicate by symbol: keep the row with the largest |stat| per symbol
r <- r[order(-abs(r$stat)), ]
r <- r[!duplicated(r$gene_symbol), ]

# Sort by Wald stat descending (strongest positive at top)
r <- r[order(-r$stat), ]

ranks <- setNames(r$stat, r$gene_symbol)
cat("Ranked list:", length(ranks), "unique symbols\n")
cat("Top 3 (strongest up)   :", paste(names(head(ranks, 3)), collapse = ", "), "\n")
cat("Bottom 3 (strongest dn):", paste(names(tail(ranks, 3)), collapse = ", "), "\n")


# ---- 4. RUN GSEA (GO:BP) ---------------------------------------
set.seed(42)
cat("\nRunning clusterProfiler::gseGO (BP)... this takes a few minutes\n")
gsea_bp <- gseGO(geneList     = ranks,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "SYMBOL",
                 ont          = "BP",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 nPermSimple  = 1000,
                 eps          = 0,
                 seed         = TRUE,
                 verbose      = FALSE)

gsea_df <- as.data.frame(gsea_bp)
cat("\nSignificant GO:BP terms (padj < 0.05):", nrow(gsea_df), "\n")

# Tag direction
gsea_df$direction <- ifelse(gsea_df$NES > 0, "Up", "Down")
cat("  Up   (NES > 0):", sum(gsea_df$direction == "Up"),   "\n")
cat("  Down (NES < 0):", sum(gsea_df$direction == "Down"), "\n")


# ---- 5. SAVE ---------------------------------------------------
write.csv(gsea_df, file.path(OUT_TBL, "GSEA_GOBP_n45_outliers_included.csv"),
          row.names = FALSE)
saveRDS(gsea_bp, file.path(OUT_RDS, "gsea_bp_n45.rds"))


# ---- 6. DOT PLOT: TOP 30 (15 up + 15 down) --------------------
top_up   <- head(gsea_df[order(-gsea_df$NES), ], 15)
top_down <- head(gsea_df[order( gsea_df$NES), ], 15)
top30    <- rbind(top_up, top_down)
top30$Description <- factor(top30$Description,
                            levels = top30$Description[order(top30$NES)])

p_gsea <- ggplot(top30,
                 aes(x = NES, y = Description,
                     size  = setSize,
                     color = -log10(p.adjust))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_viridis_c(option = "C") +
  scale_size_continuous(range = c(2, 7)) +
  labs(x = "Normalized Enrichment Score",
       y = NULL,
       size  = "Gene set size",
       color = "-log10(padj)",
       title = "GSEA: Top 30 GO:BP terms (n = 45, outliers retained)",
       subtitle = sprintf("%d significant terms at padj < 0.05  |  %d up, %d down",
                          nrow(gsea_df),
                          sum(gsea_df$direction == "Up"),
                          sum(gsea_df$direction == "Down"))) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 8),
        plot.title.position = "plot")

ggsave(file.path(OUT_FIG, "Figure4_GSEA_dotplot_n45.png"),
       p_gsea, width = 9, height = 8.5, dpi = 300)


# ---- 7. AGING-RELEVANT TERMS (for comparison to dissertation) -
# Pull terms matching themes from the original Results section 4.2
patterns <- list(
  "Immune/IFN"  = "immune|interferon|cytokine|defense response to virus|antiviral|complement",
  "Synaptic"    = "synaptic|synapse|neurotransmitter|glutamate receptor",
  "Neuronal"    = "axon|neuron|membrane potential|action potential"
)

key_terms <- do.call(rbind, lapply(names(patterns), function(n) {
  hits <- gsea_df[grepl(patterns[[n]], gsea_df$Description, ignore.case = TRUE), ]
  if (nrow(hits) == 0) return(NULL)
  hits$category <- n
  hits[, c("category", "ID", "Description", "NES", "pvalue", "p.adjust", "setSize")]
}))
key_terms <- key_terms[order(key_terms$category, -abs(key_terms$NES)), ]

write.csv(key_terms,
          file.path(OUT_TBL, "GSEA_aging_relevant_terms_n45.csv"),
          row.names = FALSE)


# ---- 8. SPOT-CHECK KEY TERMS FROM THE DISSERTATION ------------
# Original reported: chemical synaptic transmission NES = -2.41
dissertation_terms <- c(
  "innate immune response",
  "defense response to virus",
  "interferon-beta signaling|type I interferon",
  "cytokine-mediated signaling",
  "complement activation",
  "chemical synaptic transmission",
  "regulation of membrane potential",
  "neurotransmitter receptor",
  "action potential",
  "axon development"
)

cat("\n--- Dissertation-level term spot-check ---\n")
spot <- data.frame(
  Term_queried = character(),
  Matched_term = character(),
  NES = numeric(),
  padj = numeric(),
  setSize = integer(),
  stringsAsFactors = FALSE
)
for (t in dissertation_terms) {
  hit <- gsea_df[grepl(t, gsea_df$Description, ignore.case = TRUE), ]
  if (nrow(hit) == 0) {
    spot <- rbind(spot, data.frame(Term_queried = t, Matched_term = "(not found)",
                                   NES = NA, padj = NA, setSize = NA))
  } else {
    hit <- hit[order(-abs(hit$NES)), ][1, ]
    spot <- rbind(spot, data.frame(
      Term_queried = t,
      Matched_term = hit$Description,
      NES          = round(hit$NES, 3),
      padj         = signif(hit$p.adjust, 3),
      setSize      = hit$setSize
    ))
  }
}
print(spot, row.names = FALSE)
write.csv(spot, file.path(OUT_TBL, "GSEA_dissertation_term_spotcheck_n45.csv"),
          row.names = FALSE)


# ---- 9. STEP 3 SUMMARY -----------------------------------------
ORIG_TOTAL_TERMS <- 1025  # from dissertation Table 3

cat("\n================= STEP 3 SUMMARY =================\n")
cat(sprintf("%-32s %12s %12s %+10s\n", "Metric",
            "Original n=42", "New n=45", "Change"))
cat(sprintf("%-32s %12d %12d %+10d\n", "Significant GO:BP terms",
            ORIG_TOTAL_TERMS, nrow(gsea_df),
            nrow(gsea_df) - ORIG_TOTAL_TERMS))
cat(sprintf("%-32s %12s %12d %10s\n", "  Upregulated terms",
            "(not reported)", sum(gsea_df$direction == "Up"), ""))
cat(sprintf("%-32s %12s %12d %10s\n", "  Downregulated terms",
            "(not reported)", sum(gsea_df$direction == "Down"), ""))
cat("===================================================\n")
cat("Outputs written to:", OUT, "\n")
