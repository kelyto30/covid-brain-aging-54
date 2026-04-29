## ============================================================
##  step7d_volcano_labels_v2.R
##  CORRECTED volcano regeneration with gene labels.
##
##  Issue with v1 (step7c): the script symbol-mapped the FULL
##  DESeq2 result before plotting, dropping ~6% of genes with no
##  symbol annotation. The DEG counts in the figure (3,715 up /
##  2,888 down / 6,603 total) drifted from the manuscript-anchored
##  Ensembl-keyed counts (4,086 up / 2,993 down / 7,079 total).
##
##  Fix: keep the full Ensembl-keyed DESeq2 result for the plot
##  cloud and counts. Use symbol→Ensembl REVERSE mapping only for
##  the 9 genes we want to label, find them by Ensembl ID in the
##  full result, and place labels at their actual coordinates.
##
##  Output overwrites:
##    outliers_included/figures/Figure3_volcano_n45.png
##  and matches the manuscript text exactly:
##    7,079 DEGs at padj < 0.05  |  4,086 up  |  2,993 down
##
##  Run on local Mac. ~ 30 seconds.
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
stopifnot(dir.exists(rd_dir), dir.exists(fg_dir))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
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

cat("Loading res_apeglm_n45.rds ...\n")
res_apeglm <- readRDS(file.path(rd_dir, "res_apeglm_n45.rds"))

## ---- Build the FULL Ensembl-keyed plot dataframe ----------
## (keep every gene; NO symbol-based filter here)
res_df <- as.data.frame(res_apeglm)
res_df$ensembl <- rownames(res_df)

# Strip any version suffix on ensembl IDs for clean matching later
res_df$ensembl_clean <- sub("\\..*$", "", res_df$ensembl)

# Drop only rows with NA padj (these can't be plotted anyway)
res_df <- res_df %>% filter(!is.na(padj))

# Build status column from Ensembl-keyed full result
res_df <- res_df %>%
  mutate(status = case_when(
    padj < 0.05 & log2FoldChange >  0 ~ "Up",
    padj < 0.05 & log2FoldChange <  0 ~ "Down",
    TRUE                               ~ "NS"
  ))

n_up    <- sum(res_df$status == "Up")
n_down  <- sum(res_df$status == "Down")
n_total <- n_up + n_down

cat(sprintf("  Plot cloud: %d genes total\n", nrow(res_df)))
cat(sprintf("  DEGs at padj < 0.05: %d up, %d down, %d total\n",
            n_up, n_down, n_total))

## ---- Reverse map: SYMBOL -> ENSEMBL for ONLY the 9 labels --
label_up   <- c("MT2A", "SOCS3", "CXCL10", "MMP9", "BCL2A1")
label_down <- c("GRIN3A", "KIF21B", "BDNF", "LINGO1")
label_genes <- c(label_up, label_down)

cat("\nReverse-mapping label gene symbols to Ensembl IDs ...\n")
sym_to_ens <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = label_genes,
  column    = "ENSEMBL",
  keytype   = "SYMBOL",
  multiVals = "first"
)
cat("Symbol -> Ensembl:\n")
print(sym_to_ens)

# Detach AnnotationDbi to prevent dplyr namespace masking downstream
suppressWarnings({
  try(detach("package:org.Hs.eg.db",  unload = TRUE, character.only = TRUE), silent = TRUE)
  try(detach("package:AnnotationDbi", unload = TRUE, character.only = TRUE), silent = TRUE)
})
library(dplyr)

## ---- Build label dataframe by Ensembl-ID lookup ----------
lab_lookup <- data.frame(
  gene    = names(sym_to_ens),
  ensembl_clean = unname(sym_to_ens),
  stringsAsFactors = FALSE
) %>% filter(!is.na(ensembl_clean))

# Join to the full plot dataframe via cleaned Ensembl ID
lab_df <- lab_lookup %>%
  inner_join(
    res_df %>% dplyr::select(ensembl_clean, log2FoldChange, padj, status),
    by = "ensembl_clean"
  )

cat(sprintf("\nLabeled genes resolved (%d of %d):\n",
            nrow(lab_df), length(label_genes)))
print(lab_df)

missing <- setdiff(label_genes, lab_df$gene)
if (length(missing) > 0) {
  cat("WARNING — these genes could not be resolved (will not be labeled):\n")
  cat("  ", paste(missing, collapse = ", "), "\n")
}

## ---- Plot ------------------------------------------------
g <- ggplot(res_df,
            aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.55, size = 0.9) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = lab_df,
    aes(label = gene),
    color           = "black",
    size            = 4.0,
    fontface        = "bold",
    max.overlaps    = Inf,
    box.padding     = 0.6,
    point.padding   = 0.3,
    segment.color   = "grey40",
    segment.size    = 0.4,
    min.segment.length = 0,
    force           = 2,
    seed            = 42
  ) +
  scale_color_manual(
    values = c(Up = "#C0392B", Down = "#2E86C1", NS = "grey70"),
    labels = c(
      Up   = sprintf("Up (n=%d)",   n_up),
      Down = sprintf("Down (n=%d)", n_down),
      NS   = "NS"
    )
  ) +
  labs(
    title    = "Differential expression: COVID-19 vs Control",
    subtitle = sprintf(
      "%d DEGs at padj < 0.05  |  design: ~ age + sex + group  (n = 45)",
      n_total),
    x = "Log2 fold change (COVID-19 / Control) — apeglm shrunken",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  neutral_theme() +
  theme(legend.position = "right")

out_file <- file.path(fg_dir, "Figure3_volcano_n45.png")
ggsave(out_file, plot = g, width = 8, height = 5.5, dpi = 300, bg = "white")

cat("\n================================================\n")
cat("STEP 7d COMPLETE — corrected labeled volcano regenerated.\n")
cat(sprintf("  Subtitle: %d DEGs (Up=%d, Down=%d)\n",
            n_total, n_up, n_down))
cat("File: ", out_file, "\n", sep = "")
cat("================================================\n")
