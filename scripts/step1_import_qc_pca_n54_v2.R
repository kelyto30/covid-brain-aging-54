# ================================================================
# Step 1 — Data import, QC, and diagnostic PCA (n = 54)
# OUTLIERS-INCLUDED sensitivity analysis
# Dissertation: COVID-19 Brain Aging (GSE188847 frontal cortex)
# Kelechi Wisdom Elechi | UT Health San Antonio
# ----------------------------------------------------------------
# Re-runs import/QC/PCA on ALL 54 frontal-cortex samples,
# retaining COVID11/14/15 (removed in original n=51 analysis).
# PCA is diagnostic, not a filter. Library-size QC is retained.
#
# Outputs: <PROJ>/outliers_included/{figures,tables,rdata}
# ================================================================


# ---- 0. PATHS (loaded from config.R at repo root) -------------
# Loads PROJECT_ROOT, QUANTS_DIR, META_PATH, TX2GENE_PATH, and
# the OUT_* directory variables from the central config file.
# Run from the repo root (e.g. `Rscript scripts/step1_...R`) or
# set the working directory to the repo root before sourcing.
if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists("../config.R")) {
  source("../config.R")
} else {
  stop("Cannot find config.R. Set your working directory to the ",
       "repository root before running this script.")
}


# ---- 1. OUTPUT DIRECTORIES (already created by config.R) ------
# (PROJECT_ROOT, QUANTS_DIR, META_PATH, TX2GENE_PATH, OUT, OUT_FIG,
#  OUT_TBL, OUT_RDS are all defined by config.R)


# ---- 2. LIBRARIES ---------------------------------------------
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})


# ---- 3. METADATA (filter to the 54-sample disease cohort) -----
meta <- read.csv(META_PATH, stringsAsFactors = FALSE)
cat("Full metadata:", nrow(meta), "rows\n")

meta <- subset(meta, group %in% c("COVID-19", "Control", "ICU/VENT"))
cat("After filtering to frontal-cortex disease cohort:", nrow(meta), "rows\n")
cat("Group counts (expected: COVID-19=22, Control=23, ICU/VENT=9):\n")
print(table(meta$group))

if (nrow(meta) != 54) {
  stop("Expected 54 rows after filtering; got ", nrow(meta),
       ". Check the 'group' values in the metadata file.")
}


# ---- 4. MATCH EACH GSM TO ITS .gz FILE ------------------------
all_files <- list.files(QUANTS_DIR, pattern = "\\.gz$", full.names = FALSE)

meta$file_path <- NA_character_
meta$sample_id <- NA_character_

for (i in seq_len(nrow(meta))) {
  gsm <- meta$geo_accession[i]
  hits <- all_files[grepl(paste0("^", gsm, "_"), all_files)]
  if (length(hits) != 1) {
    stop("Expected 1 file for ", gsm, "; found ", length(hits))
  }
  meta$file_path[i] <- file.path(QUANTS_DIR, hits)
  # Extract the human-readable name between GSMxxxx_ and .salmon
  meta$sample_id[i] <- sub("^GSM\\d+_", "", sub("\\.salmon.*$", "", hits))
}

stopifnot(all(file.exists(meta$file_path)))

# Normalize ICU label: filenames use "ICUVENT", dissertation uses "ICU/VENT" for group
# The sample_id keeps the filename form (ICUVENT1..9) which is fine for labeling.
cat("\nSample IDs parsed from filenames (first few):\n")
print(head(meta$sample_id, 10))
cat("Originally-flagged outliers present?:\n")
print(c("COVID11", "COVID14", "COVID15") %in% meta$sample_id)

# Flag for plotting (NOT a filter)
ORIG_OUTLIERS <- c("COVID11", "COVID14", "COVID15")
meta$pca_outlier_orig <- meta$sample_id %in% ORIG_OUTLIERS


# ---- 5. tximport ----------------------------------------------
files <- setNames(meta$file_path, meta$sample_id)

tx2gene <- read.csv(TX2GENE_PATH, stringsAsFactors = FALSE)
cat("\ntx2gene loaded:", nrow(tx2gene), "transcripts\n")

cat("Running tximport (this takes ~1-2 min for 54 samples)...\n")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")])
cat("Gene-level count matrix:", nrow(txi$counts), "genes x", ncol(txi$counts), "samples\n")


# ---- 6. LIBRARY-SIZE QC ---------------------------------------
lib_M <- colSums(txi$counts) / 1e6
lib_df <- data.frame(
  sample       = names(lib_M),
  lib_size_M   = as.numeric(lib_M),
  group        = meta$group[match(names(lib_M), meta$sample_id)],
  outlier_orig = meta$pca_outlier_orig[match(names(lib_M), meta$sample_id)]
)
lib_df <- lib_df[order(lib_df$lib_size_M), ]

cat("\n--- Library size (M reads) ---\n")
print(summary(lib_df$lib_size_M))
cat("Samples below 10M threshold:", sum(lib_df$lib_size_M < 10), "\n")

write.csv(lib_df, file.path(OUT_TBL, "library_sizes_n54.csv"), row.names = FALSE)

p_lib <- ggplot(lib_df,
                aes(x = reorder(sample, lib_size_M), y = lib_size_M, fill = group)) +
  geom_col() +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  geom_text(aes(label = ifelse(outlier_orig, sample, "")),
            hjust = -0.1, size = 2.8) +
  coord_flip() +
  scale_fill_manual(values = c("COVID-19" = "#E74C3C",
                               "Control"  = "#3498DB",
                               "ICU/VENT" = "#27AE60")) +
  labs(x = NULL, y = "Library size (millions of mapped reads)",
       title = "Library size distribution — n = 54 (outliers retained)",
       fill = "Group") +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(OUT_FIG, "Figure2_library_sizes_n54.png"),
       p_lib, width = 7, height = 9, dpi = 300)


# ---- 7. DESeq2 OBJECT (all 54) --------------------------------
meta <- meta[match(colnames(txi$counts), meta$sample_id), ]
rownames(meta) <- meta$sample_id
meta$group <- factor(meta$group, levels = c("Control", "COVID-19", "ICU/VENT"))
meta$sex   <- factor(meta$sex)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ age + sex + group)

# Pre-filter low-count genes (matching original pipeline)
keep <- rowSums(counts(dds) >= 10) >= 5
dds  <- dds[keep, ]
cat("\nGenes after pre-filter (>=10 counts in >=5 samples):", nrow(dds), "\n")


# ---- 8. VST + DIAGNOSTIC PCA ----------------------------------
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE, ntop = 500)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_data$sample       <- rownames(pca_data)
pca_data$outlier_orig <- pca_data$sample %in% ORIG_OUTLIERS

# Re-apply the original ">3 SD from group centroid" rule for reporting only
pca_data$dist_sd <- NA_real_
for (g in unique(pca_data$group)) {
  idx <- pca_data$group == g
  cx  <- mean(pca_data$PC1[idx]); cy <- mean(pca_data$PC2[idx])
  sx  <- sd(pca_data$PC1[idx]);   sy <- sd(pca_data$PC2[idx])
  pca_data$dist_sd[idx] <- sqrt(((pca_data$PC1[idx] - cx)/sx)^2 +
                                ((pca_data$PC2[idx] - cy)/sy)^2)
}
pca_data$flagged_by_3sd_rule <- pca_data$dist_sd > 3

write.csv(pca_data, file.path(OUT_TBL, "pca_coordinates_n54.csv"), row.names = FALSE)

cat("\n--- PCA diagnostic (no filtering applied) ---\n")
cat("PC1 variance:", percentVar[1], "% | PC2 variance:", percentVar[2], "%\n")
cat("Samples >3 SD from group centroid (retained anyway):\n")
print(pca_data[pca_data$flagged_by_3sd_rule,
               c("sample", "group", "PC1", "PC2", "dist_sd")])
cat("\nCoordinates of originally-removed samples:\n")
print(pca_data[pca_data$outlier_orig,
               c("sample", "group", "PC1", "PC2", "dist_sd", "flagged_by_3sd_rule")])

p_pca <- ggplot(pca_data,
                aes(x = PC1, y = PC2, color = group, shape = outlier_orig)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(aes(label = ifelse(outlier_orig | flagged_by_3sd_rule,
                                     sample, "")),
                  size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Control"  = "#3498DB",
                                "COVID-19" = "#E74C3C",
                                "ICU/VENT" = "#27AE60")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "Retained in original QC",
                                `TRUE`  = "Removed in original QC"),
                     name = "Original QC status") +
  labs(x = paste0("PC1 (", percentVar[1], "%)"),
       y = paste0("PC2 (", percentVar[2], "%)"),
       title = "Diagnostic PCA - all 54 samples retained",
       subtitle = "Originally-flagged outliers (COVID11/14/15) shown as triangles",
       color = "Group") +
  theme_classic(base_size = 11)

ggsave(file.path(OUT_FIG, "Figure1_PCA_n54.png"),
       p_pca, width = 8, height = 6, dpi = 300)


# ---- 9. SAVE OBJECTS FOR STEP 2 -------------------------------
saveRDS(dds,      file.path(OUT_RDS, "dds_n54_all_groups.rds"))
saveRDS(vsd,      file.path(OUT_RDS, "vsd_n54.rds"))
saveRDS(meta,     file.path(OUT_RDS, "meta_n54.rds"))
saveRDS(pca_data, file.path(OUT_RDS, "pca_data_n54.rds"))


# ---- 10. STEP 1 SUMMARY ---------------------------------------
cat("\n================= STEP 1 SUMMARY =================\n")
cat("Total samples         :", ncol(dds), "\n")
cat("By group              :\n"); print(table(meta$group))
cat("Genes after pre-filter:", nrow(dds), "\n")
cat("Library size range (M):",
    round(min(lib_df$lib_size_M), 1), "-",
    round(max(lib_df$lib_size_M), 1), "\n")
cat("Samples <10M reads    :", sum(lib_df$lib_size_M < 10), "\n")
cat("PC1 variance          :", percentVar[1], "%\n")
cat("PC2 variance          :", percentVar[2], "%\n")
cat("Outputs written to    :", OUT, "\n")
cat("===================================================\n")
