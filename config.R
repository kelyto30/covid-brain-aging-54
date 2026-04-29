# ============================================================
# config.R — paths and global settings for the analysis pipeline
# ============================================================
# This file is sourced by every script in /scripts/.
# Edit the four paths below to match your local environment.
#
# Only PROJECT_ROOT and QUANTS_DIR need to be changed by users
# cloning this repository. The other two paths are derived
# automatically.
# ============================================================


# ---- 1. PROJECT ROOT ---------------------------------------
# The folder containing this repository (the cloned `covid-brain-aging-thesis`).
# Change this to the absolute path of your clone.
PROJECT_ROOT <- "/Users/yourname/path/to/covid-brain-aging-thesis"


# ---- 2. SALMON QUANTIFICATION DIRECTORY --------------------
# Folder containing the 54 Salmon-quantified .gz files
# downloaded from GEO accession GSE188847.
# See data/README.md for download instructions.
QUANTS_DIR <- "/Users/yourname/path/to/GSE188847_quants"


# ---- 3. METADATA AND TX2GENE PATHS -------------------------
# These come bundled with the repository under /data/ and do
# not normally need to be changed.
META_PATH      <- file.path(PROJECT_ROOT, "data", "GSE188847_metadata_clean.csv")
TX2GENE_PATH   <- file.path(PROJECT_ROOT, "data", "tx2gene.csv")


# ---- 4. OUTPUT DIRECTORIES ---------------------------------
# All scripts write to outliers_included/{figures,tables,rdata}
# inside PROJECT_ROOT. These are created on first run.
OUT     <- file.path(PROJECT_ROOT, "outliers_included")
OUT_FIG <- file.path(OUT, "figures")
OUT_TBL <- file.path(OUT, "tables")
OUT_RDS <- file.path(OUT, "rdata")

for (d in c(OUT, OUT_FIG, OUT_TBL, OUT_RDS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ---- 5. SANITY CHECKS --------------------------------------
if (!dir.exists(PROJECT_ROOT)) {
  stop("PROJECT_ROOT does not exist: ", PROJECT_ROOT,
       "\nEdit config.R and set PROJECT_ROOT to your local clone path.")
}
if (!file.exists(META_PATH)) {
  stop("Metadata file not found at: ", META_PATH,
       "\nExpected to be bundled in the repository at /data/.")
}
if (!file.exists(TX2GENE_PATH)) {
  stop("tx2gene file not found at: ", TX2GENE_PATH,
       "\nExpected to be bundled in the repository at /data/.")
}
# QUANTS_DIR is checked only by step1 (the only script that needs it).
