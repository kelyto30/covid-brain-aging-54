# Data

This folder contains two reference files used by the analysis pipeline,
plus instructions for obtaining the raw transcriptomic data which is
too large to commit to the repository.

## Files included

### `GSE188847_metadata_clean.csv`
Cleaned and harmonized sample metadata for all 54 frontal cortex donors
in the Mavrikaki et al. cohort (GEO accession GSE188847). Columns:

- `geo_accession`: GSM identifier matching the .gz filename prefix
- `age`: donor age in years at death
- `group`: one of `COVID-19`, `Control`, `ICU/VENT`
- `sex`: `M` or `F`
- `source`: tissue source (Frontal cortex of brain)
- `sex_male`: indicator (1 = male, 0 = female)
- `age_group_60`: indicator (`<60` or `>=60`)
- `diagnosis`: collapsed diagnostic label

### `tx2gene.csv` (24 MB)
Ensembl v104 (GRCh38) transcript-to-gene mapping table used by
`tximport` to summarize transcript-level Salmon quantifications to
gene-level counts. Columns:

- `TXNAME`: Ensembl transcript ID with version
- `GENEID`: Ensembl gene ID with version
- `GENENAME`: HGNC gene symbol

This file is bundled with the repository to avoid forcing users to
regenerate it. It corresponds to the Ensembl v104 release available
at <https://ftp.ensembl.org/pub/release-104/>.

---

## Salmon quantifications (NOT included)

The 54 Salmon-quantified `.gz` files (one per donor) are too large to
commit to a Git repository. Total size after download is
approximately 1–2 GB.

### Option A: Download from GEO
1. Visit the GEO accession page:
   <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188847>
2. In the **Supplementary file** section, download
   `GSE188847_RAW.tar` (contains all 54 sample files).
3. Extract the tar archive into a folder of your choice, for example
   `~/GSE188847_quants/`.
4. Edit `config.R` at the repo root and set `QUANTS_DIR` to the
   absolute path of that folder.

### Option B: Download via GEOquery in R
```r
# Requires: BiocManager::install("GEOquery")
library(GEOquery)
getGEOSuppFiles("GSE188847", baseDir = "~/")
# Then untar GSE188847_RAW.tar in ~/GSE188847/ and point QUANTS_DIR there.
```

### Verification

After download, the folder should contain 54 files matching the
pattern `GSM569XXXX_*.salmon.sf.gz`. The first script (`step1_…R`)
will fail with an informative error if any expected sample file is
missing.

---

## Original data citation

Mavrikaki M, Lee JD, Solomon IH, Slack FJ.
**Severe COVID-19 is associated with molecular signatures of aging
in the human brain.** *Nature Aging* 2:1130–1137 (2022).
<https://doi.org/10.1038/s43587-022-00321-w>

Original investigators retain all rights to the underlying data.
This repository contains analytical code only; raw transcriptomic
data must be obtained from GEO under the original investigators'
deposition terms.
