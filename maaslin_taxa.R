# =============================================================================
# MaAsLin2 — Differential abundance on MetaPhlAn taxonomic features
#
# Runs MaAsLin2 to identify taxa associated with a clinical outcome,
# controlling for covariates. Input is a species/genus-level relative
# abundance table from MetaPhlAn (samples × taxa) and a metadata table.
#
# Usage: Rscript maaslin_taxa.R
# =============================================================================

library(Maaslin2)
library(dplyr)

# =============================================================================
# CONFIG — set BASE_DIR to your project root, then adjust filenames if needed
# =============================================================================

BASE_DIR <- "/path/to/your/project"   # <-- set this to your project directory

FEATURES_FILE  <- file.path(BASE_DIR, "metaphlan_features.tsv")  # samples as rows, taxa as columns
METADATA_FILE  <- file.path(BASE_DIR, "metadata.tsv")             # samples as rows, covariates as columns
OUTPUT_DIR     <- file.path(BASE_DIR, "output", "maaslin_taxa")   # will be created if it doesn't exist

# Column names in metadata
SAMPLE_ID_COL  <- "sampleID"        # column with sample identifiers
TARGET_COL     <- "symptoms"        # main outcome variable
REFERENCE_LEVEL <- "Control"        # reference level for the outcome

# Covariates to include as fixed effects (must be columns in metadata)
FIXED_EFFECTS  <- c("symptoms", "mode_of_delivery", "visit_age_mo", "probiotics_firstyr")

# Optional: filter samples by age window (set to NULL to use all samples)
MAX_AGE_MONTHS <- 2       # keep only samples with visit_age_mo <= this value
EXCLUDE_GROUP  <- "Resolved"  # exclude this group from the outcome column (set to NULL to skip)

# MaAsLin2 filtering thresholds
MIN_ABUNDANCE  <- 0.1     # minimum mean relative abundance
MIN_PREVALENCE <- 0.1     # minimum fraction of samples where feature is present

# =============================================================================
# LOAD DATA
# =============================================================================

features <- read.table(FEATURES_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.table(METADATA_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)

# =============================================================================
# FILTER SAMPLES
# =============================================================================

if (!is.null(MAX_AGE_MONTHS)) {
  metadata <- metadata %>% filter(visit_age_mo <= MAX_AGE_MONTHS)
}

if (!is.null(EXCLUDE_GROUP)) {
  metadata <- metadata %>% filter(.data[[TARGET_COL]] != EXCLUDE_GROUP)
}

# Align features to filtered metadata
features <- features[features[[SAMPLE_ID_COL]] %in% metadata[[SAMPLE_ID_COL]], ]

# Sort both by sample ID to ensure alignment
features <- features[order(features[[SAMPLE_ID_COL]]), ]
metadata <- metadata[order(metadata[[SAMPLE_ID_COL]]), ]

if (!all(features[[SAMPLE_ID_COL]] == metadata[[SAMPLE_ID_COL]])) {
  stop("Sample IDs in features and metadata do not match after filtering.")
}

# Set sample IDs as row names
rownames(features) <- features[[SAMPLE_ID_COL]]
rownames(metadata) <- metadata[[SAMPLE_ID_COL]]
features[[SAMPLE_ID_COL]] <- NULL
metadata[[SAMPLE_ID_COL]] <- NULL

cat("Samples after filtering:", nrow(features), "\n")
cat("Features:", ncol(features), "\n")

# =============================================================================
# SET FACTOR LEVELS
# =============================================================================

metadata[[TARGET_COL]] <- factor(metadata[[TARGET_COL]])
metadata[[TARGET_COL]] <- relevel(metadata[[TARGET_COL]], ref = REFERENCE_LEVEL)

# =============================================================================
# RUN MaAsLin2
# =============================================================================

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

fit <- Maaslin2(
  input_data     = features,
  input_metadata = metadata,
  output         = OUTPUT_DIR,
  fixed_effects  = FIXED_EFFECTS,
  reference      = paste0(TARGET_COL, ",", REFERENCE_LEVEL),
  transform      = "AST",
  min_abundance  = MIN_ABUNDANCE,
  min_prevalence = MIN_PREVALENCE
)

cat("MaAsLin2 completed. Results saved to:", OUTPUT_DIR, "\n")
