# =============================================================================
# MaAsLin2 — Differential abundance on HUMAnN pathway tables
#
# Runs MaAsLin2 on one or more HUMAnN pathway/gene-family tables,
# against one or more metadata subgroups (e.g. Resolved, Symptomatic).
# Each combination of pathway file × metadata group is run separately
# and results are saved to its own output subdirectory.
#
# Input feature files are expected in HUMAnN merged-table format:
#   rows = pathways/gene families, columns = samples
#   (the script transposes them internally)
#
# Usage: Rscript maaslin_humann.R
# =============================================================================

library(Maaslin2)
library(dplyr)

# =============================================================================
# CONFIG — set BASE_DIR to your project root, then adjust filenames if needed
# =============================================================================

BASE_DIR <-   # <-- fill it with path to the dir to the samples

# HUMAnN feature tables (one or more; can be pathway, gene family, etc.)
# Each file will be run against each metadata group below.
FEATURE_FILES <- c(
  file.path(BASE_DIR, "humann", "butyrate_pathways.csv"),
  file.path(BASE_DIR, "humann", "scfa_pathways.csv"),
  file.path(BASE_DIR, "humann", "succinate_pathways.csv")
)

# Metadata files — one per subgroup (e.g. clinical state, time point)
# Names will be used to label output subdirectories.
METADATA_FILES <- list(
  Resolved    = file.path(BASE_DIR, "metadata", "metadata_Resolved.tsv"),
  PreSymptoms = file.path(BASE_DIR, "metadata", "metadata_PreSymptoms.tsv"),
  Symptomatic = file.path(BASE_DIR, "metadata", "metadata_Symptomatic.tsv")
)

# Root output directory — subdirs will be created automatically
OUTPUT_ROOT <- file.path(BASE_DIR, "output", "maaslin_humann")

# Column names
SAMPLE_ID_COL   <- "sampleID"
TARGET_COL      <- "symptoms"
REFERENCE_LEVEL <- "Control"

# Covariates to include as fixed effects
FIXED_EFFECTS <- c("symptoms", "mode_of_delivery", "visit_age_mo", "probiotics_firstyr")

# MaAsLin2 filtering thresholds
# Note: MIN_PREVALENCE can be set per-file using the pattern below
DEFAULT_MIN_PREVALENCE <- 0.0
SCFA_MIN_PREVALENCE    <- 0.05   # stricter threshold for SCFA pathways

MIN_ABUNDANCE <- 0.0

# =============================================================================
# HELPER — convert HUMAnN sample IDs to metadata-compatible format
# Adjust the regex if your sample ID format is different.
# Default pattern: "sub.110.0_initial_..." -> "SUB110-0"
# =============================================================================

convert_sample_ids <- function(ids) {
  gsub("sub\\.([0-9]+)\\.([0-9]+)_.*", "SUB\\1-\\2", ids, perl = TRUE)
}

# =============================================================================
# MAIN LOOP — iterate over all metadata groups × feature files
# =============================================================================

for (meta_name in names(METADATA_FILES)) {

  metadata_path <- METADATA_FILES[[meta_name]]
  metadata <- read.table(metadata_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)

  metadata[[TARGET_COL]] <- factor(metadata[[TARGET_COL]])
  metadata[[TARGET_COL]] <- relevel(metadata[[TARGET_COL]], ref = REFERENCE_LEVEL)
  rownames(metadata) <- metadata[[SAMPLE_ID_COL]]
  metadata[[SAMPLE_ID_COL]] <- NULL

  for (feature_path in FEATURE_FILES) {

    message("\nRunning: ", basename(feature_path), " | Metadata: ", meta_name)

    # ---- Load and transpose features ----
    features <- read.csv(feature_path, header = TRUE, stringsAsFactors = FALSE,
                         check.names = FALSE)

    features_t <- as.data.frame(t(features[, -1]))
    colnames(features_t) <- features[[1]]
    features_t[[SAMPLE_ID_COL]] <- convert_sample_ids(rownames(features_t))

    # ---- Align to metadata ----
    features_t <- features_t[features_t[[SAMPLE_ID_COL]] %in% rownames(metadata), ]
    metadata_sub <- metadata[rownames(metadata) %in% features_t[[SAMPLE_ID_COL]], , drop = FALSE]

    features_t <- features_t[order(features_t[[SAMPLE_ID_COL]]), ]
    metadata_sub <- metadata_sub[order(rownames(metadata_sub)), , drop = FALSE]

    if (!all(features_t[[SAMPLE_ID_COL]] == rownames(metadata_sub))) {
      stop("Sample ID mismatch for: ", basename(feature_path), " / ", meta_name)
    }

    rownames(features_t) <- features_t[[SAMPLE_ID_COL]]
    features_t[[SAMPLE_ID_COL]] <- NULL

    cat("  Samples:", nrow(features_t), "| Features:", ncol(features_t), "\n")

    # ---- Set min_prevalence based on file type ----
    file_label  <- tools::file_path_sans_ext(basename(feature_path))
    min_prev    <- if (grepl("scfa", file_label, ignore.case = TRUE)) SCFA_MIN_PREVALENCE else DEFAULT_MIN_PREVALENCE

    # ---- Output directory ----
    output_dir <- file.path(OUTPUT_ROOT, meta_name, file_label)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # ---- Run MaAsLin2 ----
    Maaslin2(
      input_data     = features_t,
      input_metadata = metadata_sub,
      output         = output_dir,
      fixed_effects  = FIXED_EFFECTS,
      reference      = paste0(TARGET_COL, ",", REFERENCE_LEVEL),
      transform      = "AST",
      min_abundance  = MIN_ABUNDANCE,
      min_prevalence = min_prev
    )

    cat("  Done. Results in:", output_dir, "\n")
  }
}

cat("\nAll MaAsLin2 runs completed.\n")
