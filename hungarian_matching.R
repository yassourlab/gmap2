# =============================================================================
# Hungarian Algorithm — Age-matched case-control pairing
#
# Matches each case sample (e.g. Symptomatic, Pre-symptoms, Resolved) to a
# control sample of the closest visit age, using the Hungarian algorithm
# (scipy LSAP / clue::solve_LSAP). Controls are assigned without replacement
# across groups: a control used for one group is excluded from the next.
#
# For each case group, one representative sample per infant is selected:
#   - Pre-symptoms: the LAST sample before symptom onset
#   - Symptomatic:  the FIRST symptomatic sample
#   - Resolved:     the FIRST resolved sample
#
# Outputs:
#   - One matched TSV per group (cases + matched controls)
#   - One scatter plot per group showing age alignment (PDF)
#   - One combined plot with all groups (PDF)
#
# Usage: Rscript hungarian_matching.R
# =============================================================================

library(clue)     # solve_LSAP — Hungarian algorithm
library(dplyr)
library(readr)
library(ggplot2)

# =============================================================================
# CONFIG — set these paths and parameters before running
# =============================================================================

BASE_DIR      <-    # <-- fill it with path to the dir to the samples

METADATA_FILE <- file.path(BASE_DIR, "metadata.tsv")
OUTPUT_DIR    <- file.path(BASE_DIR, "output", "hungarian_matching")

# Column names in metadata
SAMPLE_ID_COL    <- "sampleID"
SUBJECT_ID_COL   <- "SubjectID"    # per-infant identifier; derived from sampleID if absent
AGE_COL          <- "visit_age_mo"
OUTCOME_COL      <- "symptoms"

# Group labels in OUTCOME_COL
CONTROL_LABEL       <- "Control"
CASE_GROUPS <- list(
  PreSymptoms = list(label = "Pre-symptoms", select = "last"),   # last sample before onset
  Symptomatic = list(label = "Symptomatic",  select = "first"),  # first symptomatic sample
  Resolved    = list(label = "Resolved",     select = "first")   # first resolved sample
)

# Columns to keep in output matched tables
OUTPUT_COLS <- c("sampleID", "visit_age_mo", "symptoms", "mode_of_delivery", "probiotics_firstyr")

# Plot appearance
CASE_COLORS <- c(
  "Control"      = "#4393C3",
  "Pre-symptoms" = "#F4A582",
  "Symptomatic"  = "#D6604D",
  "Resolved"     = "#74C476"
)
AGE_AXIS_LIMITS <- c(-1, 14)   # months; adjust to your cohort's age range

# =============================================================================
# LOAD & PREPARE METADATA
# =============================================================================

metadata <- read_tsv(METADATA_FILE, show_col_types = FALSE)

required_cols <- c(SAMPLE_ID_COL, AGE_COL, OUTCOME_COL)
missing <- setdiff(required_cols, colnames(metadata))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

metadata[[SAMPLE_ID_COL]] <- toupper(metadata[[SAMPLE_ID_COL]])

# Derive SubjectID from sampleID if not present (expects format SUB###-visit)
if (!SUBJECT_ID_COL %in% colnames(metadata)) {
  metadata[[SUBJECT_ID_COL]] <- gsub("-.*", "", metadata[[SAMPLE_ID_COL]])
}

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# =============================================================================
# HELPERS
# =============================================================================

# Select one representative sample per infant for a given symptom group.
# select = "first" takes the earliest visit; "last" takes the latest.
select_representative <- function(df, symptom_label, select = "first") {
  df %>%
    filter(.data[[OUTCOME_COL]] == symptom_label) %>%
    arrange(.data[[AGE_COL]]) %>%
    group_by(.data[[SUBJECT_ID_COL]]) %>%
    slice(if (select == "first") 1L else n()) %>%
    ungroup()
}

# Match cases to controls by minimum age difference (Hungarian algorithm).
# Returns a combined data frame of matched cases + their paired controls.
match_by_age <- function(cases, controls) {
  if (nrow(cases) == 0 || nrow(controls) == 0) return(NULL)

  cost_matrix <- abs(outer(cases[[AGE_COL]], controls[[AGE_COL]], "-"))
  assignment  <- solve_LSAP(cost_matrix, maximum = FALSE)

  matched_controls <- controls[as.integer(assignment), ]
  bind_rows(cases, matched_controls)
}

# =============================================================================
# EXTRACT REPRESENTATIVE SAMPLES
# =============================================================================

case_samples <- lapply(CASE_GROUPS, function(grp) {
  select_representative(metadata, grp$label, grp$select)
})
names(case_samples) <- names(CASE_GROUPS)

# =============================================================================
# HUNGARIAN MATCHING (sequential — controls are not reused across groups)
# =============================================================================

available_controls <- metadata %>% filter(.data[[OUTCOME_COL]] == CONTROL_LABEL)
matched <- list()

for (grp_name in names(CASE_GROUPS)) {
  cases <- case_samples[[grp_name]]
  result <- match_by_age(cases, available_controls)

  if (!is.null(result)) {
    # Remove used controls from the pool
    used_control_ids <- result %>%
      filter(.data[[OUTCOME_COL]] == CONTROL_LABEL) %>%
      pull(.data[[SAMPLE_ID_COL]])
    available_controls <- available_controls %>%
      filter(!.data[[SAMPLE_ID_COL]] %in% used_control_ids)
  }

  matched[[grp_name]] <- result
  cat(grp_name, ": matched", nrow(cases), "cases to",
      sum(result[[OUTCOME_COL]] == CONTROL_LABEL), "controls\n")
}

# =============================================================================
# SAVE MATCHED TABLES
# =============================================================================

out_cols <- intersect(OUTPUT_COLS, colnames(metadata))

for (grp_name in names(matched)) {
  if (is.null(matched[[grp_name]])) next
  label <- CASE_GROUPS[[grp_name]]$label
  out_file <- file.path(OUTPUT_DIR, paste0("matched_", tolower(gsub("[-/ ]", "_", label)), ".tsv"))
  matched[[grp_name]] %>%
    select(all_of(out_cols)) %>%
    write_tsv(out_file)
  cat("Saved:", out_file, "\n")
}

# =============================================================================
# VISUALIZE — age alignment scatter plots
# =============================================================================

make_scatter <- function(data, title) {
  ggplot(data, aes(x = .data[[AGE_COL]], y = .data[[AGE_COL]],
                   fill = .data[[OUTCOME_COL]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(shape = 21, size = 3, alpha = 0.9, color = "black", stroke = 1.0) +
    scale_fill_manual(values = CASE_COLORS, name = "Group") +
    coord_cartesian(xlim = AGE_AXIS_LIMITS, ylim = AGE_AXIS_LIMITS) +
    labs(title = paste("Hungarian Matching —", title),
         x = "Visit age (months)", y = "Visit age (months)") +
    theme_minimal(base_size = 13) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          axis.line = element_line(color = "black"),
          legend.key = element_rect(fill = "white", color = NA))
}

# Individual plots
for (grp_name in names(matched)) {
  if (is.null(matched[[grp_name]])) next
  label <- CASE_GROUPS[[grp_name]]$label
  p <- make_scatter(matched[[grp_name]], label)
  out_pdf <- file.path(OUTPUT_DIR, paste0("matching_plot_", tolower(gsub("[-/ ]", "_", label)), ".pdf"))
  ggsave(out_pdf, plot = p, width = 6, height = 5, dpi = 300, bg = "white")
  cat("Plot saved:", out_pdf, "\n")
}

# Combined plot (all groups together)
all_matched <- bind_rows(matched)
if (nrow(all_matched) > 0) {
  p_all <- make_scatter(all_matched, "All Groups")
  ggsave(file.path(OUTPUT_DIR, "matching_plot_all.pdf"),
         plot = p_all, width = 7, height = 6, dpi = 300, bg = "white")
  cat("Combined plot saved.\n")
}

cat("\nDone. All outputs in:", OUTPUT_DIR, "\n")
