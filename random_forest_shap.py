"""
Random Forest Classification + SHAP Feature Importance
=======================================================
Trains a Random Forest classifier on microbiome abundance features to predict
a binary clinical outcome. Uses 5-fold cross-validation with infant-level
grouping (all samples from the same infant stay in the same fold), iterative
Gini-based feature selection, and SHAP to explain model predictions.

Outputs per age condition:
  - Confusion matrix (counts + normalized)
  - Per-fold and average AUC / accuracy / F1
  - SHAP summary plot (PDF)
  - SHAP values table (TSV)

Usage:
    python random_forest_shap.py
"""

import re
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shap
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import (accuracy_score, roc_auc_score,
                             precision_recall_fscore_support,
                             confusion_matrix)
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_class_weight

# =============================================================================
# CONFIG — set BASE_DIR to your project root, then adjust filenames if needed
# =============================================================================

BASE_DIR =    # <-- set this to your project directory

FEATURES_FILE   = os.path.join(BASE_DIR, "features.tsv")       # samples as rows, taxa/features as columns
METADATA_FILE   = os.path.join(BASE_DIR, "metadata.tsv")        # must contain SAMPLE_ID_COL, TARGET_COL, AGE_COL, INFANT_ID_COL
OUTPUT_DIR      = os.path.join(BASE_DIR, "output", "random_forest")  # directory for all outputs

# Column names in metadata
SAMPLE_ID_COL   = "sampleID"
TARGET_COL      = "symptoms"
CONTROL_LABEL   = "Control"       # label for the negative/reference class
CASE_LABEL      = None            # set to a string to restrict positives; None = everything that is not CONTROL_LABEL
AGE_COL         = "visit_age_mo"
INFANT_ID_COL   = "InfantID"      # used for grouped CV; set to None to derive from sample ID

# Age-based subgroups to run (set to None to run on all samples only)
# Each entry: (label, age_filter_function)
AGE_CONDITIONS = [
    ("All_samples",   lambda df: df),
    ("Age_le_4mo",    lambda df: df[df[AGE_COL] <= 4]),
    ("Age_gt_6mo",    lambda df: df[df[AGE_COL] > 6]),
]

# Feature filtering
MEAN_ABUNDANCE_THRESHOLD = 0.07   # percent; features below this are dropped
MAX_FEATURES_AFTER_SELECTION = 120  # iterative Gini selection stops here

# Random Forest hyperparameters (best params from hyperparameter search)
N_ESTIMATORS      = 150
MAX_DEPTH         = 20
MIN_SAMPLES_SPLIT = 5
N_SPLITS          = 5
RANDOM_STATE      = 42

# SHAP plot settings
SHAP_MAX_FEATURES = 30       # how many features to show in SHAP plot
SHAP_XLIM         = (-0.2, 0.2)  # fixed x-axis range for comparability across conditions

# =============================================================================
# HELPERS
# =============================================================================

def sanitize(name):
    return re.sub(r"[^\w\-_\.]", "_", name)


def load_data(features_file, metadata_file):
    features = pd.read_csv(features_file, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata_file, sep="\t", index_col=0)

    if "SubjectID" in features.columns:
        features = features.drop(columns=["SubjectID"])

    # Derive InfantID from sample index if column is missing
    if INFANT_ID_COL not in metadata.columns:
        metadata[INFANT_ID_COL] = metadata.index.str.split("-").str[0]

    return features, metadata


def filter_and_normalize(features, metadata):
    features = features.loc[metadata.index]
    features = features.select_dtypes(include=[np.number])

    # Mean abundance filter
    mean_abund = features.mean(axis=0) * 100
    features = features.loc[:, mean_abund > MEAN_ABUNDANCE_THRESHOLD]

    scaler = StandardScaler()
    features_norm = pd.DataFrame(
        scaler.fit_transform(features),
        index=features.index,
        columns=features.columns
    )
    return features_norm


def iterative_feature_selection(X, y):
    """Drop 10 least-important features at a time until MAX_FEATURES_AFTER_SELECTION remain."""
    feature_set = X.copy()
    while feature_set.shape[1] > MAX_FEATURES_AFTER_SELECTION:
        rf_temp = RandomForestClassifier(
            n_estimators=100, random_state=RANDOM_STATE
        )
        rf_temp.fit(feature_set, y)
        importances = pd.Series(rf_temp.feature_importances_, index=feature_set.columns)
        feature_set = feature_set.drop(columns=importances.nsmallest(10).index)
    return feature_set


def encode_target(metadata):
    y = metadata[TARGET_COL]
    if CASE_LABEL is not None:
        # Binary: CASE_LABEL=1, CONTROL_LABEL=0, rest dropped
        mask = y.isin([CONTROL_LABEL, CASE_LABEL])
        y = y[mask]
    y_encoded = y.map(lambda v: 0 if v == CONTROL_LABEL else 1)
    return y_encoded


def grouped_kfold_splits(metadata, y, n_splits):
    """Split by infant so all samples from one infant are in the same fold."""
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_STATE)
    infants = metadata[INFANT_ID_COL].unique()
    for train_inf_idx, test_inf_idx in kf.split(infants):
        train_infants = infants[train_inf_idx]
        test_infants  = infants[test_inf_idx]
        train_idx = metadata[metadata[INFANT_ID_COL].isin(train_infants)].index
        test_idx  = metadata[metadata[INFANT_ID_COL].isin(test_infants)].index
        yield train_idx, test_idx


# =============================================================================
# MAIN
# =============================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

features_all, metadata_all = load_data(FEATURES_FILE, METADATA_FILE)

for condition_label, age_filter in AGE_CONDITIONS:

    print(f"\n{'='*60}")
    print(f"Condition: {condition_label}")
    print(f"{'='*60}")

    # ---- Filter samples ----
    meta_cond = age_filter(metadata_all)
    y_all     = encode_target(meta_cond)
    meta_cond = meta_cond.loc[y_all.index]   # drop excluded labels if any

    # ---- Normalize features ----
    X = filter_and_normalize(features_all, meta_cond)
    y = y_all.loc[X.index]
    meta_cond = meta_cond.loc[X.index]

    print(f"Samples: {len(X)}  |  Features before selection: {X.shape[1]}")

    # ---- Iterative feature selection ----
    X_selected = iterative_feature_selection(X, y)
    print(f"Features after selection: {X_selected.shape[1]}")

    # ---- Cross-validation ----
    fold_metrics = []
    all_shap_values = []
    all_test_indices = []
    all_cm = np.zeros((2, 2), dtype=int)

    for fold, (train_idx, test_idx) in enumerate(
            grouped_kfold_splits(meta_cond, y, N_SPLITS), start=1):

        X_train = X_selected.loc[train_idx]
        X_test  = X_selected.loc[test_idx]
        y_train = y.loc[train_idx]
        y_test  = y.loc[test_idx]

        class_weights = compute_class_weight(
            "balanced", classes=np.unique(y_train), y=y_train
        )
        cw_dict = dict(enumerate(class_weights))

        rf = RandomForestClassifier(
            n_estimators=N_ESTIMATORS,
            max_depth=MAX_DEPTH,
            min_samples_split=MIN_SAMPLES_SPLIT,
            class_weight=cw_dict,
            random_state=RANDOM_STATE
        )
        rf.fit(X_train, y_train)

        y_pred  = rf.predict(X_test)
        y_prob  = rf.predict_proba(X_test)[:, 1]

        auc = roc_auc_score(y_test, y_prob) if len(np.unique(y_test)) > 1 else float("nan")
        acc = accuracy_score(y_test, y_pred)
        prec, rec, f1, _ = precision_recall_fscore_support(
            y_test, y_pred, average="binary", zero_division=0
        )
        fold_metrics.append({"Fold": fold, "AUC": auc, "Accuracy": acc,
                              "Precision": prec, "Recall": rec, "F1": f1})
        all_cm += confusion_matrix(y_test, y_pred, labels=[0, 1])

        # SHAP
        explainer   = shap.TreeExplainer(rf)
        shap_vals   = explainer.shap_values(X_test)
        if isinstance(shap_vals, list):
            shap_vals = shap_vals[1]
        elif shap_vals.ndim == 3:
            shap_vals = shap_vals[:, :, 1]

        all_shap_values.append(shap_vals)
        all_test_indices.append(X_test.index)

    # ---- Metrics summary ----
    metrics_df = pd.DataFrame(fold_metrics)
    print(metrics_df.to_string(index=False))
    print("\nMean:")
    print(metrics_df.drop(columns="Fold").mean().to_string())

    metrics_df.to_csv(
        os.path.join(OUTPUT_DIR, f"metrics_{sanitize(condition_label)}.tsv"),
        sep="\t", index=False
    )

    # ---- Confusion matrices ----
    cm_norm = all_cm.astype(float) / all_cm.sum(axis=1, keepdims=True)
    pd.DataFrame(all_cm,
                 index=[f"True_{CONTROL_LABEL}", "True_Case"],
                 columns=[f"Pred_{CONTROL_LABEL}", "Pred_Case"]
    ).to_csv(os.path.join(OUTPUT_DIR, f"cm_counts_{sanitize(condition_label)}.tsv"), sep="\t")

    pd.DataFrame(cm_norm,
                 index=[f"True_{CONTROL_LABEL}", "True_Case"],
                 columns=[f"Pred_{CONTROL_LABEL}", "Pred_Case"]
    ).to_csv(os.path.join(OUTPUT_DIR, f"cm_norm_{sanitize(condition_label)}.tsv"), sep="\t")

    # ---- SHAP summary ----
    combined_indices = np.concatenate([idx for idx in all_test_indices])
    # Map string indices to positional for iloc
    pos_map = {idx: i for i, idx in enumerate(X_selected.index)}
    combined_pos = [pos_map[i] for i in combined_indices]

    combined_shap = np.concatenate(all_shap_values, axis=0)
    X_test_combined = X_selected.iloc[combined_pos]

    # Save SHAP values
    shap_df = pd.DataFrame(combined_shap, columns=X_selected.columns,
                           index=X_test_combined.index)
    shap_df.to_csv(
        os.path.join(OUTPUT_DIR, f"shap_values_{sanitize(condition_label)}.tsv"),
        sep="\t"
    )

    # SHAP plot — top features by mean absolute SHAP
    mean_abs_shap = np.abs(combined_shap).mean(axis=0)
    top_idx = np.argsort(mean_abs_shap)[::-1][:SHAP_MAX_FEATURES]
    top_features = X_selected.columns[top_idx]

    shap_plot_vals = combined_shap[:, top_idx]
    X_plot = X_test_combined[top_features]

    # Clean up feature labels (strip "s__" prefix if MetaPhlAn format)
    clean_labels = [f.split("s__")[-1] if "s__" in f else f for f in top_features]
    X_plot_labeled = X_plot.copy()
    X_plot_labeled.columns = clean_labels

    plt.figure(figsize=(10, 8))
    shap.summary_plot(
        shap_plot_vals, X_plot_labeled,
        max_display=SHAP_MAX_FEATURES,
        show=False, plot_type="dot",
        plot_size=(10, 8), color_bar=True
    )
    ax = plt.gca()
    ax.set_xlim(SHAP_XLIM)
    plt.title(f"SHAP — {condition_label}", fontsize=14)
    plt.tight_layout()

    shap_plot_path = os.path.join(
        OUTPUT_DIR, f"shap_plot_{sanitize(condition_label)}.pdf"
    )
    plt.savefig(shap_plot_path, bbox_inches="tight")
    plt.close()
    print(f"SHAP plot saved: {shap_plot_path}")

print("\nAll conditions done. Outputs in:", OUTPUT_DIR)
