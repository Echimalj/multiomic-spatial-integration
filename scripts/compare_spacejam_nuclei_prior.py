"""
Compare baseline and nuclei-informed SpaceJam posterior abundances.

This script compares:

    results/spacejam/
    results/spacejam_nuclei_informed/

for ADCAA and CTRL independently and jointly.

Outputs are written to:

    results/spacejam_nuclei_prior_comparison/

The script does not modify either SpaceJam result directory.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable

import matplotlib

# Use a non-interactive backend for HPC execution.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, spearmanr


COHORTS = ("ADCAA", "CTRL")

ROI_COLUMN_CANDIDATES = (
    "ROI_ID",
    "roi_id",
    "ROI",
    "roi",
    "segment",
    "Segment",
    "segment_id",
    "SegmentID",
    "roi_label",
    "ROI_label",
)

NUCLEI_COLUMN_CANDIDATES = (
    "nuclei_count",
    "Nuclei_Count",
    "nuc_count",
    "NucCount",
    "nuclei",
    "Nuclei",
    "nucleus_count",
    "Nucleus_Count",
    "count",
)

MANIFEST_VALUE_CANDIDATES = {
    "roi": (
        "ROI_ID",
        "roi_id",
        "ROI",
        "roi",
        "spot",
        "spot_id",
        "index",
        "label",
        "name",
    ),
    "factor": (
        "factor",
        "Factor",
        "factor_name",
        "celltype",
        "cell_type",
        "cell type",
        "label",
        "name",
        "index",
    ),
}


# ============================================================
# General helpers
# ============================================================


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare baseline and nuclei-informed SpaceJam "
            "posterior abundance estimates."
        )
    )

    parser.add_argument(
        "--project-dir",
        type=Path,
        default=Path(
            "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
            "Human_GeoMx_Sep2025/multiomic-spatial-integration"
        ),
        help="Repository root.",
    )

    parser.add_argument(
        "--baseline-dir",
        type=Path,
        default=None,
        help=(
            "Baseline SpaceJam result directory. Defaults to "
            "<project-dir>/results/spacejam."
        ),
    )

    parser.add_argument(
        "--nuclei-informed-dir",
        type=Path,
        default=None,
        help=(
            "Nuclei-informed result directory. Defaults to "
            "<project-dir>/results/spacejam_nuclei_informed."
        ),
    )

    parser.add_argument(
        "--nuclei-file",
        type=Path,
        default=None,
        help=(
            "ROI nuclei-count CSV. Defaults to "
            "<project-dir>/data/nuc_count.csv."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Comparison output directory. Defaults to "
            "<project-dir>/results/"
            "spacejam_nuclei_prior_comparison."
        ),
    )

    parser.add_argument(
        "--heatmap-top-factors",
        type=int,
        default=20,
        help=(
            "Number of factors with the largest mean absolute "
            "relative-abundance change to include in heatmaps."
        ),
    )

    parser.add_argument(
        "--heatmap-top-rois",
        type=int,
        default=100,
        help=(
            "Maximum number of ROIs with the largest overall "
            "relative-abundance change to include in heatmaps."
        ),
    )

    return parser.parse_args()


def ensure_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file does not exist: {path}")


def ensure_directory(path: Path) -> None:
    if not path.is_dir():
        raise FileNotFoundError(
            f"Required directory does not exist: {path}"
        )


def clean_label(value: Any) -> str:
    """Convert a manifest identifier into a clean string."""

    if pd.isna(value):
        raise ValueError("Manifest contains a missing identifier.")

    label = str(value).strip()

    # Remove the common artifact produced when integer identifiers
    # are read as floating-point values.
    if re.fullmatch(r"-?\d+\.0", label):
        label = label[:-2]

    if not label:
        raise ValueError("Manifest contains an empty identifier.")

    return label


def select_manifest_column(
    frame: pd.DataFrame,
    manifest_type: str,
) -> str:
    """
    Identify the column containing ROI or factor labels.

    The exported manifests may contain an unnamed index column plus a
    named label column. Named candidate columns are preferred.
    """

    candidates = MANIFEST_VALUE_CANDIDATES[manifest_type]

    for candidate in candidates:
        if candidate in frame.columns:
            return candidate

    non_unnamed = [
        col
        for col in frame.columns
        if not str(col).startswith("Unnamed:")
    ]

    if len(non_unnamed) == 1:
        return non_unnamed[0]

    if len(frame.columns) == 1:
        return frame.columns[0]

    object_like = [
        col
        for col in non_unnamed
        if (
            pd.api.types.is_object_dtype(frame[col])
            or pd.api.types.is_string_dtype(frame[col])
        )
    ]

    if len(object_like) == 1:
        return object_like[0]

    raise ValueError(
        f"Could not identify the {manifest_type} label column. "
        f"Columns were: {frame.columns.tolist()}"
    )


def read_manifest(
    path: Path,
    manifest_type: str,
) -> list[str]:
    ensure_file(path)

    frame = pd.read_csv(path)
    column = select_manifest_column(frame, manifest_type)

    values = [clean_label(value) for value in frame[column]]

    duplicated = pd.Series(values).duplicated(keep=False)

    if duplicated.any():
        duplicate_values = sorted(
            pd.Series(values)[duplicated].unique().tolist()
        )
        raise ValueError(
            f"Duplicate labels in {path}: {duplicate_values[:10]}"
        )

    return values


def load_tensor(path: Path) -> np.ndarray:
    ensure_file(path)

    loaded = torch.load(
        path,
        map_location="cpu",
        weights_only=False,
    )

    array = (
        torch.as_tensor(loaded)
        .detach()
        .cpu()
        .numpy()
        .astype(np.float64, copy=False)
    )

    if array.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional tensor in {path}, "
            f"found shape {array.shape}."
        )

    if not np.isfinite(array).all():
        raise ValueError(
            f"Non-finite values detected in tensor: {path}"
        )

    return array


def safe_pearson(
    x: Iterable[float],
    y: Iterable[float],
) -> tuple[float, float]:
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)

    valid = np.isfinite(x_array) & np.isfinite(y_array)
    x_array = x_array[valid]
    y_array = y_array[valid]

    if len(x_array) < 3:
        return math.nan, math.nan

    if np.std(x_array) == 0 or np.std(y_array) == 0:
        return math.nan, math.nan

    result = pearsonr(x_array, y_array)

    return float(result.statistic), float(result.pvalue)


def safe_spearman(
    x: Iterable[float],
    y: Iterable[float],
) -> tuple[float, float]:
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)

    valid = np.isfinite(x_array) & np.isfinite(y_array)
    x_array = x_array[valid]
    y_array = y_array[valid]

    if len(x_array) < 3:
        return math.nan, math.nan

    if len(np.unique(x_array)) < 2:
        return math.nan, math.nan

    if len(np.unique(y_array)) < 2:
        return math.nan, math.nan

    result = spearmanr(x_array, y_array)

    return float(result.statistic), float(result.pvalue)


def summarize_numeric(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]

    if len(array) == 0:
        return {
            "n": 0,
            "min": math.nan,
            "q25": math.nan,
            "median": math.nan,
            "mean": math.nan,
            "q75": math.nan,
            "max": math.nan,
        }

    return {
        "n": int(len(array)),
        "min": float(np.min(array)),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "mean": float(np.mean(array)),
        "q75": float(np.quantile(array, 0.75)),
        "max": float(np.max(array)),
    }


def json_safe(value: Any) -> Any:
    """Convert NumPy and non-finite values into JSON-safe objects."""

    if isinstance(value, dict):
        return {
            str(key): json_safe(item)
            for key, item in value.items()
        }

    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]

    if isinstance(value, np.generic):
        value = value.item()

    if isinstance(value, float) and not np.isfinite(value):
        return None

    return value


# ============================================================
# Data loading
# ============================================================


def load_cohort(
    directory: Path,
    cohort: str,
) -> dict[str, Any]:
    rois = read_manifest(
        directory / f"{cohort}_manifest_rois.csv",
        manifest_type="roi",
    )

    factors = read_manifest(
        directory / f"{cohort}_manifest_factors.csv",
        manifest_type="factor",
    )

    absolute = load_tensor(
        directory / f"{cohort}_spot_factors_abs.pt"
    )

    relative = load_tensor(
        directory / f"{cohort}_spot_factors_rel.pt"
    )

    expected_shape = (len(rois), len(factors))

    if absolute.shape != expected_shape:
        raise ValueError(
            f"{cohort} absolute tensor shape {absolute.shape} does "
            f"not match manifest shape {expected_shape} in "
            f"{directory}."
        )

    if relative.shape != expected_shape:
        raise ValueError(
            f"{cohort} relative tensor shape {relative.shape} does "
            f"not match manifest shape {expected_shape} in "
            f"{directory}."
        )

    relative_sums = relative.sum(axis=1)

    if not np.allclose(relative_sums, 1.0, atol=1e-4):
        raise ValueError(
            f"{cohort} relative rows do not sum to one in "
            f"{directory}. Range: "
            f"{relative_sums.min()} to {relative_sums.max()}."
        )

    return {
        "rois": rois,
        "factors": factors,
        "absolute": absolute,
        "relative": relative,
    }


def compare_manifests(
    baseline: dict[str, Any],
    informed: dict[str, Any],
    cohort: str,
) -> None:
    if baseline["rois"] != informed["rois"]:
        raise ValueError(
            f"{cohort} ROI manifests differ between baseline and "
            "nuclei-informed outputs."
        )

    if baseline["factors"] != informed["factors"]:
        raise ValueError(
            f"{cohort} factor manifests differ between baseline and "
            "nuclei-informed outputs."
        )


def find_column(
    frame: pd.DataFrame,
    candidates: tuple[str, ...],
    description: str,
) -> str:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate

    lower_map = {
        str(column).strip().lower(): column
        for column in frame.columns
    }

    for candidate in candidates:
        match = lower_map.get(candidate.lower())
        if match is not None:
            return match

    raise ValueError(
        f"Could not identify the {description} column. "
        f"Available columns: {frame.columns.tolist()}"
    )


def read_nuclei_counts(path: Path) -> pd.DataFrame:
    ensure_file(path)

    frame = pd.read_csv(path)

    roi_column = find_column(
        frame,
        ROI_COLUMN_CANDIDATES,
        "ROI identifier",
    )

    nuclei_column = find_column(
        frame,
        NUCLEI_COLUMN_CANDIDATES,
        "nuclei-count",
    )

    output = frame[[roi_column, nuclei_column]].copy()
    output.columns = ["roi", "nuclei_count"]

    output["roi"] = output["roi"].map(clean_label)
    output["nuclei_count"] = pd.to_numeric(
        output["nuclei_count"],
        errors="coerce",
    )

    if output["roi"].duplicated().any():
        duplicates = sorted(
            output.loc[
                output["roi"].duplicated(keep=False),
                "roi",
            ].unique()
        )
        raise ValueError(
            "Duplicate ROI identifiers in nuclei-count file: "
            f"{duplicates[:10]}"
        )

    if output["nuclei_count"].isna().any():
        bad_rois = output.loc[
            output["nuclei_count"].isna(),
            "roi",
        ].tolist()

        raise ValueError(
            "Missing or non-numeric nuclei counts for ROI(s): "
            f"{bad_rois[:10]}"
        )

    if (output["nuclei_count"] <= 0).any():
        bad_rois = output.loc[
            output["nuclei_count"] <= 0,
            "roi",
        ].tolist()

        raise ValueError(
            "Nuclei counts must be positive. Invalid ROI(s): "
            f"{bad_rois[:10]}"
        )

    return output


# ============================================================
# Comparison tables
# ============================================================


def build_roi_metrics(
    cohort: str,
    rois: list[str],
    baseline_absolute: np.ndarray,
    informed_absolute: np.ndarray,
    baseline_relative: np.ndarray,
    informed_relative: np.ndarray,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []

    absolute_difference = (
        informed_absolute - baseline_absolute
    )
    relative_difference = (
        informed_relative - baseline_relative
    )

    for index, roi in enumerate(rois):
        absolute_pearson, absolute_pearson_p = safe_pearson(
            baseline_absolute[index, :],
            informed_absolute[index, :],
        )
        absolute_spearman, absolute_spearman_p = safe_spearman(
            baseline_absolute[index, :],
            informed_absolute[index, :],
        )

        relative_pearson, relative_pearson_p = safe_pearson(
            baseline_relative[index, :],
            informed_relative[index, :],
        )
        relative_spearman, relative_spearman_p = safe_spearman(
            baseline_relative[index, :],
            informed_relative[index, :],
        )

        baseline_total = float(
            baseline_absolute[index, :].sum()
        )
        informed_total = float(
            informed_absolute[index, :].sum()
        )

        total_difference = informed_total - baseline_total
        total_percent_change = (
            100.0 * total_difference / baseline_total
            if baseline_total > 0
            else math.nan
        )

        records.append(
            {
                "cohort": cohort,
                "roi": roi,
                "absolute_pearson_r": absolute_pearson,
                "absolute_pearson_p": absolute_pearson_p,
                "absolute_spearman_rho": absolute_spearman,
                "absolute_spearman_p": absolute_spearman_p,
                "relative_pearson_r": relative_pearson,
                "relative_pearson_p": relative_pearson_p,
                "relative_spearman_rho": relative_spearman,
                "relative_spearman_p": relative_spearman_p,
                "baseline_total_absolute": baseline_total,
                "informed_total_absolute": informed_total,
                "total_absolute_difference": total_difference,
                "total_absolute_percent_change": (
                    total_percent_change
                ),
                "absolute_mean_absolute_difference": float(
                    np.mean(
                        np.abs(absolute_difference[index, :])
                    )
                ),
                "absolute_root_mean_squared_difference": float(
                    np.sqrt(
                        np.mean(
                            absolute_difference[index, :] ** 2
                        )
                    )
                ),
                "relative_mean_absolute_difference": float(
                    np.mean(
                        np.abs(relative_difference[index, :])
                    )
                ),
                "relative_root_mean_squared_difference": float(
                    np.sqrt(
                        np.mean(
                            relative_difference[index, :] ** 2
                        )
                    )
                ),
                "relative_max_absolute_difference": float(
                    np.max(
                        np.abs(relative_difference[index, :])
                    )
                ),
                "relative_l1_distance": float(
                    np.sum(
                        np.abs(relative_difference[index, :])
                    )
                ),
            }
        )

    return pd.DataFrame.from_records(records)


def build_factor_metrics(
    cohort: str,
    factors: list[str],
    baseline_absolute: np.ndarray,
    informed_absolute: np.ndarray,
    baseline_relative: np.ndarray,
    informed_relative: np.ndarray,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []

    for index, factor in enumerate(factors):
        baseline_abs = baseline_absolute[:, index]
        informed_abs = informed_absolute[:, index]

        baseline_rel = baseline_relative[:, index]
        informed_rel = informed_relative[:, index]

        absolute_difference = informed_abs - baseline_abs
        relative_difference = informed_rel - baseline_rel

        absolute_pearson, absolute_pearson_p = safe_pearson(
            baseline_abs,
            informed_abs,
        )
        absolute_spearman, absolute_spearman_p = safe_spearman(
            baseline_abs,
            informed_abs,
        )

        relative_pearson, relative_pearson_p = safe_pearson(
            baseline_rel,
            informed_rel,
        )
        relative_spearman, relative_spearman_p = safe_spearman(
            baseline_rel,
            informed_rel,
        )

        baseline_abs_mean = float(np.mean(baseline_abs))
        informed_abs_mean = float(np.mean(informed_abs))

        baseline_rel_mean = float(np.mean(baseline_rel))
        informed_rel_mean = float(np.mean(informed_rel))

        absolute_mean_percent_change = (
            100.0
            * (informed_abs_mean - baseline_abs_mean)
            / baseline_abs_mean
            if baseline_abs_mean > 0
            else math.nan
        )

        relative_mean_percent_change = (
            100.0
            * (informed_rel_mean - baseline_rel_mean)
            / baseline_rel_mean
            if baseline_rel_mean > 0
            else math.nan
        )

        records.append(
            {
                "cohort": cohort,
                "factor": factor,
                "n_rois": int(len(baseline_abs)),
                "absolute_pearson_r": absolute_pearson,
                "absolute_pearson_p": absolute_pearson_p,
                "absolute_spearman_rho": absolute_spearman,
                "absolute_spearman_p": absolute_spearman_p,
                "relative_pearson_r": relative_pearson,
                "relative_pearson_p": relative_pearson_p,
                "relative_spearman_rho": relative_spearman,
                "relative_spearman_p": relative_spearman_p,
                "baseline_absolute_mean": baseline_abs_mean,
                "informed_absolute_mean": informed_abs_mean,
                "absolute_mean_difference": float(
                    np.mean(absolute_difference)
                ),
                "absolute_mean_absolute_difference": float(
                    np.mean(np.abs(absolute_difference))
                ),
                "absolute_median_difference": float(
                    np.median(absolute_difference)
                ),
                "absolute_mean_percent_change": (
                    absolute_mean_percent_change
                ),
                "baseline_relative_mean": baseline_rel_mean,
                "informed_relative_mean": informed_rel_mean,
                "relative_mean_difference": float(
                    np.mean(relative_difference)
                ),
                "relative_mean_absolute_difference": float(
                    np.mean(np.abs(relative_difference))
                ),
                "relative_median_difference": float(
                    np.median(relative_difference)
                ),
                "relative_mean_percent_change": (
                    relative_mean_percent_change
                ),
                "fraction_rois_relative_increased": float(
                    np.mean(relative_difference > 0)
                ),
                "fraction_rois_relative_decreased": float(
                    np.mean(relative_difference < 0)
                ),
                "relative_difference_q05": float(
                    np.quantile(relative_difference, 0.05)
                ),
                "relative_difference_q95": float(
                    np.quantile(relative_difference, 0.95)
                ),
            }
        )

    return pd.DataFrame.from_records(records)


def build_long_abundance_table(
    cohort: str,
    rois: list[str],
    factors: list[str],
    baseline_absolute: np.ndarray,
    informed_absolute: np.ndarray,
    baseline_relative: np.ndarray,
    informed_relative: np.ndarray,
) -> pd.DataFrame:
    roi_values = np.repeat(rois, len(factors))
    factor_values = np.tile(factors, len(rois))

    baseline_abs_flat = baseline_absolute.reshape(-1)
    informed_abs_flat = informed_absolute.reshape(-1)

    baseline_rel_flat = baseline_relative.reshape(-1)
    informed_rel_flat = informed_relative.reshape(-1)

    return pd.DataFrame(
        {
            "cohort": cohort,
            "roi": roi_values,
            "factor": factor_values,
            "baseline_absolute": baseline_abs_flat,
            "informed_absolute": informed_abs_flat,
            "absolute_difference": (
                informed_abs_flat - baseline_abs_flat
            ),
            "baseline_relative": baseline_rel_flat,
            "informed_relative": informed_rel_flat,
            "relative_difference": (
                informed_rel_flat - baseline_rel_flat
            ),
        }
    )


def add_nuclei_counts(
    roi_metrics: pd.DataFrame,
    nuclei_counts: pd.DataFrame,
) -> pd.DataFrame:
    merged = roi_metrics.merge(
        nuclei_counts,
        on="roi",
        how="left",
        validate="one_to_one",
    )

    missing = merged.loc[
        merged["nuclei_count"].isna(),
        "roi",
    ].tolist()

    if missing:
        raise ValueError(
            "The following SpaceJam ROIs were not found in the "
            f"nuclei-count file: {missing[:20]}"
        )

    merged["baseline_total_per_nucleus"] = (
        merged["baseline_total_absolute"]
        / merged["nuclei_count"]
    )

    merged["informed_total_per_nucleus"] = (
        merged["informed_total_absolute"]
        / merged["nuclei_count"]
    )

    merged["baseline_minus_nuclei"] = (
        merged["baseline_total_absolute"]
        - merged["nuclei_count"]
    )

    merged["informed_minus_nuclei"] = (
        merged["informed_total_absolute"]
        - merged["nuclei_count"]
    )

    merged["baseline_absolute_error_vs_nuclei"] = np.abs(
        merged["baseline_minus_nuclei"]
    )

    merged["informed_absolute_error_vs_nuclei"] = np.abs(
        merged["informed_minus_nuclei"]
    )

    merged["absolute_error_improvement_vs_nuclei"] = (
        merged["baseline_absolute_error_vs_nuclei"]
        - merged["informed_absolute_error_vs_nuclei"]
    )

    return merged


# ============================================================
# Plotting
# ============================================================


def save_figure(
    figure: plt.Figure,
    output_base: Path,
) -> None:
    figure.savefig(
        output_base.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
    )

    figure.savefig(
        output_base.with_suffix(".pdf"),
        bbox_inches="tight",
    )

    plt.close(figure)


def add_identity_line(
    axis: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
) -> None:
    finite = np.isfinite(x) & np.isfinite(y)

    if not finite.any():
        return

    lower = float(min(np.min(x[finite]), np.min(y[finite])))
    upper = float(max(np.max(x[finite]), np.max(y[finite])))

    axis.plot(
        [lower, upper],
        [lower, upper],
        linestyle="--",
        linewidth=1,
    )


def plot_total_baseline_vs_informed(
    frame: pd.DataFrame,
    output_dir: Path,
) -> None:
    figure, axis = plt.subplots(figsize=(7, 6))

    for cohort in COHORTS:
        subset = frame.loc[frame["cohort"] == cohort]

        axis.scatter(
            subset["baseline_total_absolute"],
            subset["informed_total_absolute"],
            alpha=0.75,
            label=cohort,
        )

    x = frame["baseline_total_absolute"].to_numpy()
    y = frame["informed_total_absolute"].to_numpy()

    add_identity_line(axis, x, y)

    pearson_r, _ = safe_pearson(x, y)
    spearman_rho, _ = safe_spearman(x, y)

    axis.set_xlabel("Baseline total inferred abundance")
    axis.set_ylabel("Nuclei-informed total inferred abundance")
    axis.set_title(
        "Total inferred abundance: baseline vs nuclei-informed\n"
        f"Pearson r = {pearson_r:.3f}; "
        f"Spearman ρ = {spearman_rho:.3f}"
    )
    axis.legend(frameon=False)

    save_figure(
        figure,
        output_dir / "total_abundance_baseline_vs_informed",
    )


def plot_totals_vs_nuclei(
    frame: pd.DataFrame,
    output_dir: Path,
) -> None:
    for method, column, label in (
        (
            "baseline",
            "baseline_total_absolute",
            "Baseline",
        ),
        (
            "nuclei_informed",
            "informed_total_absolute",
            "Nuclei-informed",
        ),
    ):
        figure, axis = plt.subplots(figsize=(7, 6))

        for cohort in COHORTS:
            subset = frame.loc[frame["cohort"] == cohort]

            axis.scatter(
                subset["nuclei_count"],
                subset[column],
                alpha=0.75,
                label=cohort,
            )

        x = frame["nuclei_count"].to_numpy()
        y = frame[column].to_numpy()

        add_identity_line(axis, x, y)

        pearson_r, _ = safe_pearson(x, y)
        spearman_rho, _ = safe_spearman(x, y)

        axis.set_xlabel("Observed nuclei count")
        axis.set_ylabel(f"{label} total inferred abundance")
        axis.set_title(
            f"{label} inferred abundance vs observed nuclei\n"
            f"Pearson r = {pearson_r:.3f}; "
            f"Spearman ρ = {spearman_rho:.3f}"
        )
        axis.legend(frameon=False)

        save_figure(
            figure,
            output_dir / f"{method}_total_vs_nuclei",
        )


def plot_roi_correlation_distributions(
    frame: pd.DataFrame,
    output_dir: Path,
) -> None:
    metrics = (
        (
            "relative_spearman_rho",
            "ROI relative-abundance Spearman correlation",
        ),
        (
            "relative_pearson_r",
            "ROI relative-abundance Pearson correlation",
        ),
        (
            "absolute_spearman_rho",
            "ROI absolute-abundance Spearman correlation",
        ),
        (
            "absolute_pearson_r",
            "ROI absolute-abundance Pearson correlation",
        ),
    )

    for column, title in metrics:
        figure, axis = plt.subplots(figsize=(7, 5))

        values = [
            frame.loc[
                frame["cohort"] == cohort,
                column,
            ].dropna().to_numpy()
            for cohort in COHORTS
        ]

        axis.boxplot(
            values,
            labels=list(COHORTS),
            showfliers=True,
        )

        axis.set_ylabel("Correlation")
        axis.set_title(title)
        axis.axhline(
            0,
            linestyle="--",
            linewidth=1,
        )

        save_figure(
            figure,
            output_dir / f"{column}_by_cohort",
        )


def plot_factor_shift(
    factor_metrics: pd.DataFrame,
    output_dir: Path,
    value_column: str,
    filename: str,
    x_label: str,
    top_n: int = 20,
) -> None:
    summary = (
        factor_metrics
        .groupby("factor", as_index=False)
        .agg(
            mean_value=(value_column, "mean"),
            mean_abs_value=(
                value_column,
                lambda values: float(
                    np.mean(np.abs(values))
                ),
            ),
        )
        .sort_values(
            "mean_abs_value",
            ascending=False,
        )
        .head(top_n)
        .sort_values("mean_value")
    )

    figure_height = max(
        6.0,
        0.34 * len(summary),
    )

    figure, axis = plt.subplots(
        figsize=(8, figure_height)
    )

    axis.barh(
        summary["factor"],
        summary["mean_value"],
    )

    axis.axvline(
        0,
        linestyle="--",
        linewidth=1,
    )

    axis.set_xlabel(x_label)
    axis.set_ylabel("Factor")
    axis.set_title(
        f"Top {len(summary)} factors by absolute shift"
    )

    save_figure(
        figure,
        output_dir / filename,
    )


def plot_factor_correlations(
    factor_metrics: pd.DataFrame,
    output_dir: Path,
) -> None:
    metrics = (
        (
            "relative_spearman_rho",
            "Factor-level relative-abundance correlation",
            "factor_relative_spearman",
        ),
        (
            "absolute_spearman_rho",
            "Factor-level absolute-abundance correlation",
            "factor_absolute_spearman",
        ),
    )

    for column, title, filename in metrics:
        figure, axis = plt.subplots(figsize=(8, 6))

        for cohort in COHORTS:
            subset = factor_metrics.loc[
                factor_metrics["cohort"] == cohort
            ]

            axis.scatter(
                subset["baseline_relative_mean"]
                if "relative" in column
                else subset["baseline_absolute_mean"],
                subset[column],
                alpha=0.8,
                label=cohort,
            )

        axis.axhline(
            0,
            linestyle="--",
            linewidth=1,
        )

        axis.set_xlabel(
            "Baseline mean relative abundance"
            if "relative" in column
            else "Baseline mean absolute abundance"
        )
        axis.set_ylabel("Spearman correlation")
        axis.set_title(title)
        axis.legend(frameon=False)

        save_figure(
            figure,
            output_dir / filename,
        )


def plot_difference_heatmap(
    cohort: str,
    relative_difference: np.ndarray,
    rois: list[str],
    factors: list[str],
    output_dir: Path,
    top_factor_count: int,
    top_roi_count: int,
) -> None:
    """
    Plot the factors and ROIs most affected by the prior change.

    Factors are ranked by their mean absolute difference across ROIs.
    ROIs are ranked by their total absolute difference across factors.
    """

    n_factors = min(
        max(top_factor_count, 1),
        len(factors),
    )

    n_rois = min(
        max(top_roi_count, 1),
        len(rois),
    )

    factor_scores = np.mean(
        np.abs(relative_difference),
        axis=0,
    )

    selected_factor_indices = np.argsort(
        factor_scores
    )[-n_factors:]

    roi_scores = np.sum(
        np.abs(
            relative_difference[
                :,
                selected_factor_indices,
            ]
        ),
        axis=1,
    )

    selected_roi_indices = np.argsort(
        roi_scores
    )[-n_rois:]

    matrix = relative_difference[
        np.ix_(
            selected_roi_indices,
            selected_factor_indices,
        )
    ]

    selected_rois = [
        rois[index]
        for index in selected_roi_indices
    ]

    selected_factors = [
        factors[index]
        for index in selected_factor_indices
    ]

    max_abs = float(np.max(np.abs(matrix)))

    if max_abs == 0:
        max_abs = 1.0

    figure_width = max(
        9.0,
        0.42 * len(selected_factors),
    )

    figure_height = max(
        7.0,
        0.13 * len(selected_rois),
    )

    figure, axis = plt.subplots(
        figsize=(figure_width, figure_height)
    )

    image = axis.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        vmin=-max_abs,
        vmax=max_abs,
    )

    axis.set_xticks(
        np.arange(len(selected_factors))
    )
    axis.set_xticklabels(
        selected_factors,
        rotation=90,
    )

    # ROI labels become unreadable for large heatmaps. Display
    # regularly spaced labels while retaining all rows.
    max_y_labels = 40
    step = max(
        1,
        int(
            math.ceil(
                len(selected_rois) / max_y_labels
            )
        ),
    )

    y_positions = np.arange(
        0,
        len(selected_rois),
        step,
    )

    axis.set_yticks(y_positions)
    axis.set_yticklabels(
        [
            selected_rois[position]
            for position in y_positions
        ]
    )

    axis.set_xlabel("Factor")
    axis.set_ylabel("ROI")
    axis.set_title(
        f"{cohort}: nuclei-informed minus baseline "
        "relative abundance"
    )

    colorbar = figure.colorbar(
        image,
        ax=axis,
        fraction=0.025,
        pad=0.02,
    )
    colorbar.set_label("Relative-abundance difference")

    save_figure(
        figure,
        output_dir
        / f"{cohort}_relative_difference_heatmap",
    )


# ============================================================
# Summary statistics
# ============================================================


def nuclei_association_summary(
    frame: pd.DataFrame,
) -> dict[str, Any]:
    output: dict[str, Any] = {}

    for cohort_label in (*COHORTS, "COMBINED"):
        if cohort_label == "COMBINED":
            subset = frame
        else:
            subset = frame.loc[
                frame["cohort"] == cohort_label
            ]

        nuclei = subset["nuclei_count"].to_numpy()

        baseline = subset[
            "baseline_total_absolute"
        ].to_numpy()

        informed = subset[
            "informed_total_absolute"
        ].to_numpy()

        baseline_pearson, baseline_pearson_p = safe_pearson(
            nuclei,
            baseline,
        )
        baseline_spearman, baseline_spearman_p = safe_spearman(
            nuclei,
            baseline,
        )

        informed_pearson, informed_pearson_p = safe_pearson(
            nuclei,
            informed,
        )
        informed_spearman, informed_spearman_p = safe_spearman(
            nuclei,
            informed,
        )

        baseline_mae = float(
            np.mean(np.abs(baseline - nuclei))
        )
        informed_mae = float(
            np.mean(np.abs(informed - nuclei))
        )

        baseline_rmse = float(
            np.sqrt(
                np.mean((baseline - nuclei) ** 2)
            )
        )
        informed_rmse = float(
            np.sqrt(
                np.mean((informed - nuclei) ** 2)
            )
        )

        output[cohort_label] = {
            "n_rois": int(len(subset)),
            "baseline": {
                "pearson_r": baseline_pearson,
                "pearson_p": baseline_pearson_p,
                "spearman_rho": baseline_spearman,
                "spearman_p": baseline_spearman_p,
                "mae_vs_nuclei": baseline_mae,
                "rmse_vs_nuclei": baseline_rmse,
                "median_total_per_nucleus": float(
                    np.median(
                        baseline / nuclei
                    )
                ),
            },
            "nuclei_informed": {
                "pearson_r": informed_pearson,
                "pearson_p": informed_pearson_p,
                "spearman_rho": informed_spearman,
                "spearman_p": informed_spearman_p,
                "mae_vs_nuclei": informed_mae,
                "rmse_vs_nuclei": informed_rmse,
                "median_total_per_nucleus": float(
                    np.median(
                        informed / nuclei
                    )
                ),
            },
            "change": {
                "pearson_r_difference": (
                    informed_pearson - baseline_pearson
                ),
                "spearman_rho_difference": (
                    informed_spearman - baseline_spearman
                ),
                "mae_reduction": (
                    baseline_mae - informed_mae
                ),
                "rmse_reduction": (
                    baseline_rmse - informed_rmse
                ),
            },
        }

    return output


def print_key_summary(
    roi_metrics: pd.DataFrame,
    factor_metrics: pd.DataFrame,
    nuclei_summary: dict[str, Any],
) -> None:
    print("\n" + "=" * 72)
    print("SPACEJAM NUCLEI-INFORMED PRIOR COMPARISON")
    print("=" * 72)

    for cohort in COHORTS:
        roi_subset = roi_metrics.loc[
            roi_metrics["cohort"] == cohort
        ]

        factor_subset = factor_metrics.loc[
            factor_metrics["cohort"] == cohort
        ]

        print(f"\n{cohort}")
        print("-" * len(cohort))

        print(
            "Median ROI relative Spearman correlation:",
            f"{roi_subset['relative_spearman_rho'].median():.4f}",
        )

        print(
            "Median ROI relative L1 distance:",
            f"{roi_subset['relative_l1_distance'].median():.6f}",
        )

        print(
            "Median total-abundance percent change:",
            f"{roi_subset['total_absolute_percent_change'].median():.2f}%",
        )

        print(
            "Median factor relative Spearman correlation:",
            f"{factor_subset['relative_spearman_rho'].median():.4f}",
        )

        cohort_nuclei = nuclei_summary[cohort]

        print(
            "Baseline total vs nuclei Spearman:",
            f"{cohort_nuclei['baseline']['spearman_rho']:.4f}",
        )

        print(
            "Informed total vs nuclei Spearman:",
            f"{cohort_nuclei['nuclei_informed']['spearman_rho']:.4f}",
        )

        print(
            "Change in total-vs-nuclei Spearman:",
            f"{cohort_nuclei['change']['spearman_rho_difference']:.4f}",
        )

    combined = nuclei_summary["COMBINED"]

    print("\nCOMBINED")
    print("--------")

    print(
        "Baseline total vs nuclei Pearson:",
        f"{combined['baseline']['pearson_r']:.4f}",
    )

    print(
        "Informed total vs nuclei Pearson:",
        f"{combined['nuclei_informed']['pearson_r']:.4f}",
    )

    print(
        "Baseline total vs nuclei Spearman:",
        f"{combined['baseline']['spearman_rho']:.4f}",
    )

    print(
        "Informed total vs nuclei Spearman:",
        f"{combined['nuclei_informed']['spearman_rho']:.4f}",
    )

    print(
        "MAE reduction relative to nuclei:",
        f"{combined['change']['mae_reduction']:.4f}",
    )

    print(
        "RMSE reduction relative to nuclei:",
        f"{combined['change']['rmse_reduction']:.4f}",
    )


# ============================================================
# Main
# ============================================================


def main() -> None:
    args = parse_args()

    project_dir = args.project_dir.resolve()

    baseline_dir = (
        args.baseline_dir.resolve()
        if args.baseline_dir is not None
        else project_dir / "results" / "spacejam"
    )

    informed_dir = (
        args.nuclei_informed_dir.resolve()
        if args.nuclei_informed_dir is not None
        else (
            project_dir
            / "results"
            / "spacejam_nuclei_informed"
        )
    )

    nuclei_file = (
        args.nuclei_file.resolve()
        if args.nuclei_file is not None
        else project_dir / "data" / "nuc_count.csv"
    )

    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else (
            project_dir
            / "results"
            / "spacejam_nuclei_prior_comparison"
        )
    )

    table_dir = output_dir / "tables"
    figure_dir = output_dir / "figures"

    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    ensure_directory(baseline_dir)
    ensure_directory(informed_dir)
    ensure_file(nuclei_file)

    print("Project directory:", project_dir)
    print("Baseline directory:", baseline_dir)
    print("Nuclei-informed directory:", informed_dir)
    print("Nuclei-count file:", nuclei_file)
    print("Output directory:", output_dir)

    nuclei_counts = read_nuclei_counts(nuclei_file)

    all_roi_metrics: list[pd.DataFrame] = []
    all_factor_metrics: list[pd.DataFrame] = []
    all_long_tables: list[pd.DataFrame] = []

    relative_differences: dict[str, np.ndarray] = {}
    cohort_manifests: dict[str, dict[str, list[str]]] = {}

    for cohort in COHORTS:
        print(f"\nLoading {cohort}...")

        baseline = load_cohort(
            baseline_dir,
            cohort,
        )

        informed = load_cohort(
            informed_dir,
            cohort,
        )

        compare_manifests(
            baseline,
            informed,
            cohort,
        )

        rois = baseline["rois"]
        factors = baseline["factors"]

        print(
            f"{cohort}: {len(rois)} ROIs × "
            f"{len(factors)} factors"
        )

        roi_metrics = build_roi_metrics(
            cohort=cohort,
            rois=rois,
            baseline_absolute=baseline["absolute"],
            informed_absolute=informed["absolute"],
            baseline_relative=baseline["relative"],
            informed_relative=informed["relative"],
        )

        factor_metrics = build_factor_metrics(
            cohort=cohort,
            factors=factors,
            baseline_absolute=baseline["absolute"],
            informed_absolute=informed["absolute"],
            baseline_relative=baseline["relative"],
            informed_relative=informed["relative"],
        )

        long_table = build_long_abundance_table(
            cohort=cohort,
            rois=rois,
            factors=factors,
            baseline_absolute=baseline["absolute"],
            informed_absolute=informed["absolute"],
            baseline_relative=baseline["relative"],
            informed_relative=informed["relative"],
        )

        all_roi_metrics.append(roi_metrics)
        all_factor_metrics.append(factor_metrics)
        all_long_tables.append(long_table)

        relative_differences[cohort] = (
            informed["relative"]
            - baseline["relative"]
        )

        cohort_manifests[cohort] = {
            "rois": rois,
            "factors": factors,
        }

    roi_metrics = pd.concat(
        all_roi_metrics,
        ignore_index=True,
    )

    roi_metrics = add_nuclei_counts(
        roi_metrics,
        nuclei_counts,
    )

    factor_metrics = pd.concat(
        all_factor_metrics,
        ignore_index=True,
    )

    long_abundance = pd.concat(
        all_long_tables,
        ignore_index=True,
    )

    # --------------------------------------------------------
    # Save tabular outputs
    # --------------------------------------------------------

    roi_metrics.to_csv(
        table_dir / "roi_level_comparison.csv",
        index=False,
    )

    factor_metrics.to_csv(
        table_dir / "factor_level_comparison.csv",
        index=False,
    )

    long_abundance.to_csv(
        table_dir / "roi_factor_abundance_changes.csv.gz",
        index=False,
        compression="gzip",
    )

    roi_metrics.sort_values(
        "relative_l1_distance",
        ascending=False,
    ).to_csv(
        table_dir
        / "rois_ranked_by_relative_composition_change.csv",
        index=False,
    )

    factor_metrics.sort_values(
        "relative_mean_absolute_difference",
        ascending=False,
    ).to_csv(
        table_dir
        / "factors_ranked_by_relative_change.csv",
        index=False,
    )

    factor_metrics.sort_values(
        "absolute_mean_absolute_difference",
        ascending=False,
    ).to_csv(
        table_dir
        / "factors_ranked_by_absolute_change.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    nuclei_summary = nuclei_association_summary(
        roi_metrics
    )

    summary = {
        "inputs": {
            "baseline_dir": str(baseline_dir),
            "nuclei_informed_dir": str(informed_dir),
            "nuclei_file": str(nuclei_file),
        },
        "dimensions": {
            cohort: {
                "n_rois": len(
                    cohort_manifests[cohort]["rois"]
                ),
                "n_factors": len(
                    cohort_manifests[cohort]["factors"]
                ),
            }
            for cohort in COHORTS
        },
        "roi_relative_spearman": {
            cohort: summarize_numeric(
                roi_metrics.loc[
                    roi_metrics["cohort"] == cohort,
                    "relative_spearman_rho",
                ]
            )
            for cohort in COHORTS
        },
        "roi_relative_l1_distance": {
            cohort: summarize_numeric(
                roi_metrics.loc[
                    roi_metrics["cohort"] == cohort,
                    "relative_l1_distance",
                ]
            )
            for cohort in COHORTS
        },
        "roi_total_absolute_percent_change": {
            cohort: summarize_numeric(
                roi_metrics.loc[
                    roi_metrics["cohort"] == cohort,
                    "total_absolute_percent_change",
                ]
            )
            for cohort in COHORTS
        },
        "factor_relative_spearman": {
            cohort: summarize_numeric(
                factor_metrics.loc[
                    factor_metrics["cohort"] == cohort,
                    "relative_spearman_rho",
                ]
            )
            for cohort in COHORTS
        },
        "nuclei_association": nuclei_summary,
    }

    with open(
        output_dir / "comparison_summary.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            json_safe(summary),
            handle,
            indent=2,
        )

    # Save a compact nuclei-association table.
    nuclei_rows: list[dict[str, Any]] = []

    for cohort, values in nuclei_summary.items():
        nuclei_rows.append(
            {
                "cohort": cohort,
                "n_rois": values["n_rois"],
                "baseline_pearson_r": (
                    values["baseline"]["pearson_r"]
                ),
                "informed_pearson_r": (
                    values["nuclei_informed"]["pearson_r"]
                ),
                "pearson_r_difference": (
                    values["change"][
                        "pearson_r_difference"
                    ]
                ),
                "baseline_spearman_rho": (
                    values["baseline"]["spearman_rho"]
                ),
                "informed_spearman_rho": (
                    values["nuclei_informed"][
                        "spearman_rho"
                    ]
                ),
                "spearman_rho_difference": (
                    values["change"][
                        "spearman_rho_difference"
                    ]
                ),
                "baseline_mae_vs_nuclei": (
                    values["baseline"]["mae_vs_nuclei"]
                ),
                "informed_mae_vs_nuclei": (
                    values["nuclei_informed"][
                        "mae_vs_nuclei"
                    ]
                ),
                "mae_reduction": (
                    values["change"]["mae_reduction"]
                ),
                "baseline_rmse_vs_nuclei": (
                    values["baseline"]["rmse_vs_nuclei"]
                ),
                "informed_rmse_vs_nuclei": (
                    values["nuclei_informed"][
                        "rmse_vs_nuclei"
                    ]
                ),
                "rmse_reduction": (
                    values["change"]["rmse_reduction"]
                ),
            }
        )

    pd.DataFrame(nuclei_rows).to_csv(
        table_dir / "total_abundance_vs_nuclei_summary.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Figures
    # --------------------------------------------------------

    plot_total_baseline_vs_informed(
        roi_metrics,
        figure_dir,
    )

    plot_totals_vs_nuclei(
        roi_metrics,
        figure_dir,
    )

    plot_roi_correlation_distributions(
        roi_metrics,
        figure_dir,
    )

    plot_factor_shift(
        factor_metrics=factor_metrics,
        output_dir=figure_dir,
        value_column="relative_mean_difference",
        filename="factor_mean_relative_abundance_shift",
        x_label=(
            "Mean relative-abundance difference "
            "(nuclei-informed − baseline)"
        ),
    )

    plot_factor_shift(
        factor_metrics=factor_metrics,
        output_dir=figure_dir,
        value_column="absolute_mean_difference",
        filename="factor_mean_absolute_abundance_shift",
        x_label=(
            "Mean absolute-abundance difference "
            "(nuclei-informed − baseline)"
        ),
    )

    plot_factor_correlations(
        factor_metrics,
        figure_dir,
    )

    for cohort in COHORTS:
        plot_difference_heatmap(
            cohort=cohort,
            relative_difference=(
                relative_differences[cohort]
            ),
            rois=cohort_manifests[cohort]["rois"],
            factors=(
                cohort_manifests[cohort]["factors"]
            ),
            output_dir=figure_dir,
            top_factor_count=args.heatmap_top_factors,
            top_roi_count=args.heatmap_top_rois,
        )

    print_key_summary(
        roi_metrics,
        factor_metrics,
        nuclei_summary,
    )

    print("\nOutputs written to:")
    print(output_dir)
    print("\nMain tables:")
    print(
        table_dir / "roi_level_comparison.csv"
    )
    print(
        table_dir / "factor_level_comparison.csv"
    )
    print(
        table_dir
        / "total_abundance_vs_nuclei_summary.csv"
    )
    print("\nSummary:")
    print(output_dir / "comparison_summary.json")


if __name__ == "__main__":
    main()

