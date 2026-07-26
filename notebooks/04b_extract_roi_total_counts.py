# 04b — Extract Per-ROI Total Counts and Size Factors

# Separate from 04_extract_cell_proportions.py: that notebook needs both
# trained param stores; this one only needs the raw h5ad, so it can run
# independently and sooner. Produces the per-ROI total-count offset used by
# the absolute-abundance contrast pipeline (R/absolute_abundance_utils.R)
# and the calibration gate (R/absolute_calibration_utils.R).


# This script reads the raw GeoMx AnnData object, computes ROI-level total-count
# information using `compute_roi_total_counts()`, and saves the resulting table
# for the absolute-abundance contrast and calibration workflows.
# """

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import scanpy as sc

from python.abundance_extraction_utils import compute_roi_total_counts

INPUT_PATH = REPO_ROOT / "data" / "CAA-AD_AnnData.h5ad"
OUTPUT_PATH = (
    REPO_ROOT
    / "results"
    / "cell_proportions"
    / "roi_total_counts.csv"
)


def main() -> None:
    if not INPUT_PATH.exists():
        raise FileNotFoundError(
            f"Input AnnData file not found: {INPUT_PATH}"
        )

    print(f"Reading AnnData: {INPUT_PATH}")
    adata = sc.read_h5ad(INPUT_PATH)

    if adata.n_obs == 0:
        raise ValueError("The AnnData object contains no ROIs.")

    if adata.n_vars == 0:
        raise ValueError("The AnnData object contains no genes.")

    roi_total_counts = compute_roi_total_counts(adata)

    if roi_total_counts.empty:
        raise ValueError(
            "compute_roi_total_counts() returned an empty table."
        )

    required_columns = {"n_genes_used_for_size_factor"}
    missing_columns = required_columns.difference(
        roi_total_counts.columns
    )

    if missing_columns:
        raise ValueError(
            "ROI total-count table is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    n_genes_used = roi_total_counts[
        "n_genes_used_for_size_factor"
    ].dropna().unique()

    if len(n_genes_used) != 1:
        raise ValueError(
            "Expected one consistent value for "
            "'n_genes_used_for_size_factor', but found: "
            f"{n_genes_used.tolist()}"
        )

    print("\nROI total counts preview:")
    print(roi_total_counts.head())

    print("\nSummary statistics:")
    print(
        roi_total_counts[
            [
                "total_counts",
                "n_genes_detected",
                "size_factor_mor",
            ]
        ].describe()
    )

    print(
        "\nn_genes_used_for_size_factor:",
        int(n_genes_used[0]),
        "/",
        adata.n_vars,
    )

    output_dir = REPO_ROOT / "results" / "cell_proportions"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "roi_total_counts.csv"

    roi_total_counts.to_csv(
        output_file,
        index=False,
    )

    print(f"\nSaved: {output_file}")


if __name__ == "__main__":
    main()
