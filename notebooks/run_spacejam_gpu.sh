#!/bin/bash
#SBATCH -J spacejam_nuclei
#SBATCH -p gpu
#SBATCH -A r01604
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration/results/Script_results/spacejam_nuclei_%j.out
#SBATCH -e /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration/results/Script_results/spacejam_nuclei_%j.err

set -eo pipefail

# Module/Conda initialization may expect PS1 in a batch shell.
export PS1=""

echo "====================================="
echo "Nuclei-informed SpaceJam sensitivity run"
echo "Job started on: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-unknown}"
echo "Start time: $(date)"
echo "====================================="

PROJECT_DIR="/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration"

SCRIPT_RESULTS_DIR="${PROJECT_DIR}/results/Script_results"
SPACEJAM_DIR="${PROJECT_DIR}/results/spacejam"
REGRESSION_DIR="${PROJECT_DIR}/results/regression_model"
NUCLEI_FILE="${PROJECT_DIR}/data/nuc_count.csv"

mkdir -p "${SCRIPT_RESULTS_DIR}"
mkdir -p "${SPACEJAM_DIR}"

cd "${PROJECT_DIR}"

# Make the output location available to the Python runner.
export SPACEJAM_RESULTS_DIR="${SPACEJAM_DIR}"

# ------------------------------------------------------------
# Environment
# ------------------------------------------------------------

module purge
module load gnu/12.2.0
module load sqlite/3.35.5
module load python/gpu/3.10.10

module list

source /N/soft/rhel8/conda/26.3.2/etc/profile.d/conda.sh

conda activate \
  /N/u/echimal/Quartz/.conda/envs/integration_env

export PYTHONPATH="${PROJECT_DIR}"

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16

echo
echo "Python executable: $(command -v python)"
python --version
echo "Conda environment: ${CONDA_PREFIX:-unknown}"
echo "Analysis output: ${SPACEJAM_DIR}"
echo

# ------------------------------------------------------------
# GPU check
# ------------------------------------------------------------

python - <<'PY'
import torch

print("CUDA available:", torch.cuda.is_available())

if not torch.cuda.is_available():
    raise RuntimeError(
        "GPU was requested, but CUDA is not available."
    )

print("GPU:", torch.cuda.get_device_name(0))
print("CUDA device count:", torch.cuda.device_count())
PY

# ------------------------------------------------------------
# Syntax and import checks
# ------------------------------------------------------------

echo
echo "Checking model and runner syntax..."

python -m py_compile \
  models/LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified.py \
  notebooks/03_spacejam_cell2location.py

python - <<'PY'
from models.LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified import (
    LocationModelPyro,
)

print("LocationModelPyro import successful.")
PY

# ------------------------------------------------------------
# Input checks
# ------------------------------------------------------------

echo
echo "Checking required inputs..."

test -f "${REGRESSION_DIR}/AD+CAA_inferred_signatures.csv"
test -f "${REGRESSION_DIR}/Control_inferred_signatures.csv"
test -f "${NUCLEI_FILE}"

python - <<'PY'
from pathlib import Path

import pandas as pd

project = Path(
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/multiomic-spatial-integration"
)

signature_dir = project / "results" / "regression_model"
nuclei_file = project / "data" / "nuc_count.csv"

ad = pd.read_csv(
    signature_dir / "AD+CAA_inferred_signatures.csv",
    index_col=0,
)

ctrl = pd.read_csv(
    signature_dir / "Control_inferred_signatures.csv",
    index_col=0,
)

print("AD+CAA signature shape:", ad.shape)
print("Control signature shape:", ctrl.shape)
print(
    "Same factor set:",
    set(ad.columns) == set(ctrl.columns),
)
print(
    "Same factor order:",
    ad.columns.tolist() == ctrl.columns.tolist(),
)

if ad.columns.tolist() != ctrl.columns.tolist():
    raise RuntimeError(
        "AD+CAA and Control signature columns are not "
        "in the same order."
    )

if len(ad.columns) != 46:
    raise RuntimeError(
        f"Expected 46 factors, found {len(ad.columns)}."
    )

nuclei = pd.read_csv(nuclei_file)

required = {"ROI_ID", "nuclei_count"}

missing = required - set(nuclei.columns)

if missing:
    raise RuntimeError(
        "Nuclei file is missing columns: "
        f"{sorted(missing)}"
    )

if nuclei["ROI_ID"].duplicated().any():
    raise RuntimeError(
        "Duplicated ROI_ID values were found in nuc_count.csv."
    )

if nuclei["nuclei_count"].isna().any():
    raise RuntimeError(
        "Missing nuclei_count values were found."
    )

if (nuclei["nuclei_count"] <= 0).any():
    raise RuntimeError(
        "All nuclei counts must be positive."
    )

print("Nuclei rows:", len(nuclei))
print(
    "Nuclei range:",
    nuclei["nuclei_count"].min(),
    "to",
    nuclei["nuclei_count"].max(),
)
print(
    "Median nuclei count:",
    nuclei["nuclei_count"].median(),
)

print("Input validation passed.")
PY

# ------------------------------------------------------------
# Protect against accidental overwrite
# ------------------------------------------------------------

if find "${SPACEJAM_DIR}" \
    -mindepth 1 \
    -maxdepth 1 \
    -print \
    -quit | grep -q .
then

    STAMP=$(date +%Y%m%d_%H%M%S)

    BACKUP_DIR="${PROJECT_DIR}/results/archive/spacejam_nuclei_informed_${STAMP}"

    mkdir -p "${BACKUP_DIR}"

    cp -a \
      "${SPACEJAM_DIR}/." \
      "${BACKUP_DIR}/"

    echo "Existing sensitivity outputs backed up to:"
    echo "${BACKUP_DIR}"
fi

# ------------------------------------------------------------
# Training
# ------------------------------------------------------------

echo
echo "Starting nuclei-informed SpaceJam training..."
echo "Start time: $(date)"

python \
  notebooks/03_spacejam_cell2location.py

# ------------------------------------------------------------
# Output checks
# ------------------------------------------------------------

echo
echo "Checking expected outputs..."

expected_files=(
    "${SPACEJAM_DIR}/ADCAA_spot_factors_abs.pt"
    "${SPACEJAM_DIR}/ADCAA_spot_factors_rel.pt"
    "${SPACEJAM_DIR}/CTRL_spot_factors_abs.pt"
    "${SPACEJAM_DIR}/CTRL_spot_factors_rel.pt"
    "${SPACEJAM_DIR}/ADCAA_param_store.pt"
    "${SPACEJAM_DIR}/CTRL_param_store.pt"
    "${SPACEJAM_DIR}/ADCAA_training_loss.csv"
    "${SPACEJAM_DIR}/CTRL_training_loss.csv"
    "${SPACEJAM_DIR}/ADCAA_manifest_rois.csv"
    "${SPACEJAM_DIR}/ADCAA_manifest_factors.csv"
    "${SPACEJAM_DIR}/ADCAA_manifest_experiments.csv"
    "${SPACEJAM_DIR}/CTRL_manifest_rois.csv"
    "${SPACEJAM_DIR}/CTRL_manifest_factors.csv"
    "${SPACEJAM_DIR}/CTRL_manifest_experiments.csv"
    "${SPACEJAM_DIR}/ADCAA_run_metadata.json"
    "${SPACEJAM_DIR}/CTRL_run_metadata.json"
)

for file in "${expected_files[@]}"; do
    if [[ ! -f "${file}" ]]; then
        echo "Missing expected output: ${file}" >&2
        exit 1
    fi
done

echo
echo "All expected nuclei-informed SpaceJam outputs were found."

echo "====================================="
echo "Job finished successfully"
echo "End time: $(date)"
echo "Results: ${SPACEJAM_DIR}"
echo "====================================="
