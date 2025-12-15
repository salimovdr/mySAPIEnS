#!/usr/bin/env bash
set -euo pipefail

# =========================
# celltypist_cicero_scopen.sh
#   Re-run many times: CellTypist -> core regions -> intersect -> Cicero coaccess expansion
#   -> reduce matrix -> scOpen -> clustering -> metrics
#
# It DOES NOT recompute Cicero baseline.
# It assumes Cicero results already exist in:
#   ~/Datasets/{DATA_NAME}/output/cicero/...
#
# Usage:
#   ./celltypist_cicero_scopen.sh DATA_NAME MODEL_NAME TOP_N_PER_CLASS WINDOW COACC_THRESH
#
# Example:
#   ./celltypist_cicero_scopen.sh PBMC10k Healthy_COVID19_PBMC 300 2000 0.1
# =========================

DATA_NAME="${1:?Usage: $0 DATA_NAME MODEL_NAME TOP_N_PER_CLASS WINDOW COACC_THRESH}"
MODEL_NAME="${2:?Usage: $0 DATA_NAME MODEL_NAME TOP_N_PER_CLASS WINDOW COACC_THRESH}"
TOP_N="${3:?Usage: $0 DATA_NAME MODEL_NAME TOP_N_PER_CLASS WINDOW COACC_THRESH}"
WINDOW="${4:?Usage: $0 DATA_NAME MODEL_NAME TOP_N_PER_CLASS WINDOW COACC_THRESH}"
COACC_THRESH="${5:?Usage: $0 DATA_NAME MODEL_NAME TOP_N_PER_CLASS WINDOW COACC_THRESH}"

HOME_DIR="${HOME:?HOME is not set}"
DATA_DIR="${HOME_DIR}/Datasets/${DATA_NAME}"
INPUT_DIR="${DATA_DIR}/input"
OUTPUT_DIR="${DATA_DIR}/output"

SAPIENS_DIR="${HOME_DIR}/mySAPIEnS"
IMP_DIR="${SAPIENS_DIR}/imputation"
CLUSTER_DIR="${SAPIENS_DIR}/clustering"
PREP_DIR="${SAPIENS_DIR}/preprocessing"
METRICS_SCRIPT="${CLUSTER_DIR}/make_metrics_cicero_scopen.py"

# Your python scripts
CELLTYPIST_DIR="${SAPIENS_DIR}/celltypist"
CORE_SCRIPT="${CELLTYPIST_DIR}/core_genes_and_regions.py"  
FILTER_SCRIPT="${CELLTYPIST_DIR}/cicero_filtering.py"

# Intersect script you already use
INTERSECT_SCRIPT="${CELLTYPIST_DIR}/intersect_core_peaks.sh"

LABELS_TSV="${INPUT_DIR}/labels.tsv"

LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/celltypist_cicero_scopen_${DATA_NAME}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[INFO] DATA_NAME=${DATA_NAME}"
echo "[INFO] MODEL_NAME=${MODEL_NAME}"
echo "[INFO] TOP_N=${TOP_N}"
echo "[INFO] WINDOW=${WINDOW}"
echo "[INFO] COACC_THRESH=${COACC_THRESH}"
echo "[INFO] LOG_FILE=${LOG_FILE}"

need_file() { [[ -s "$1" ]] || { echo "[ERROR] Missing/empty file: $1"; exit 2; }; }
need_dir()  { [[ -d "$1" ]] || { echo "[ERROR] Missing dir: $1"; exit 2; }; }

activate_conda() {
  local env="$1"
  command -v conda >/dev/null 2>&1 || { echo "[ERROR] conda not found"; exit 2; }
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$env"
  echo "[INFO] conda env: $env"
}

# ---- checks ----
need_dir "${OUTPUT_DIR}"
need_file "${LABELS_TSV}"
need_dir "${SAPIENS_DIR}"
need_dir "${CLUSTER_DIR}"
need_dir "${PREP_DIR}"
need_file "${METRICS_SCRIPT}"
need_dir "${IMP_DIR}"

# must exist from baseline Cicero run
CICERO_DIR="${OUTPUT_DIR}/cicero"
need_dir "${CICERO_DIR}"
need_file "${CICERO_DIR}/peaks.txt"
need_file "${CICERO_DIR}/matrix.mtx"
need_file "${CICERO_DIR}/barcodes.tsv"

# scripts
need_file "${CORE_SCRIPT}"
need_file "${FILTER_SCRIPT}"
need_file "${INTERSECT_SCRIPT}"

mkdir -p "${CELLTYPIST_DIR}"
mkdir -p "${OUTPUT_DIR}/celltypist_cicero"

# ============================================================
# 1) CellTypist: core genes + regions around TSS (BED)
# ============================================================
activate_conda "celltypist"
echo "[CMD] python '${CORE_SCRIPT}' --data_name '${DATA_NAME}' --model_name '${MODEL_NAME}' --top_n_genes_per_class '${TOP_N}' --window '${WINDOW}'"
python "${CORE_SCRIPT}" \
  --data_name "${DATA_NAME}" \
  --model_name "${MODEL_NAME}" \
  --top_n_genes_per_class "${TOP_N}" \
  --window "${WINDOW}"

# ============================================================
# 2) Intersect: core regions BED -> intersected peaks (txt)
#    (your existing script)
# ============================================================
echo "[CMD] '${INTERSECT_SCRIPT}' '${DATA_NAME}' '${WINDOW}'"
bash "${INTERSECT_SCRIPT}" "${DATA_NAME}" "${WINDOW}"

CORE_INTERSECT_TXT="${OUTPUT_DIR}/celltypist/core_genes_regions_${WINDOW}bp_intersected_peaks.txt"
need_file "${CORE_INTERSECT_TXT}"

# ============================================================
# 3) Cicero filtering/expansion by coaccess + reduce matrix for scOpen
# ============================================================
echo "[CMD] python '${FILTER_SCRIPT}' --data_name '${DATA_NAME}' --window '${WINDOW}' --coacc_thresh '${COACC_THRESH}'"
python "${FILTER_SCRIPT}" \
  --data_name "${DATA_NAME}" \
  --window "${WINDOW}" \
  --coacc_thresh "${COACC_THRESH}"

# must produce scOpen input here:
CELLTYPIST_CICERO_DIR="${OUTPUT_DIR}/celltypist_cicero"
need_dir "${CELLTYPIST_CICERO_DIR}"
need_file "${CELLTYPIST_CICERO_DIR}/matrix.mtx"
need_file "${CELLTYPIST_CICERO_DIR}/barcodes.tsv"
need_file "${CELLTYPIST_CICERO_DIR}/peaks.bed"

# ============================================================
# 4) scOpen (method = celltypist_cicero)
# ============================================================
activate_conda "scopen"
(
  cd "${IMP_DIR}"
  echo "[CMD] ./02_run_scopen_full.sh '${DATA_NAME}' celltypist_cicero"
  ./02_run_scopen_full.sh "${DATA_NAME}" celltypist_cicero
)

# ============================================================
# 5) clustering + metrics
# ============================================================
activate_conda "clustering_sapiens"
CELLT_SCOPEN_DIR="${OUTPUT_DIR}/celltypist_cicero_scopen"
MATRICES_DIR="${CELLT_SCOPEN_DIR}/matrices"
need_dir "${CELLT_SCOPEN_DIR}"
need_dir "${MATRICES_DIR}"

(
  cd "${CLUSTER_DIR}"
  echo "[CMD] ./run_scopen_clustering.sh '${MATRICES_DIR}' '${CELLT_SCOPEN_DIR}/' '${LABELS_TSV}'"
  ./run_scopen_clustering.sh "${MATRICES_DIR}" "${CELLT_SCOPEN_DIR}/" "${LABELS_TSV}"
)

echo "[CMD] python '${METRICS_SCRIPT}' --base_dir '${CELLT_SCOPEN_DIR}/'"
python "${METRICS_SCRIPT}" --base_dir "${CELLT_SCOPEN_DIR}/"

echo "[DONE] CellTypist branch finished: ${CELLT_SCOPEN_DIR}"

