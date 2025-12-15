#!/usr/bin/env bash
set -euo pipefail

set -x
trap 'echo "[ERROR] Failed at line $LINENO: $BASH_COMMAND" >&2' ERR


# =========================
# cicero_scopen.sh
#   Run ONCE: Cicero -> scOpen -> clustering -> metrics
#
# Usage:
#   ./cicero_scopen.sh DATA_NAME ORGANISM CICERO_TOPN
# Example:
#   ./cicero_scopen.sh PBMC5k human 10000
# =========================

DATA_NAME="${1:?Usage: $0 DATA_NAME ORGANISM CICERO_TOPN}"
ORGANISM="${2:?Usage: $0 DATA_NAME ORGANISM CICERO_TOPN}"
CICERO_TOPN="${3:?Usage: $0 DATA_NAME ORGANISM CICERO_TOPN}"

HOME_DIR="${HOME:?HOME is not set}"
DATA_DIR="${HOME_DIR}/Datasets/${DATA_NAME}"
INPUT_DIR="${DATA_DIR}/input"
OUTPUT_DIR="${DATA_DIR}/output"

SAPIENS_DIR="${HOME_DIR}/mySAPIEnS"
PREP_DIR="${SAPIENS_DIR}/preprocessing"
IMP_DIR="${SAPIENS_DIR}/imputation"
CLUSTER_DIR="${SAPIENS_DIR}/clustering"
METRICS_SCRIPT="${CLUSTER_DIR}/make_metrics_cicero_scopen.py"

LABELS_TSV="${INPUT_DIR}/labels.tsv"

LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/cicero_scopen_${DATA_NAME}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[INFO] DATA_NAME=${DATA_NAME}"
echo "[INFO] ORGANISM=${ORGANISM}"
echo "[INFO] CICERO_TOPN=${CICERO_TOPN}"
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
need_dir "${INPUT_DIR}"
need_dir "${OUTPUT_DIR}"
need_file "${LABELS_TSV}"
need_dir "${SAPIENS_DIR}"
need_dir "${IMP_DIR}"
need_dir "${PREP_DIR}"
need_dir "${CLUSTER_DIR}"
need_file "${METRICS_SCRIPT}"

# ============================================================
# 1) Cicero
# ============================================================
activate_conda "sapiens"

(
  cd "${PREP_DIR}"
  echo "[CMD] ./run_cicero_pipeline.sh '${DATA_DIR}' '${ORGANISM}' '${CICERO_TOPN}'"
  ./run_cicero_pipeline.sh "${DATA_DIR}" "${ORGANISM}" "${CICERO_TOPN}"
)

# prepare 10X-like files for scOpen in output/cicero
CICERO_DIR="${OUTPUT_DIR}/cicero"
need_dir "${CICERO_DIR}"
need_file "${CICERO_DIR}/counts.mtx"
need_file "${CICERO_DIR}/barcodes.txt"
need_file "${CICERO_DIR}/peaks.txt"

(
  cd "${CICERO_DIR}"
  echo "[CMD] cp -f counts.mtx matrix.mtx"
  cp -f counts.mtx matrix.mtx

  echo "[CMD] cp -f barcodes.txt barcodes.tsv"
  cp -f barcodes.txt barcodes.tsv

  echo "[CMD] peaks.txt -> peaks.bed (underscore to tab)"
  awk -F'_' 'BEGIN{OFS="\t"} {print $1,$2,$3}' peaks.txt > peaks.bed
)

need_file "${CICERO_DIR}/matrix.mtx"
need_file "${CICERO_DIR}/barcodes.tsv"
need_file "${CICERO_DIR}/peaks.bed"

# ============================================================
# 2) scOpen (baseline method = cicero)
# ============================================================
activate_conda "scopen"
(
  cd "${IMP_DIR}"
  echo "[CMD] ./02_run_scopen_full.sh '${DATA_NAME}' cicero"
  ./02_run_scopen_full.sh "${DATA_NAME}" cicero
)

# ============================================================
# 3) clustering + metrics
# ============================================================
activate_conda "clustering_sapiens"
CICERO_SCOPEN_DIR="${OUTPUT_DIR}/cicero_scopen"
MATRICES_DIR="${CICERO_SCOPEN_DIR}/matrices"
need_dir "${CICERO_SCOPEN_DIR}"
need_dir "${MATRICES_DIR}"

(
  cd "${CLUSTER_DIR}"
  echo "[CMD] ./run_scopen_clustering.sh '${MATRICES_DIR}' '${CICERO_SCOPEN_DIR}/' '${LABELS_TSV}'"
  ./run_scopen_clustering.sh "${MATRICES_DIR}" "${CICERO_SCOPEN_DIR}/" "${LABELS_TSV}"
)

echo "[CMD] python '${METRICS_SCRIPT}' --base_dir '${CICERO_SCOPEN_DIR}/'"
python "${METRICS_SCRIPT}" --base_dir "${CICERO_SCOPEN_DIR}/"

echo "[DONE] Baseline finished: ${CICERO_SCOPEN_DIR}"

