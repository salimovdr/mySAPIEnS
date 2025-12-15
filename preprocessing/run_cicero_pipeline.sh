#!/usr/bin/env bash
#SBATCH --partition=comet
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

set -euo pipefail

set -x
trap 'echo "[ERROR] run_cicero_pipeline.sh failed at line $LINENO: $BASH_COMMAND" >&2' ERR

DATASET_DIR="${1:?Usage: $0 DATASET_DIR [ORGANISM] [NUM_REMAIN]}"
ORGANISM="${2:-human}"
NUM_REMAIN="${3:-50000}"

# Expected layout:
#   ${DATASET_DIR}/input
#   ${DATASET_DIR}/output
INPUT_DIR="${DATASET_DIR}/input"
OUTPUT_DIR="${DATASET_DIR}/output"

CICERO_FOLDER="${OUTPUT_DIR}/cicero"

# Resolve script directory to make relative calls robust
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${CICERO_FOLDER}"

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "[ERROR] INPUT_DIR not found: ${INPUT_DIR}" >&2
  exit 2
fi

echo "[INFO] DATASET_DIR=${DATASET_DIR}"
echo "[INFO] INPUT_DIR=${INPUT_DIR}"
echo "[INFO] OUTPUT_DIR=${OUTPUT_DIR}"
echo "[INFO] CICERO_FOLDER=${CICERO_FOLDER}"
echo "[INFO] ORGANISM=${ORGANISM}"
echo "[INFO] NUM_REMAIN=${NUM_REMAIN}"

python "${SCRIPT_DIR}/extract_cicero_regions_original.py" \
  --folder "${INPUT_DIR}" \
  --output "${CICERO_FOLDER}/peaks_dumped.tsv"

bash "${SCRIPT_DIR}/split_dataset.sh" \
  "${CICERO_FOLDER}/peaks_dumped.tsv" \
  "${CICERO_FOLDER}/peaks" \
  "${ORGANISM}"

mkdir -p "${CICERO_FOLDER}/filtered"

bash "${SCRIPT_DIR}/run_cicero.sh" \
  "${CICERO_FOLDER}/peaks_" \
  "${CICERO_FOLDER}/filtered/peaks_" \
  "${HOME}/Datasets/hg38.chrom.sizes"

python "${SCRIPT_DIR}/filter_cells_by_coaccess_count_mtx.py" \
  --input "${INPUT_DIR}" \
  --prefix "${CICERO_FOLDER}/filtered/peaks_" \
  --organism "${ORGANISM}" \
  --output "${CICERO_FOLDER}" \
  --remain "${NUM_REMAIN}"

