#!/usr/bin/env bash
set -euo pipefail

# ==========================
# Проверка аргументов
# ==========================
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <DATA_NAME> <WINDOW>"
    echo "Example: $0 PBMC10k 2000"
    exit 1
fi

DATA_NAME="$1"
WINDOW="$2"

# ==========================
# HOME
# ==========================
if [[ -z "${HOME:-}" ]]; then
    echo "[ERROR] HOME is not set"
    exit 1
fi

DATASET_DIR="${HOME}/Datasets/${DATA_NAME}"

echo "[INFO] Dataset: ${DATA_NAME}"
echo "[INFO] Window:  ${WINDOW} bp"
echo "[INFO] Dataset dir: ${DATASET_DIR}"

# ==========================
# Пути
# ==========================
PEAKS="${DATASET_DIR}/input/peaks.bed"
CORE_BED="${DATASET_DIR}/output/celltypist/core_genes_regions_${WINDOW}bp.bed"

OUT_BED="${DATASET_DIR}/output/celltypist/core_genes_regions_${WINDOW}bp_intersected_peaks.bed"
OUT_TXT="${DATASET_DIR}/output/celltypist/core_genes_regions_${WINDOW}bp_intersected_peaks.txt"

# ==========================
# Проверка файлов
# ==========================
if [[ ! -f "${PEAKS}" ]]; then
    echo "[ERROR] peaks file not found: ${PEAKS}"
    exit 1
fi

if [[ ! -f "${CORE_BED}" ]]; then
    echo "[ERROR] core gene regions file not found: ${CORE_BED}"
    exit 1
fi

mkdir -p "$(dirname "${OUT_BED}")"

# ==========================
# BEDTOOLS INTERSECT
# ==========================
echo "[INFO] Running bedtools intersect"
echo "  -a ${PEAKS}"
echo "  -b ${CORE_BED}"
echo "  -> ${OUT_BED}"

bedtools intersect \
    -a "${PEAKS}" \
    -b "${CORE_BED}" \
    -wa -u > "${OUT_BED}"

# ==========================
# BED → TXT (underscore-separated)
# ==========================
echo "[INFO] Converting BED to TXT (tab → underscore)"
awk '{gsub(/\t/, "_"); print}' "${OUT_BED}" > "${OUT_TXT}"

# ==========================
# DONE
# ==========================
echo "[DONE] Saved:"
echo "  BED: ${OUT_BED}"
echo "  TXT: ${OUT_TXT}"

