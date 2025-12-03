#!/bin/bash
set -e

# ==========================
# Проверка аргументов
# ==========================
if [ $# -ne 2 ]; then
    echo "Usage: $0 <DATA_NAME> <WINDOW>"
    echo "Example: $0 PBMC10k 2000"
    exit 1
fi

DATA_NAME="$1"
WINDOW="$2"

echo "Dataset: $DATA_NAME"
echo "Window:  ${WINDOW}bp"

# ==========================
# Пути
# ==========================
PEAKS="../../Datasets/${DATA_NAME}/input/peaks.bed"
CORE="../../Datasets/${DATA_NAME}/output/celltypist/core_genes_regions_${WINDOW}bp.bed"
OUT="../../Datasets/${DATA_NAME}/output/celltypist/core_genes_regions_${WINDOW}bp_intersected_peaks.bed"

# ==========================
# Проверка файлов
# ==========================
if [ ! -f "$PEAKS" ]; then
    echo "ERROR: peaks file not found: $PEAKS"
    exit 1
fi

if [ ! -f "$CORE" ]; then
    echo "ERROR: core gene regions file not found: $CORE"
    exit 1
fi

# ==========================
# Выполняем пересечение bedtools
# ==========================
echo "Running bedtools intersect..."
echo "  -a: $PEAKS"
echo "  -b: $CORE"
echo "  ->  $OUT"

bedtools intersect \
    -a "$PEAKS" \
    -b "$CORE" \
    -wa -u > "$OUT"

echo "Done."
echo "Saved: $OUT"

