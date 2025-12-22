#!/bin/bash


DATA_NAME="$1"
INPUT_DIR="../../Datasets/$DATA_NAME/"
BED_FILE1="$INPUT_DIR/input/peaks.bed"
OUT_FILE1="$INPUT_DIR/input/peaks.txt"

#BED_FILE2="$INPUT_DIR/output/raw/peaks.bed"
#OUT_FILE2="$INPUT_DIR/output/raw/peaks.txt"

if [ ! -f "$BED_FILE1" ]; then
    echo "Ошибка: файл $BED_FILE1 не найден."
    exit 1
fi

#if [ ! -f "$BED_FILE2" ]; then
#    echo "Ошибка: файл $BED_FILE2 не найден."
#    exit 1
#fi


#echo "Конвертирую $BED_FILE1 и $BED_FILE2 в $OUT_FILE1 и $OUT_FILE2 ..."
echo "Конвертирую $BED_FILE1 в $OUT_FILE1..."

awk 'BEGIN{OFS="_"} {print $1, $2, $3}' "$BED_FILE1" > "$OUT_FILE1"

#awk 'BEGIN{OFS="_"} {print $1, $2, $3}' "$BED_FILE2" > "$OUT_FILE2"


echo "Готово."

