#!/usr/bin/env bash
set -euo pipefail

input_file="${1:?missing input_file}"
peaks_prefix="${2:?missing peaks_prefix}"
organism="${3:-human}"

if [[ ! -s "${input_file}" ]]; then
  echo "[ERROR] input_file missing/empty: ${input_file}" >&2
  exit 2
fi

# Определяем диапазон хромосом
max_chr=22
if [[ "${organism}" == "mouse" ]]; then
  max_chr=19
fi

out_dir="$(dirname "${peaks_prefix}")"
mkdir -p "${out_dir}"

# Заранее создаём пустые файлы (чтобы downstream не ломался, даже если нет строк)
for i in $(seq 1 "${max_chr}"); do
  : > "${peaks_prefix}_chr${i}.tsv"
done
: > "${peaks_prefix}_chrX.tsv"
: > "${peaks_prefix}_chrY.tsv"

# Один проход: по первому полю (peak_id вида chr1_10109_10357) выбираем chr
awk -v pref="${peaks_prefix}" -v max_chr="${max_chr}" '
BEGIN { FS="\t"; OFS="\t" }
{
  # peak_id в первом поле: chrN_start_end
  # извлекаем "chrN" как подстроку до первого "_"
  split($1, a, "_")
  chr = a[1]          # "chr1", "chrX", "chrY", ...
  sub(/^chr/, "", chr) # chr -> "1"/"X"/"Y"

  # фильтруем хромосомы: 1..max_chr, X, Y
  if (chr == "X" || chr == "Y") {
    print $0 >> (pref "_chr" chr ".tsv")
  } else if (chr ~ /^[0-9]+$/) {
    c = chr + 0
    if (c >= 1 && c <= max_chr) {
      print $0 >> (pref "_chr" c ".tsv")
    }
  }
}
' "${input_file}"

echo "[INFO] Split done: ${peaks_prefix}_chr*.tsv (1..${max_chr}, X, Y)"

