#!/bin/bash

folder=$1
output_folder=$2
labels=$3

mkdir -p "$output_folder"

MAX_JOBS=16

running=0

for file in $(ls "$folder"); do
    echo "Запускаю $file"

    input_folder="$folder/$file"
    output_metrics_folder="$output_folder/$file"
    mkdir -p "$output_metrics_folder"

    bash run_clustering.sh \
        "$input_folder" \
        "$labels" \
        "$output_metrics_folder/metrics.json" \
        dense matrix.csv &

    running=$((running+1))

    if [[ $running -ge $MAX_JOBS ]]; then
        wait
        running=0
    fi
done

wait
echo "Закончил работу."

