#!/bin/bash
#SBATCH --partition=comet
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G


sbatch_command=$(sbatch --help)
result=$?

dataset=$1
method=$2
DEFAULT_FOLDER=../../Datasets/$dataset/output
if [ $result -eq 0 ] 
then
    echo "Sbatch found. Running pipeline on cluster"
    sbatch run_scopen.sh $DEFAULT_FOLDER/cicero 10X $DEFAULT_FOLDER/scopen_cicero
#    sbatch run_scopen.sh $DEFAULT_FOLDER/threshold 10X $DEFAULT_FOLDER/scopen_threshold
else
    echo "Sbatch not found. Please, notify your users about hard tasks. In 10 seconds it will start"
    sleep 1
    bash run_scopen.sh $DEFAULT_FOLDER/$method 10X $DEFAULT_FOLDER/${method}_scopen
#    bash run_scopen.sh $DEFAULT_FOLDER/threshold 10X $DEFAULT_FOLDER/scopen_threshold
fi
