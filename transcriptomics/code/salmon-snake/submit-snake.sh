#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=180GB
#SBATCH --time=1000

jobname=$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2)
outdir=$(cat config.yaml | grep outputDIR | cut -d ":" -f 2 | cut -d " " -f 2)

SNAKEFILE=$1
JOBS=$2 # should default to 500
RERUN_INCOMPLETE=$3
DRYRUN=$4

JOBS=100

snakemake  \
    -s $SNAKEFILE $RERUN_INCOMPLETE $DRYRUN --jobs $JOBS --use-conda \
    --cluster "sbatch --time=1000 --job-name=salmon-phys --mem=180GB"
