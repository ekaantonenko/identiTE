#!/bin/sh
#
#
# Script: Precise annotation of transposable elements across a genomes dataset
# Description: This script performs the following steps:
#   1. 
#   2.
# Usage: 
# Author: Katia Antonenko (ekaantonenko.github.io)
# Date: 16/04/2026

#set -euo pipefail


set -e

CONFIG_FILE="config_default.yaml"
run_step1=false
run_step2=false
run_step3=false
run_step4=false

# Parse arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    --config) CONFIG_FILE="$2"; shift ;;
    --step1) run_step1=true ;;
    --step2) run_step2=true ;;
    --step3) run_step3=true ;;
    --step4) run_step4=true ;;
    --all)
      run_step1=true
      run_step2=true
      run_step3=true
      run_step4=true
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
  shift
done

if [ ! -f "$CONFIG_FILE" ]; then
  echo "Config file not found: $CONFIG_FILE"
  exit 1
fi

if ! $run_step1 && ! $run_step2 && ! $run_step3 && ! $run_step4; then
  echo "No steps specified, running all"
  run_step1=true
  run_step2=true
  run_step3=true
  run_step4=true
fi

########################################################################################################################
# LOAD CONFIG VARIABLES
########################################################################################################################


batch_size=$(yq '.batch_size' "$CONFIG_FILE")
batch_first=$(yq '.batch_first' "$CONFIG_FILE")
batch_last=$(yq '.batch_last' "$CONFIG_FILE")

library_path=$(yq -r '.library_path' "$CONFIG_FILE")
library_format=$(yq -r '.library_format' "$CONFIG_FILE")
ids_path=$(yq -r '.ids_path' "$CONFIG_FILE")
gene_path=$(yq -r '.gene_path' "$CONFIG_FILE")
gene_ref_path=$(yq -r '.gene_ref_path' "$CONFIG_FILE")
seq_path=$(yq -r '.seq_path' "$CONFIG_FILE")``
seq_ref_path=$(yq -r '.seq_ref_path' "$CONFIG_FILE")
peri_path=$(yq -r '.peri_path' "$CONFIG_FILE")
out_path=$(yq -r '.out_path' "$CONFIG_FILE")

gene_suffix=$(yq -r '.gene_suffix' "$CONFIG_FILE")
gene_format=$(yq -r '.gene_format' "$CONFIG_FILE")
seq_suffix=$(yq -r '.seq_suffix' "$CONFIG_FILE")

returnTE=$(yq '.returnTE' "$CONFIG_FILE")
returnFlanks=$(yq '.returnFlanks' "$CONFIG_FILE")

meth_path=$(yq -r '.meth_path' "$CONFIG_FILE")
meth_format=$(yq -r '.meth_format' "$CONFIG_FILE")
meth_suffix=$(yq -r '.meth_suffix' "$CONFIG_FILE")
threads=$(yq '.threads' "$CONFIG_FILE")


mkdir -p "$out_path"
JOB_ID=${SLURM_JOB_ID:-manual}

exec > >(tee $out_path/annotate.${JOB_ID}.log) 2>&1
echo "Using config: $CONFIG_FILE"
cp "$CONFIG_FILE" "$out_path/config_used.yaml"
echo "Running steps:"
$run_step1 && echo "  - Step 1"
$run_step2 && echo "  - Step 2"
$run_step3 && echo "  - Step 3"
$run_step4 && echo "  - Step 4"

mkdir -p $out_path

########################################################################################################################
# Step 1: blast TE sequences to the genomes sequences
########################################################################################################################

if [ "$run_step1" = true ]; then
        echo "Step 1: blast TE sequences to the genomes sequences"
        echo "Blasting started : $(date)"
        mkdir -p ${out_path}/blast
        mkdir -p ${out_path}/blastedSummary

        python3 -u ./modules/1_Blast_TEtoGenomes_flanks_and_TE.py \
                --library_path=$library_path \
                --library_format=$library_format \
                --ids_path=$ids_path \
                --gene_path=$gene_path \
                --gene_ref_path=$gene_ref_path \
                --gene_suffix=$gene_suffix \
                --gene_format=$gene_format \
                --seq_path=$seq_path \
                --seq_ref_path=$seq_ref_path \
                --seq_suffix=$seq_suffix \
                --peri_path=$peri_path \
                --out_path=$out_path  \
                --batch_first=$batch_first \
                --batch_last=$batch_last \
                --batch_size=$batch_size 
                # --n_threads=$n_threads

        rm -r ${out_path}/blast
        echo "Blasting completed : $(date)"
fi

########################################################################################################################
# Step 2: splitting the blasted TEs and flanking regions into individual genomes
########################################################################################################################

blasted_path=$out_path/blastedSummary/blastedSummary.csv
returnTE=True
returnFlanks=True

if [ "$run_step2" = true ]; then
        echo "Step 2: combining blasted TEs"
        echo "Combining started : $(date)"

        python3 -u ./modules/2_CombiningBlastedTEs.py \
                --blasted_path=$blasted_path \
                --ids_path=$ids_path \
                --out_path=$out_path \
                --returnTE=$returnTE \
                --returnFlanks=$returnFlanks 
                # --n_threads=$n_threads


        echo "Combining completed : $(date)"
fi

########################################################################################################################
# Step 3: collecting the methylation levels + nearest genes annotation for the Transposable Elements
########################################################################################################################

te_path="$out_path/TEperGenome"

if [ "$run_step3" = true ]; then

        echo "Step 3: adding genes and methylation to TEs"
        echo "Started : $(date)"

        bash ./modules/3_addGenesMeth_toTE.sh \
        --ids    "$ids_path" \
        --te "$te_path" \
        --gene "$gene_path" \
        --gene_format "$gene_format" \
        --gene_suffix "$gene_suffix" \
        --meth   "$meth_path" \
        --meth_format "$meth_format" \
        --meth_suffix "$meth_suffix" \
        --out    "$out_path" \
        --threads "$threads"
        
        echo "Completed : $(date)"
fi

########################################################################################################################
# Step 4 (optional): collecting the methylation levels of the flanks
########################################################################################################################

flanks_path="$out_path/FlanksPerGenome"


if [ "$run_step4" = true ]; then

        echo "Step 4: adding methylation to flanks"
        echo "Methylation (flanks) started : $(date)"

        bash ./modules/4_addMeth_to_flanks.sh \
        --ids    "$ids_path" \
        --flanks "$flanks_path" \
        --meth   "$meth_path" \
        --meth_format "$meth_format" \
        --meth_suffix "$meth_suffix" \
        --out    "$out_path" \
        --threads "$threads"
        
        # Source instead of bash so set -euo pipefail and exit codes propagate
        #bash "$(dirname "$0")/modules/_addMeth_to_flanks.sh"
        
        echo "Methylation (flanks) completed : $(date)"
fi