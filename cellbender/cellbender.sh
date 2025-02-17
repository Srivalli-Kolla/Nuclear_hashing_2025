#!/bin/bash

# Define base directories
base_input_dir="/home/guest/srivalli/Github/Nuclear_hashing_2025/data"
base_output_dir="/home/guest/srivalli/Github/Nuclear_hashing_2025/data/cellbender_processed_data"

# Define the FPR values
fpr_values=(0.01 0.03 0.05 0.07 0.1)

# Run CellBender with the specified model type
run_cellbender_full() {
    local input_dir=$1
    local fpr=$2

    local input_file="$input_dir/filtered_feature_bc_matrix"
    local output_dir="$base_output_dir/${fpr}_full/"

    mkdir -p "$output_dir"
    
    echo "Running full model for $input_file with FPR $fpr..."

    cellbender remove-background \
        --input "$input_file" \
        --output "$output_dir/${fpr}_after_cb.h5" \
        --fpr "$fpr" \
        --model full \
        --cpu-threads 48 
}

# Loop through FPR values and run CellBender
for fpr in "${fpr_values[@]}"; do
    run_cellbender_full "$base_input_dir" "$fpr"
done