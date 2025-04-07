#!/bin/bash

# Initialize Conda
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate cellbender

cellbender remove-background --model 'full' --fpr '0.1' --input '/home/gruengroup/srivalli/Github/Nuclear_hashing_2025/data/demultiplexed_HTODemux_cellbender' --output '/home/gruengroup/srivalli/Github/Nuclear_hashing_2025/data/cellbender_processed_demultiplexed_HTODemux.h5' 