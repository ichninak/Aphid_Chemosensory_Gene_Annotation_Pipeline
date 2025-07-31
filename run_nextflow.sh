#!/bin/bash
#SBATCH --job-name=nextflow_job
#SBATCH --cpus-per-task=50       
#SBATCH --mem=20G                
#SBATCH --time=24:00:00          
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

# Load necessary modules
source /local/miniconda3/etc/profile.d/conda.sh
conda activate aphid-nf

# Variables

echo "Start at $(date)"


# Launch the pipeline
nextflow run main.nf --OR true -resume -profile slurm

echo "end at $(date)"
