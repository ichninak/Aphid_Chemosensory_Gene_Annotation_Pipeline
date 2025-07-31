#!/bin/bash
#SBATCH --job-name=nextflow_job
#SBATCH --cpus-per-task=4        # Juste pour Nextflow principal
#SBATCH --mem=8G                 # Juste pour Nextflow principal                
#SBATCH --time=24:00:00          
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

# Load necessary modules
source /local/miniconda3/etc/profile.d/conda.sh
conda activate aphid-nf

# Variables

echo "Start at $(date)"


# Launch the pipeline
nextflow run main.nf --OR true -profile slurm

echo "end at $(date)"
