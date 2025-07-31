#!/bin/bash
#SBATCH --job-name=nextflow_job
#SBATCH --cpus-per-task=50       
#SBATCH --mem=100G                
#SBATCH --time=24:00:00          
#SBATCH --output=nextflow_%j.out
#SBATCH --error=nextflow_%j.err

# Load necessary modules
conda activate aphid-nf

# Variables

echo "Start at $(date)"


# Launch the pipeline
nextflow run main.nf --OR true -resume

echo "end at $(date)"
