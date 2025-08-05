#!/bin/bash
#SBATCH --job-name=snakemake_aphid
#SBATCH --cpus-per-task=2        # CPUs pour le processus Snakemake principal
#SBATCH --mem=4G                 # Mémoire pour Snakemake principal
#SBATCH --time=48:00:00          # Temps maximum pour tout le pipeline
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err

# Chargement de l'environnement
source /local/miniconda3/etc/profile.d/conda.sh
conda activate VOTRE_NOUVEL_ENVIRONNEMENT

# Variables de configuration
WORKFLOW_TYPE=${1:-"OR"}  # OR, GR, ou ALL par défaut

echo "Démarrage du pipeline Snakemake à $(date)"
echo "Type de workflow: $WORKFLOW_TYPE"

# Fonction pour exécuter Snakemake avec SLURM
run_snakemake() {
    local config_flags="$1"
    
    snakemake \
        --jobs 100 \
        --cluster-config cluster.yaml \
        --cluster "sbatch --job-name={cluster.name} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --time={cluster.time} --output={cluster.output} --error={cluster.error}" \
        --latency-wait 60 \
        --rerun-incomplete \
        --keep-going \
        $config_flags \
        --printshellcmds \
        --reason
}

# Exécution selon le type de workflow
case $WORKFLOW_TYPE in
    "OR")
        echo "Exécution du workflow OR uniquement"
        run_snakemake "--config OR=true"
        ;;
    "GR")
        echo "Exécution du workflow GR uniquement"
        run_snakemake "--config GR=true"
        ;;
    "ALL")
        echo "Exécution des workflows OR et GR"
        run_snakemake "--config OR=true GR=true"
        ;;
    *)
        echo "Type de workflow invalide: $WORKFLOW_TYPE"
        echo "Utilisation: $0 [OR|GR|ALL]"
        exit 1
        ;;
esac

echo "Fin du pipeline Snakemake à $(date)"
