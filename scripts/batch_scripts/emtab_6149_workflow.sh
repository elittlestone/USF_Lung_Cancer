#!/bin/bash
#SBATCH --partition=snsm_itn19
#SBATCH --qos=openaccess
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=emtab_6149
#SBATCH --output=emtab_6149.out
#SBATCH --error=emtab_6149.err

# Load GPU
module load apps/cuda/12.4

cd ~/myprojects/USF_Lung_Cancer

# Activate pixi
eval "$(pixi shell-hook)"

snakemake -s scripts/snakemake_scripts/e-mtab-6149_workflow.smk -p --cores $SLURM_CPUS_ON_NODE
