#!/usr/bin/env bash
#SBATCH --job-name=jules_driver
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=co_moilab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=/global/scratch/projects/fc_moilab/julesperez/snakemake_pipeline/panNrefNkmer/workflow/driver_%j.out

# Snakemake driver on a compute node (avoids stale Lustre handles on login nodes).
# Submit from workflow/:   sbatch scripts/run_driver.sbatch
# Resubmit after it ends; Snakemake resumes from what is finished.
# Stop it with:            scancel <jobid>   then run the --unlock line below.

set -uo pipefail

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake

cd /global/scratch/projects/fc_moilab/julesperez/snakemake_pipeline/panNrefNkmer/workflow

TARGETS="jules_only jules_roh_only jules_msmc2_only"

snakemake --sdm conda --profile profiles/default $TARGETS --unlock

snakemake --sdm conda --profile profiles/default $TARGETS \
    --keep-going --rerun-incomplete --rerun-triggers mtime
