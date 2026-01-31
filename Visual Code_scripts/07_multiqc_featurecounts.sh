#!/usr/bin/bash
# Script: 07_multiqc_featureCounts.sh
# Purpose: Feature counts MultiQC
# Usage:  sbatch 07_multiqc_featureCounts.sh

#####################
# SLURM RESOURCES
#####################
#SBATCH --job-name=hisat2_lung
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00            
#SBATCH --output=/data/users/etosun/DE_Analysis_Project/logs/FeatureCount_%j.o
#SBATCH --error=/data/users/etosun/DE_Analysis_Project/logs/FeatureCount_%j.e
#SBATCH --mail-user=elif.tosun@unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

#Load modules
module purge
module load MultiQC

#Directories
PROJECT_DIR=/data/users/etosun/DE_Analysis_Project
OUTDIR="${PROJECT_DIR}/qc/featureCounts"

mkdir -p "${OUTDIR}"

# Find featureCounts summary automatically
FC_SUMMARY=$(find "${PROJECT_DIR}" -name "*counts*.summary" -type f | head -1)
echo "Found: ${FC_SUMMARY}"

if [ -z "${FC_SUMMARY}" ]; then
  echo "No *counts*.summary found! List files:"
  ls -la "${PROJECT_DIR}"/*summary*
  exit 1
fi

# Run MultiQC
multiqc "$(dirname "${FC_SUMMARY}")" \
    -o "${OUTDIR}" \
    --filename multiqc_featurecounts.html \
    --force

echo "Success: ${OUTDIR}/multiqc_featurecounts.html"
