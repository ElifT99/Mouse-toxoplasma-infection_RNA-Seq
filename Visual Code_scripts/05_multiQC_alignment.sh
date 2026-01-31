#!/usr/bin/bash
# Script: 05_multiQC_alignment.sh
# Purpose: QC of hisat2 alignment and summary
# Usage:  sbatch 05_multiQC_alignment.sh

#####################
# SLURM RESOURCES
#####################
#SBATCH --job-name=hisat2_lung
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00            
#SBATCH --output=/data/users/etosun/DE_Analysis_Project/logs/hisat2_%j.o
#SBATCH --error=/data/users/etosun/DE_Analysis_Project/logs/hisat2_%j.e
#SBATCH --mail-user=elif.tosun@unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

#Load modules
module purge
module load MultiQC

#Directories
BASE="/data/users/etosun/DE_Analysis_Project"
ALIGN="${BASE}/alignment"
OUTDIR="${BASE}/qc/alignment/multiqc"

mkdir -p "${OUTDIR}"

#Run MultiQC
multiqc "${ALIGN}" \
    -o "${OUTDIR}" \
    --filename multiqc_alignment.html

echo "MultiQC report written to ${OUTDIR}/multiqc_alignment_report.html"