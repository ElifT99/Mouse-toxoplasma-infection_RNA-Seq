#!/usr/bin/bash
# Script: 06_featureCounts.sh
# Purpose: Feature counts
# Usage:  sbatch 06_featureCounts.sh

#####################
# SLURM RESOURCES
#####################
#SBATCH --job-name=hisat2_lung
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00            
#SBATCH --output=/data/users/etosun/DE_Analysis_Project/logs/FeatureCount_%j.o
#SBATCH --error=/data/users/etosun/DE_Analysis_Project/logs/FeatureCount_%j.e
#SBATCH --mail-user=elif.tosun@unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

#Load module
module purge
module load Subread/2.0.3-GCC-10.3.0
module load SAMtools

#Directories
BASE="/data/users/etosun/DE_Analysis_Project"
ALIGN="${BASE}/alignment"
COUNTS="${BASE}/counts"
REF="${BASE}/ref"
LOGS="${BASE}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${COUNTS}" "${LOGS}"

#Files
GTF="${REF}/gencode.vM30.annotation.gtf"
OUT="${COUNTS}/gene_counts.txt"

#Run FeatureCount
featureCounts \
    -T "${THREADS}" \
    -p \
    -B \
    -C \
    -a "${GTF}" \
    -o "${OUT}" \
    -g gene_id \
    -t exon \
    "${ALIGN}"/*.bam

#-p for paired-end reads