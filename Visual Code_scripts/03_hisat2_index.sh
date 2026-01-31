#!/usr/bin/bash
# Script: 03_hisat2_index.sh
# Purpose: Hisat2 index generation
# Usage:  sbatch 03_hisat2_index.sh

#####################
# SLURM RESOURCES
#####################
#SBATCH --job-name=hisat2_lung
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=04:00:00            
#SBATCH --output=/data/users/etosun/DE_Analysis_Project/logs/hisat2_%j.o
#SBATCH --error=/data/users/etosun/DE_Analysis_Project/logs/hisat2_%j.e
#SBATCH --mail-user=elif.tosun@unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

#Load module
module purge
module load Biopython/1.79-foss-2021a
module load HISAT2/2.2.1-gompi-2021a
module load SAMtools

#Project folder structure
BASE="/data/users/etosun/DE_Analysis_Project"
REF="${BASE}/ref"
IDX="${BASE}/index/hisat2/mm39_index"

mkdir -p "${IDX}" "${BASE}/logs"
cd "${IDX}"

#hisat2 version check
which hisat2-build
hisat2-build --version

#Splice sites, reads GTF --> extract exon-exon junction --> splicesites.txt
[[ -s splicesites.txt ]] || \
hisat2_extract_splice_sites.py "${REF}/gencode.vM30.annotation.gtf" > splicesites.txt

#Build Index
hisat2-build -p 8 \
    --ss splicesites.txt \
    "${REF}/mm39.fa" \
    mm39_index

echo "Index files:"
ls -lh mm39_index.*.ht2