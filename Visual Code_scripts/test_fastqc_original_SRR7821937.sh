#!/usr/bin/bash
# Script: test_fastqc_original_SRR7821937.sh
# Purpose: To solve problem with the raw data of SRR7821937_1
# Usage:  sbatch test_fastqc_original_SRR7821937.sh

#####################
# SLURM RESOURCES
#####################
#SBATCH --job-name=qc_lung
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=02:00:00            
#SBATCH --output=/data/users/etosun/DE_Analysis_Project/logs/qc_%j.o
#SBATCH --error=/data/users/etosun/DE_Analysis_Project/logs/qc_%j.e
#SBATCH --mail-user=elif.tosun@unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

#While performing fastqc on my copied data I got an error message for the run SRR7821937_1.fastq.gz, to see if the error is caused by an copy error I will run fastqc on the original data
module load FastQC

ORIG_DIR=/data/courses/rnaseq_course/toxoplasma_de/reads_Lung
OUT_DIR=/data/users/etosun/DE_Analysis_Project/qc/fastqc_original

mkdir -p "${OUT_DIR}"

fastqc -t ${SLURM_CPUS_PER_TASK} \
  -o "${OUT_DIR}" \
  "${ORIG_DIR}/SRR7821937_1.fastq.gz"