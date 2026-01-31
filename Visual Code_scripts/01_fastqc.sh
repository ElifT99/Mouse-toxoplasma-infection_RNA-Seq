#!/usr/bin/bash
# Script: 01_fastqc_multiqc.sh
# Purpose: Run FastQC on all lung FASTQ files and summarize with MultiQC
# Usage:  sbatch 01_fastqc_multiqc.sh

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

# Base project directory
PROJECT_DIR=/data/users/etosun/DE_Analysis_Project

#Input and output directories
READS_DIR=${PROJECT_DIR}/raw_data
FASTQC_DIR=${PROJECT_DIR}/qc/fastqc
MULTIQC_DIR=${PROJECT_DIR}/qc

#Load modules
module load FastQC
module load MultiQC

# 1. Run FASTQC
fastqc \
  -t ${SLURM_CPUS_PER_TASK} \
  -o "${FASTQC_DIR}" \
  "${READS_DIR}"/*.fastq.gz

# -t ${SLURM_CPUS_PER_TASK} determines the CPU used for the job

# 2. Run MULTIQC
multiqc "${FASTQC_DIR}" -o "${MULTIQC_DIR}"
