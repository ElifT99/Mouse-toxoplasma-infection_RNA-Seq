#!/usr/bin/bash
# Script: 02_multiqc.sh
# Purpose: Run FastQC on all lung FASTQ files and summarize with MultiQC
# Usage:  sbatch 02_multiqc.sh

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

# Base project directory, input and output directory
PROJECT_DIR=/data/users/etosun/DE_Analysis_Project
FASTQC_DIR=${PROJECT_DIR}/qc/fastqc
ORIG_DIR=${PROJECT_DIR}/qc/fastqc_original
MULTIQC_DIR=${PROJECT_DIR}/qc

#Move the fastqc report of the re-run sample SRR7821937_1.fastq.gz in the overall fastqc folder
cd /data/users/etosun/DE_Analysis_Project
mv qc/fastqc_original/* qc/fastqc/

#Load MultiQC module
module load MultiQC

#Perform MultiQC
multiqc "${FASTQC_DIR}" -o "${MULTIQC_DIR}"