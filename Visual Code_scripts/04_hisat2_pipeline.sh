#!/usr/bin/bash
# Script: 04_hisat2_pipeline.sh
# Purpose: Hisat2 pipeline
# Usage:  sbatch 04_hisat2_pipeline.sh

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

#Load module
module purge
module load Biopython/1.79-foss-2021a
module load HISAT2/2.2.1-gompi-2021a
module load SAMtools

which python hisat2_extract_splice_sites.py hisat2-build

#Directories
BASE="/data/users/etosun/DE_Analysis_Project"
REF="${BASE}/ref"
FASTQ="${BASE}/raw_data"
ALIGN="${BASE}/alignment"
LOGS="${BASE}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Mouse reference genome sequence (FASTA file) and annotation file
FA="${REF}/mm39.fa"
GTF="${REF}/gencode.vM30.annotation.gtf"

# Index location
IDX="${BASE}/index/hisat2/mm39_index/mm39_index"

mkdir -p "${ALIGN}" "${LOGS}"

# Splice sites
SPLICE="${BASE}/index/hisat2/mm39_index/splicesites.txt"
[[ -s "${SPLICE}" ]] || hisat2_extract_splice_sites.py "${GTF}" > "${SPLICE}"

# Samples list
SAMPLES="${BASE}/samples.txt" #sample list creation
[[ -s "${SAMPLES}" ]] || {
    ls "${FASTQ}"/*_1.fastq.gz | sed 's|.*/||; s/_1.*//' > "${SAMPLES}"
    echo "Created samples.txt: $(wc -l < "${SAMPLES}") samples"
}

echo "Mapping $(wc -l < "${SAMPLES}") samples..." #loops through each line in samples.txt

while read SAMPLE; do
    [[ -z "${SAMPLE}" ]] && continue  # skip empty lines

    R1="${FASTQ}/${SAMPLE}_1.fastq.gz" #because of paired-end RNA-Seq, files are named R1 and R2
    R2="${FASTQ}/${SAMPLE}_2.fastq.gz"

    [[ ! -f "${R1}" || ! -f "${R2}" ]] && {
        echo "Missing FASTQ files for ${SAMPLE}, skipping" >> "${LOGS}/missing_fastq.log"
        continue
    }
    #FASTQ integrity check
    if ! gzip -t "${R1}" || ! gzip -t "${R2}"; then
        echo "Corrupted FASTQ for ${SAMPLE}, skipping" >> "${LOGS}/corrupted_fastq.log"
        continue
    fi

    BAM="${ALIGN}/${SAMPLE}.bam"
    LOG="${ALIGN}/${SAMPLE}.hisat2.log"

    #skip samples that already aligned successfully
    if [[ -s "${BAM}" ]] && samtools quickcheck "${BAM}" 2>/dev/null; then
        echo "[$(date)] ${SAMPLE}: BAM already exists, skipping"
        continue
    fi
    #clean up leftover temp files from previous crash
    rm -f "${ALIGN}/${SAMPLE}.bam.tmp."*

    echo "[$(date)] Mapping ${SAMPLE}"

    #hisat2 aligment
    hisat2 -p "${THREADS}" --dta --known-splicesite-infile "${SPLICE}" \
        -x "${IDX}" -1 "${R1}" -2 "${R2}" -S /dev/stdout 2> "${LOG}" | \
    samtools sort -@ "${THREADS}" -o "${BAM}" -

    #Post-processing, checking of alignemnt quality and indexing
    samtools index "${BAM}"
    samtools flagstat "${BAM}" > "${ALIGN}/${SAMPLE}.flagstat.txt"
done < "${SAMPLES}"