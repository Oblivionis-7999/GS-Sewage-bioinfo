#!/bin/bash

# ==============================================================================
# 0. DYNAMIC BATCH SETUP & AUTOMATED LOGGING
# ==============================================================================
BATCH_ID="$1"

if [ -z "$BATCH_ID" ]; then
    echo "Error: You must provide a batch date!"
    echo "Usage: ./52week.sh 261001"
    exit 1
fi

# Define Master Variables
BATCH_DIR="/Users/bartoncheng/projects/52-week-wbe/batch${BATCH_ID}"
REF_GENOME="/Users/bartoncheng/projects/SARS-CoV-2_reference_genome/reference-NC_045512.fasta"
BED_FILE="/Users/bartoncheng/projects/SARS-CoV-2_reference_genome/artic_532.bed"

# Sub-directories
RAW_DIR="${BATCH_DIR}/raw"
CLEAN_DIR="${BATCH_DIR}/fastp"
REPORT_DIR="${BATCH_DIR}/fastp_reports"
ALIGNED_DIR="${BATCH_DIR}/aligned"
TRIMMED_DIR="${BATCH_DIR}/trimmed"
VARIANTS_DIR="${BATCH_DIR}/variants_files"
DEPTH_DIR="${BATCH_DIR}/depth_files"
DEMIX_DIR="${BATCH_DIR}/demix_files"

# Create all folders upfront
mkdir -p "$CLEAN_DIR" "$REPORT_DIR" "$ALIGNED_DIR" "$TRIMMED_DIR" "$VARIANTS_DIR" "$DEPTH_DIR" "$DEMIX_DIR"

# Export variables for xargs
export BATCH_DIR REF_GENOME BED_FILE RAW_DIR CLEAN_DIR REPORT_DIR ALIGNED_DIR TRIMMED_DIR VARIANTS_DIR DEPTH_DIR DEMIX_DIR

# Start automated logging
LOG_FILE="${BATCH_DIR}/run_log_${BATCH_ID}.txt"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================================================"
echo "Starting Pipeline for Batch: $BATCH_ID"
echo "========================================================================"

# Initialize conda for the bash script
eval "$(conda shell.bash hook)"


# ==============================================================================
# STAGE 1: fastp
# ==============================================================================
echo "Starting Stage 1: fastp"
conda activate fastp

ls "$RAW_DIR"/*_R1.fastq.gz | xargs -n 1 -P 4 sh -c '
    r1="$1"
    name=$(basename "$r1" _R1.fastq.gz)
    r2="${RAW_DIR}/${name}_R2.fastq.gz"
    
    echo "Running fastp QC & Trimming for: $name"
    fastp -i "$r1" -I "$r2" \
          -o "${CLEAN_DIR}/${name}_R1.fastq.gz" -O "${CLEAN_DIR}/${name}_R2.fastq.gz" \
          -h "${REPORT_DIR}/${name}_fastp.html" -j "${REPORT_DIR}/${name}_fastp.json" \
          --detect_adapter_for_pe -q 20 -l 50 -w 2 >/dev/null 2>&1
    
    echo "Finished $name!"
' _

conda deactivate
echo "Stage 1 Complete!"
echo "------------------------------------------------------------------------"


# ==============================================================================
# STAGE 2: BWA-Mem (reference mapping)
# ==============================================================================
echo "Starting Stage 2: BWA Alignment"
conda activate bwa

for r1_file in "$CLEAN_DIR"/*_R1.fastq.gz; do
    name=$(basename "$r1_file" _R1.fastq.gz)
    r2_file="${CLEAN_DIR}/${name}_R2.fastq.gz"
    
    echo "Aligning Sample: $name"
    bwa mem -t 8 "$REF_GENOME" "$r1_file" "$r2_file" > "${ALIGNED_DIR}/${name}_aligned.sam" 2>/dev/null
done

conda deactivate
echo "Stage 2 Complete!"
echo "------------------------------------------------------------------------"


# ==============================================================================
# STAGE 3: Samtools (Sorting & Compressing)
# ==============================================================================
echo "Starting Stage 3: Samtools Sort"
# UPDATED: Now uses your new combined 'ivar' environment
conda activate ivar

for sam_file in "$ALIGNED_DIR"/*_aligned.sam; do
    name=$(basename "$sam_file" _aligned.sam)
    
    echo "Sorting & Compressing Sample: $name"
    samtools sort -@ 8 -o "${ALIGNED_DIR}/${name}_aligned_sorted.bam" "$sam_file"
    rm "$sam_file"
done

conda deactivate
echo "Stage 3 Complete!"
echo "------------------------------------------------------------------------"


# ==============================================================================
# STAGE 4: iVar (Primer Trimming)
# ==============================================================================
echo "Starting Stage 4: iVar Trim"
# Uses the exact same 'ivar' environment as Stage 3
conda activate ivar

ls "$ALIGNED_DIR"/*_aligned_sorted.bam | xargs -n 1 -P 4 sh -c '
    bam="$1"
    name=$(basename "$bam" _aligned_sorted.bam)
    
    echo "Trimming primers for: $name"
    
    # Step A: Run ivar trim
    ivar trim -i "$bam" -b "$BED_FILE" -p "${TRIMMED_DIR}/${name}_trimmed" -e -q 20 >/dev/null 2>&1
    
    # Step B: Sort and index the newly trimmed file
    samtools sort -@ 2 -o "${TRIMMED_DIR}/${name}_trimmed_sorted.bam" "${TRIMMED_DIR}/${name}_trimmed.bam"
    samtools index "${TRIMMED_DIR}/${name}_trimmed_sorted.bam"
    
    # Step C: Clean up the intermediate unsorted file
    rm "${TRIMMED_DIR}/${name}_trimmed.bam"
    
    echo "Finished trimming $name!"
' _

conda deactivate
echo "Stage 4 Complete!"
echo "------------------------------------------------------------------------"


# ==============================================================================
# STAGE 5: Freyja (Variants and Depths)
# ==============================================================================
echo "Starting Stage 5: Freyja"
conda activate freyja

echo "Updating Freyja databases..."
freyja update

ls "$TRIMMED_DIR"/*_trimmed_sorted.bam | xargs -n 1 -P 8 sh -c '
    bam="$1"
    name=$(basename "$bam" _trimmed_sorted.bam)

    echo "Processing variants & demix for: $name"

    freyja variants "$bam" \
        --variants "${VARIANTS_DIR}/${name}.variants.tsv" \
        --depths "${DEPTH_DIR}/${name}.depth.tsv" \
        --ref "$REF_GENOME"

    freyja demix "${VARIANTS_DIR}/${name}.variants.tsv" "${DEPTH_DIR}/${name}.depth.tsv" \
        --output "${DEMIX_DIR}/${name}.demix.tsv" --confirmedonly

    echo "Finished $name!"
' _

conda deactivate

echo "========================================="
echo "Aggregating all demix results..."
echo "========================================="
conda activate freyja

freyja aggregate "$DEMIX_DIR/" --output "${BATCH_DIR}/Master_Freyja_Summary.tsv"

conda deactivate

echo "========================================================================"
echo "PIPELINE FULLY COMPLETE!"
echo "Master summary saved to: ${BATCH_DIR}/Master_Freyja_Summary.tsv"
echo "Run log saved to: $LOG_FILE"
echo "========================================================================"
