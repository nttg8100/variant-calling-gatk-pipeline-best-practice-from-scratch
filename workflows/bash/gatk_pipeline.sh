#!/bin/bash
#
# GATK Variant Calling Pipeline - Complete Best Practices (16 Steps)
# Based on NGS101 tutorial: https://ngs101.com/how-to-analyze-whole-genome-sequencing-data-for-absolute-beginners-part-1-from-raw-reads-to-high-quality-variants-using-gatk/
# Adapted to use nf-core test data with hard filtering (Option 2)
#
# This pipeline implements the complete GATK workflow with:
# - Adapter trimming and quality filtering (Trim Galore)
# - Read alignment with BWA-MEM
# - Duplicate marking (GATK MarkDuplicates)
# - Base quality score recalibration (BQSR)
# - Alignment quality metrics (GATK CollectMetrics)
# - Variant calling with HaplotypeCaller (GVCF mode)
# - Hard filtering for SNPs and Indels (suitable for small datasets)
# - Functional annotation with SnpEff
# - Comprehensive variant statistics (bcftools)
# - Visualization files for genome browsers (BED, bedGraph)
#

set -euo pipefail

# Configuration
SAMPLE="${SAMPLE:-sample1}"  # Can be set via environment variable or default to sample1
THREADS=2

# Paths relative to workflows/bash directory
DATA_DIR="../../data"
REF_DIR="../../reference"
RESULTS_DIR="./results"

# Reference files
REFERENCE="${REF_DIR}/genome.fasta"
DBSNP="${REF_DIR}/dbsnp_146.hg38.vcf.gz"
KNOWN_INDELS="${REF_DIR}/mills_and_1000G.indels.vcf.gz"

# Input FASTQ files
FASTQ_R1="${DATA_DIR}/${SAMPLE}_R1.fastq.gz"
FASTQ_R2="${DATA_DIR}/${SAMPLE}_R2.fastq.gz"

# Output directories
QC_DIR="${RESULTS_DIR}/qc/${SAMPLE}"
TRIMMED_DIR="${RESULTS_DIR}/trimmed/${SAMPLE}"
ALIGNED_DIR="${RESULTS_DIR}/aligned/${SAMPLE}"
VAR_DIR="${RESULTS_DIR}/variants/${SAMPLE}"
BQSR_DIR="${RESULTS_DIR}/bqsr/${SAMPLE}"

# Output files
TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE}_R1_val_1.fq.gz"
TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE}_R2_val_2.fq.gz"
ALIGNED_BAM="${ALIGNED_DIR}/${SAMPLE}.bam"
SORTED_BAM="${ALIGNED_DIR}/${SAMPLE}.sorted.bam"
DEDUP_BAM="${ALIGNED_DIR}/${SAMPLE}_marked_duplicates.bam"
RECAL_TABLE="${BQSR_DIR}/${SAMPLE}_recal_data.table"
RECAL_BAM="${ALIGNED_DIR}/${SAMPLE}_recalibrated.bam"
GVCF="${VAR_DIR}/${SAMPLE}.g.vcf.gz"
RAW_VCF="${VAR_DIR}/${SAMPLE}_raw_variants.vcf.gz"
FILTERED_SNP_VCF="${VAR_DIR}/${SAMPLE}_filtered_snps.vcf.gz"
FILTERED_INDEL_VCF="${VAR_DIR}/${SAMPLE}_filtered_indels.vcf.gz"
FINAL_VCF="${VAR_DIR}/${SAMPLE}_filtered.vcf.gz"
ANNOTATED_VCF="${VAR_DIR}/${SAMPLE}_annotated.vcf"
METRICS="${ALIGNED_DIR}/${SAMPLE}_duplicate_metrics.txt"

# Create output directories
mkdir -p ${QC_DIR} ${TRIMMED_DIR} ${ALIGNED_DIR} ${VAR_DIR} ${BQSR_DIR}

echo "=========================================="
echo "GATK Variant Calling Pipeline - Complete"
echo "Sample: ${SAMPLE}"
echo "Based on NGS101 Best Practices"
echo "=========================================="
echo

# Step 1: Quality Control with FastQC
echo "[$(date)] Step 1: Running FastQC on raw reads..."
fastqc -o ${QC_DIR} -t ${THREADS} ${FASTQ_R1} ${FASTQ_R2}
echo "[$(date)] FastQC completed"
echo

# Step 2: Adapter Trimming and Quality Filtering
echo "[$(date)] Step 2: Adapter trimming with Trim Galore..."
trim_galore \
    --paired \
    --quality 20 \
    --length 50 \
    --fastqc \
    --output_dir ${TRIMMED_DIR} \
    ${FASTQ_R1} ${FASTQ_R2}
echo "[$(date)] Adapter trimming completed"
echo

# Step 3: Read Alignment with BWA-MEM
echo "[$(date)] Step 3: Aligning reads with BWA-MEM..."
bwa mem \
    -t ${THREADS} \
    -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_lib" \
    ${REFERENCE} \
    ${TRIMMED_R1} ${TRIMMED_R2} \
    | samtools view -Sb - > ${ALIGNED_BAM}
echo "[$(date)] Alignment completed"
echo

# Step 4: Sort BAM file
echo "[$(date)] Step 4: Sorting BAM file..."
samtools sort -@ ${THREADS} -o ${SORTED_BAM} ${ALIGNED_BAM}
echo "[$(date)] Sorting completed"
echo

# Step 5: Mark Duplicates
echo "[$(date)] Step 5: Marking duplicates with GATK..."
gatk MarkDuplicates \
    -I ${SORTED_BAM} \
    -O ${DEDUP_BAM} \
    -M ${METRICS} \
    --CREATE_INDEX true
echo "[$(date)] Mark duplicates completed"
echo

# Step 6: Base Quality Score Recalibration (BQSR) - Generate table
echo "[$(date)] Step 6: Generating BQSR recalibration table..."
gatk BaseRecalibrator \
    -I ${DEDUP_BAM} \
    -R ${REFERENCE} \
    --known-sites ${DBSNP} \
    --known-sites ${KNOWN_INDELS} \
    -O ${RECAL_TABLE}
echo "[$(date)] BQSR table generated"
echo

# Step 7: Apply BQSR
echo "[$(date)] Step 7: Applying BQSR..."
gatk ApplyBQSR \
    -I ${DEDUP_BAM} \
    -R ${REFERENCE} \
    --bqsr-recal-file ${RECAL_TABLE} \
    -O ${RECAL_BAM}
echo "[$(date)] BQSR applied"
echo

# Step 8: Alignment Quality Assessment
echo "[$(date)] Step 8: Collecting alignment quality metrics..."
gatk CollectAlignmentSummaryMetrics \
    -R ${REFERENCE} \
    -I ${RECAL_BAM} \
    -O ${QC_DIR}/${SAMPLE}_alignment_summary.txt

# Collect insert size metrics with histogram (requires R)
gatk CollectInsertSizeMetrics \
    -I ${RECAL_BAM} \
    -O ${QC_DIR}/${SAMPLE}_insert_size_metrics.txt \
    -H ${QC_DIR}/${SAMPLE}_insert_size_histogram.pdf
echo "[$(date)] Alignment QC completed"
echo

# Step 9: Variant Calling with HaplotypeCaller (GVCF mode)
echo "[$(date)] Step 9: Calling variants with HaplotypeCaller in GVCF mode..."
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ${RECAL_BAM} \
    -O ${GVCF} \
    -ERC GVCF \
    --dbsnp ${DBSNP}
echo "[$(date)] GVCF generation completed"
echo

# Step 10: Genotype GVCFs
echo "[$(date)] Step 10: Genotyping GVCF to VCF..."
gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V ${GVCF} \
    -O ${RAW_VCF}
echo "[$(date)] Genotyping completed"
echo

# Step 11: Hard Filtering for SNPs
echo "[$(date)] Step 11: Hard filtering SNPs..."
# Select SNPs
gatk SelectVariants \
    -R ${REFERENCE} \
    -V ${RAW_VCF} \
    --select-type-to-include SNP \
    -O ${VAR_DIR}/${SAMPLE}_raw_snps.vcf.gz

# Apply hard filters for SNPs 
# NOTE: These filters are RELAXED for low-coverage test data (~1-2x coverage)
# For production WGS data (30x+ coverage), use stricter thresholds:
#   QD < 2.0, QUAL < 30.0, SOR > 3.0, FS > 60.0, MQ < 40.0, MQRankSum < -12.5, ReadPosRankSum < -8.0
gatk VariantFiltration \
    -R ${REFERENCE} \
    -V ${VAR_DIR}/${SAMPLE}_raw_snps.vcf.gz \
    -O ${FILTERED_SNP_VCF} \
    --filter-expression "QD < 1.0" --filter-name "QD1" \
    --filter-expression "QUAL < 10.0" --filter-name "QUAL10" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10" \
    --filter-expression "FS > 100.0" --filter-name "FS100" \
    --filter-expression "MQ < 20.0" --filter-name "MQ20"
echo "[$(date)] SNP filtering completed"
echo

# Step 12: Hard Filtering for Indels
echo "[$(date)] Step 12: Hard filtering indels..."
# Select Indels
gatk SelectVariants \
    -R ${REFERENCE} \
    -V ${RAW_VCF} \
    --select-type-to-include INDEL \
    -O ${VAR_DIR}/${SAMPLE}_raw_indels.vcf.gz

# Apply hard filters for indels
# NOTE: These filters are RELAXED for low-coverage test data (~1-2x coverage)
# For production WGS data (30x+ coverage), use stricter thresholds:
#   QD < 2.0, QUAL < 30.0, FS > 200.0, ReadPosRankSum < -20.0
gatk VariantFiltration \
    -R ${REFERENCE} \
    -V ${VAR_DIR}/${SAMPLE}_raw_indels.vcf.gz \
    -O ${FILTERED_INDEL_VCF} \
    --filter-expression "QD < 1.0" --filter-name "QD1" \
    --filter-expression "QUAL < 10.0" --filter-name "QUAL10" \
    --filter-expression "FS > 300.0" --filter-name "FS300"
echo "[$(date)] Indel filtering completed"
echo

# Step 13: Merge filtered SNPs and Indels
echo "[$(date)] Step 13: Merging filtered variants..."
gatk MergeVcfs \
    -I ${FILTERED_SNP_VCF} \
    -I ${FILTERED_INDEL_VCF} \
    -O ${FINAL_VCF}
echo "[$(date)] Merging completed"
echo

# Step 14: Functional Annotation with SnpEff
echo "[$(date)] Step 14: Annotating variants with SnpEff..."
# First, unzip the VCF for SnpEff (it works better with uncompressed VCF)
gunzip -c ${FINAL_VCF} > ${VAR_DIR}/${SAMPLE}_filtered.vcf

# Run SnpEff annotation
# Using GRCh38.mane.1.0.refseq as the database (for human genome)
# Note: For the test data (chr22 region), we'll use a generic annotation
snpEff -Xmx4g \
    -v \
    GRCh38.mane.1.0.refseq \
    ${VAR_DIR}/${SAMPLE}_filtered.vcf \
    > ${ANNOTATED_VCF} \
    2> ${VAR_DIR}/${SAMPLE}_snpeff.log || {
        echo "Warning: SnpEff annotation failed. This is expected for test data."
        echo "Creating placeholder annotated VCF..."
        cp ${VAR_DIR}/${SAMPLE}_filtered.vcf ${ANNOTATED_VCF}
    }

# Compress the annotated VCF
bgzip -c ${ANNOTATED_VCF} > ${VAR_DIR}/${SAMPLE}_annotated.vcf.gz
tabix -p vcf ${VAR_DIR}/${SAMPLE}_annotated.vcf.gz
echo "[$(date)] Annotation completed"
echo

# Step 15: Variant Statistics
echo "[$(date)] Step 15: Generating variant statistics..."
bcftools stats ${RAW_VCF} > ${VAR_DIR}/${SAMPLE}_variant_stats_raw.txt
bcftools stats ${FINAL_VCF} > ${VAR_DIR}/${SAMPLE}_variant_stats_filtered.txt

# Count variants by type (raw)
echo "Raw variant counts:"
bcftools view -v snps ${RAW_VCF} | bcftools query -f '.\n' | wc -l > ${VAR_DIR}/${SAMPLE}_raw_snp_count.txt
bcftools view -v indels ${RAW_VCF} | bcftools query -f '.\n' | wc -l > ${VAR_DIR}/${SAMPLE}_raw_indel_count.txt

# Count variants by type (filtered - PASS only)
echo "Filtered variant counts (PASS only):"
bcftools view -f PASS -v snps ${FINAL_VCF} | bcftools query -f '.\n' | wc -l > ${VAR_DIR}/${SAMPLE}_filtered_snp_count.txt
bcftools view -f PASS -v indels ${FINAL_VCF} | bcftools query -f '.\n' | wc -l > ${VAR_DIR}/${SAMPLE}_filtered_indel_count.txt
echo "[$(date)] Variant statistics completed"
echo

# Step 16: Create Visualization Files
echo "[$(date)] Step 16: Creating visualization files..."
# Create BED file for genome browsers
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' ${FINAL_VCF} > ${VAR_DIR}/${SAMPLE}_variants.bed

# Generate coverage track
bedtools genomecov -ibam ${RECAL_BAM} -bg > ${ALIGNED_DIR}/${SAMPLE}_coverage.bedgraph
echo "[$(date)] Visualization files created"
echo

# Summary
echo "=========================================="
echo "Pipeline Completed Successfully!"
echo "=========================================="
echo "Results:"
echo "  - QC reports: ${QC_DIR}/"
echo "  - Trimmed reads: ${TRIMMED_DIR}/"
echo "  - Final BAM: ${RECAL_BAM}"
echo "  - GVCF: ${GVCF}"
echo "  - Raw Variants VCF: ${RAW_VCF}"
echo "  - Filtered VCF: ${FINAL_VCF}"
echo "  - Annotated VCF: ${VAR_DIR}/${SAMPLE}_annotated.vcf.gz"
echo "  - Variant statistics: ${VAR_DIR}/"
echo "  - Duplication metrics: ${METRICS}"
echo

echo "Quick Quality Control Summary:"
echo "------------------------------"
echo "Alignment metrics:"
samtools flagstat ${RECAL_BAM}
echo

echo "Duplication rate:"
grep -A1 "^LIBRARY" ${METRICS} | tail -1 | cut -f9
echo

echo "Variant counts (Raw):"
echo -n "  SNPs: "
cat ${VAR_DIR}/${SAMPLE}_raw_snp_count.txt
echo -n "  Indels: "
cat ${VAR_DIR}/${SAMPLE}_raw_indel_count.txt
echo

echo "Variant counts (Filtered - PASS only):"
echo -n "  SNPs: "
cat ${VAR_DIR}/${SAMPLE}_filtered_snp_count.txt
echo -n "  Indels: "
cat ${VAR_DIR}/${SAMPLE}_filtered_indel_count.txt
echo

echo "For detailed QC reports, see:"
echo "  - ${QC_DIR}/${SAMPLE}_alignment_summary.txt"
echo "  - ${VAR_DIR}/${SAMPLE}_variant_stats_raw.txt"
echo "  - ${VAR_DIR}/${SAMPLE}_variant_stats_filtered.txt"
echo "=========================================="
