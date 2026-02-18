#!/bin/bash
#
# Download test data for GATK variant calling pipeline
# This script downloads reference genome, known sites, and test FASTQ files from nf-core
#

set -euo pipefail

echo "=========================================="
echo "Downloading GATK Pipeline Test Data"
echo "=========================================="
echo

# Create directories
echo "Creating directories..."
mkdir -p data reference
echo "✓ Directories created"
echo

# Download reference genome
echo "Downloading reference genome..."
cd reference

curl -L -o genome.fasta \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta
echo "✓ genome.fasta downloaded"

curl -L -o genome.fasta.fai \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta.fai
echo "✓ genome.fasta.fai downloaded"

curl -L -o genome.dict \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.dict
echo "✓ genome.dict downloaded"

# Download known sites for BQSR
echo ""
echo "Downloading known sites VCF files..."

curl -L -o dbsnp_146.hg38.vcf.gz \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz
echo "✓ dbsnp_146.hg38.vcf.gz downloaded"

curl -L -o dbsnp_146.hg38.vcf.gz.tbi \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz.tbi
echo "✓ dbsnp_146.hg38.vcf.gz.tbi downloaded"

curl -L -o mills_and_1000G.indels.vcf.gz \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/vcf/mills_and_1000G.indels.vcf.gz
echo "✓ mills_and_1000G.indels.vcf.gz downloaded"

curl -L -o mills_and_1000G.indels.vcf.gz.tbi \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/vcf/mills_and_1000G.indels.vcf.gz.tbi
echo "✓ mills_and_1000G.indels.vcf.gz.tbi downloaded"

cd ..

# Index the reference genome
echo ""
echo "Creating BWA index for reference genome..."
echo "This may take a few moments..."
bwa index reference/genome.fasta
echo "✓ BWA index created"

# Download test FASTQ files
echo ""
echo "Downloading test FASTQ files (3 samples)..."
cd data

# Sample 1 (test data from nf-core)
echo "Downloading sample1..."
curl -L -o sample1_R1.fastq.gz \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz
echo "✓ sample1_R1.fastq.gz downloaded"

curl -L -o sample1_R2.fastq.gz \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_2.fastq.gz
echo "✓ sample1_R2.fastq.gz downloaded"

# Sample 2 (copy of sample1 for multi-sample testing)
echo "Creating sample2 (copy of sample1 for testing)..."
cp sample1_R1.fastq.gz sample2_R1.fastq.gz
cp sample1_R2.fastq.gz sample2_R2.fastq.gz
echo "✓ sample2 created"

# Sample 3 (copy of sample1 for multi-sample testing)
echo "Creating sample3 (copy of sample1 for testing)..."
cp sample1_R1.fastq.gz sample3_R1.fastq.gz
cp sample1_R2.fastq.gz sample3_R2.fastq.gz
echo "✓ sample3 created"

cd ..

# Summary
echo ""
echo "=========================================="
echo "Data Download Completed Successfully!"
echo "=========================================="
echo ""
echo "Downloaded files:"
echo "Reference genome:"
ls -lh reference/genome.fasta* reference/genome.dict 2>/dev/null || true
echo ""
echo "BWA indices:"
ls -lh reference/genome.fasta.{amb,ann,bwt,pac,sa} 2>/dev/null || true
echo ""
echo "Known sites VCFs:"
ls -lh reference/*.vcf.gz* 2>/dev/null || true
echo ""
echo "Test FASTQ files (3 samples):"
ls -lh data/*.fastq.gz 2>/dev/null || true
echo ""
echo "=========================================="
echo "You can now run the pipelines:"
echo ""
echo "Bash pipeline (single sample):"
echo "  cd workflows/bash"
echo "  pixi run bash gatk_pipeline.sh"
echo ""
echo "Nextflow pipeline (single sample):"
echo "  cd workflows/nextflow"
echo "  nextflow run main.nf -profile singularity,test -resume"
echo ""
echo "Nextflow pipeline (multi-sample):"
echo "  cd workflows/nextflow"
echo "  nextflow run main.nf -profile singularity --input samplesheet.csv -resume"
echo "=========================================="
