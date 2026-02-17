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
echo "Downloading test FASTQ files..."
cd data

curl -L -o test_1.fastq.gz \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz
echo "✓ test_1.fastq.gz downloaded"

curl -L -o test_2.fastq.gz \
  https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_2.fastq.gz
echo "✓ test_2.fastq.gz downloaded"

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
echo "Test FASTQ files:"
ls -lh data/*.fastq.gz 2>/dev/null || true
echo ""
echo "=========================================="
echo "You can now run the pipeline with:"
echo "  pixi run bash gatk_pipeline.sh"
echo "=========================================="
