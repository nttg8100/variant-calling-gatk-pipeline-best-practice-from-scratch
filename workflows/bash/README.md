# Bash GATK Variant Calling Pipeline

This directory contains the **bash implementation** of the 16-step GATK germline variant calling workflow. This serves as the educational baseline and validation reference for the Nextflow implementation.

## Purpose

✅ **Learning tool** - Understand GATK workflow step-by-step  
✅ **Validation baseline** - Establish ground truth for Nextflow migration  
✅ **Quick prototyping** - Test parameters on 1-3 samples before scaling  
✅ **Debugging reference** - Troubleshoot issues by comparing with Nextflow

## Files

- `gatk_pipeline.sh` - Complete 16-step GATK workflow
- `download_data.sh` - Download test data and reference genome (deprecated, use `../../scripts/download_data.sh`)

## Quick Start

### 1. Install Dependencies

From the **repository root** directory:

```bash
cd /path/to/variant-calling-gatk-pipeline-best-practice-from-scratch

# Install tools with Pixi
pixi install

# Verify installation
pixi run bwa 2>&1 | head -3
pixi run gatk --version
```

### 2. Download Test Data

```bash
# From repository root
pixi run bash scripts/download_data.sh
```

This creates:
- `data/` - Sample FASTQ files (sample1/2/3)
- `reference/` - hg38 chr22 subset + BWA indices + known sites

### 3. Run the Pipeline

```bash
cd workflows/bash

# Run single sample (sample1 by default)
pixi run bash gatk_pipeline.sh 2>&1 | tee pipeline.log

# Run with custom sample
SAMPLE=sample2 pixi run bash gatk_pipeline.sh 2>&1 | tee pipeline_sample2.log
```

**Expected runtime**: ~3-4 minutes per sample on modern hardware

## Pipeline Steps

The script implements all 16 GATK best practice steps:

### Pre-processing (Steps 1-8)
1. **FastQC** - Quality control metrics
2. **Trim Galore** - Adapter trimming + quality filtering
3. **BWA-MEM** - Read alignment to reference genome
4. **SAMtools sort** - Sort BAM by genomic coordinate
5. **GATK MarkDuplicates** - Mark PCR/optical duplicates
6. **GATK BaseRecalibrator** - Model systematic errors (BQSR)
7. **GATK ApplyBQSR** - Apply recalibration to base qualities
8. **GATK CollectMetrics** - Alignment QC + insert size histogram (requires R)

### Variant Calling (Steps 9-10)
9. **GATK HaplotypeCaller** - Call variants in GVCF mode
10. **GATK GenotypeGVCFs** - Joint genotyping across samples

### Variant Filtering (Steps 11-13)
11. **SelectVariants + VariantFiltration (SNPs)** - Extract and filter SNPs
12. **SelectVariants + VariantFiltration (Indels)** - Extract and filter indels
13. **GATK MergeVcfs** - Merge filtered SNPs + indels

### Annotation & Statistics (Steps 14-16)
14. **SnpEff** - Functional annotation (gene, transcript, effect)
15. **bcftools stats** - Variant statistics (raw/filtered counts)
16. **Visualization** - BED file + bedGraph coverage track

## Output Structure

After running the pipeline, results are saved in `results/`:

```
results/
├── aligned/
│   ├── sample1_aligned.bam        # Raw aligned BAM
│   ├── sample1_sorted.bam         # Coordinate-sorted BAM
│   ├── sample1_dedup.bam + .bai   # Deduplicated BAM
│   ├── sample1_recal.bam + .bai   # Final analysis-ready BAM
│   ├── sample1_metrics.txt        # Duplication metrics
│   └── sample1_coverage.bedgraph  # Coverage track
├── bqsr/
│   └── sample1_recal.table        # BQSR recalibration table
├── qc/
│   ├── sample1_R1_fastqc.html/zip # FastQC reports
│   ├── sample1_R2_fastqc.html/zip
│   ├── sample1_alignment_summary.txt    # Alignment metrics
│   ├── sample1_insert_size_metrics.txt  # Insert size stats
│   └── sample1_insert_size_histogram.pdf # Insert size plot
├── trimmed/
│   ├── sample1_R1_val_1.fq.gz     # Trimmed forward reads
│   └── sample1_R2_val_2.fq.gz     # Trimmed reverse reads
└── variants/
    ├── sample1.g.vcf.gz + .tbi    # GVCF (HaplotypeCaller output)
    ├── sample1_raw.vcf.gz + .tbi  # Raw variants (GenotypeGVCFs)
    ├── sample1_raw_snps.vcf.gz    # Raw SNPs only
    ├── sample1_filtered_snps.vcf.gz # Filtered SNPs (PASS variants)
    ├── sample1_raw_indels.vcf.gz  # Raw indels only
    ├── sample1_filtered_indels.vcf.gz # Filtered indels
    ├── sample1_filtered.vcf.gz    # Final merged VCF (SNPs + indels)
    ├── sample1_annotated.vcf      # SnpEff annotated VCF
    ├── sample1_variant_stats_raw.txt      # bcftools stats (raw)
    ├── sample1_variant_stats_filtered.txt # bcftools stats (filtered)
    ├── sample1_raw_snp_count.txt          # SNP count
    ├── sample1_raw_indel_count.txt        # Indel count
    ├── sample1_filtered_snp_count.txt     # Filtered SNP count
    └── sample1_filtered_indel_count.txt   # Filtered indel count
```

## Configuration

Edit `gatk_pipeline.sh` to customize:

```bash
# Sample name (or set via environment variable)
SAMPLE="${SAMPLE:-sample1}"

# Thread count
THREADS=2

# Paths (relative to repository root)
DATA_DIR="../../data"
REF_DIR="../../reference"
RESULTS_DIR="./results"

# Reference files
REFERENCE="${REF_DIR}/genome.fasta"
DBSNP="${REF_DIR}/dbsnp_146.hg38.vcf.gz"
MILLS="${REF_DIR}/mills_and_1000G.indels.vcf.gz"
```

## Validation

Check outputs after pipeline completion:

```bash
# Check BAM statistics
pixi run samtools flagstat results/aligned/sample1_recal.bam

# View duplication metrics
cat results/aligned/sample1_metrics.txt | grep -A 2 "LIBRARY"

# Count variants
pixi run bash -c "bcftools view -H results/variants/sample1_filtered.vcf.gz | wc -l"

# View variant statistics
cat results/variants/sample1_variant_stats_filtered.txt
```

**Expected for test data (chr22 subset):**
- Total reads: ~533k
- Duplicates: ~0.7%
- Mapping rate: ~1.9% (low because reference is only 40kb)
- Variants: 0 (expected for small test region with low coverage)

## Performance Benchmarking

| Step | Duration | % of Total |
|------|----------|------------|
| FastQC | 8s | 6.6% |
| Trim Galore | 12s | 9.9% |
| BWA-MEM | 23s | 18.9% |
| SAMtools sort | 4s | 3.3% |
| MarkDuplicates | 15s | 12.3% |
| BaseRecalibrator | 21s | 17.2% |
| ApplyBQSR | 24s | 19.7% |
| CollectMetrics | 8s | 6.6% |
| HaplotypeCaller | 13s | 10.7% |
| GenotypeGVCFs | 3s | 2.5% |
| Filtering | 5s | 4.1% |
| Annotation | 4s | 3.3% |
| Statistics | 2s | 1.6% |
| Visualization | 1s | 0.8% |
| **Total** | **~3m 38s** | **100%** |

## Limitations

### ⚠️ Known Limitations

1. **No parallelization** - Processes one sample at a time
   - 3 samples = ~11 minutes (serial)
   - Nextflow = 3m 26s (3.2x faster)

2. **No resume capability** - Must restart entire pipeline on failure
   - Nextflow supports `-resume` to restart from last successful step

3. **Hard-coded resources** - Fixed thread count (THREADS=2)
   - Nextflow dynamically allocates resources based on `process.cpus`

4. **Manual scaling** - Requires custom loops for multiple samples
   - Nextflow automatically parallelizes via samplesheet CSV

5. **Limited error handling** - Basic `set -euo pipefail`
   - Nextflow provides detailed error messages and work directory debugging

6. **No container support** - Relies on Pixi environment
   - Nextflow uses Singularity/Docker for full reproducibility

### When to Use Bash vs. Nextflow

**Use Bash when:**
- Learning GATK workflow step-by-step
- Prototyping/testing on 1-3 samples
- Establishing validation baseline
- Running on local machine without containers

**Use Nextflow when:**
- Processing 10+ samples in production
- Deploying to HPC clusters (SLURM) or cloud (AWS/GCP)
- Need resume capability after failures
- Require full reproducibility with containers

## Troubleshooting

### Issue: "Command not found" errors

**Cause**: Tools not installed or Pixi environment not activated

**Solution**:
```bash
# From repository root
pixi install
pixi run bash workflows/bash/gatk_pipeline.sh
```

### Issue: "File not found" errors for reference/data files

**Cause**: Test data not downloaded

**Solution**:
```bash
# From repository root
pixi run bash scripts/download_data.sh
```

### Issue: GATK CollectMetrics fails with "RScript not found"

**Cause**: R runtime not available in environment

**Solution**:
```bash
# Ensure R is installed via Pixi (should be in pixi.toml)
pixi run R --version

# If missing, add to pixi.toml:
# dependencies = { r-base = ">=4.0" }
```

### Issue: Zero variants in output VCF

**Cause**: Test data is small chr22 subset with low coverage

**Expected**: This is normal for test data. Real data (30x coverage, full genome) will produce thousands of variants.

### Issue: Low mapping rate (~1.9%)

**Cause**: Reference genome is only 40kb region of chr22

**Expected**: This is normal for test data. Most reads map outside this small region.

## Comparing with Nextflow

To validate that Bash and Nextflow produce equivalent results:

```bash
# Run bash pipeline
cd workflows/bash
pixi run bash gatk_pipeline.sh 2>&1 | tee bash_pipeline.log

# Run nextflow pipeline
cd ../nextflow
nextflow run main.nf -profile singularity,test -resume

# Compare variant statistics
diff bash/results/variants/sample1_variant_stats_filtered.txt \
     nextflow/results_nextflow/variants/sample1/sample1_variant_stats_filtered.txt

# Expected: Variant counts should be identical
# MD5 checksums will differ due to timestamps in BAM/VCF headers (this is normal)
```

## References

- **Blog post**: [Part 1 - Building a Production-Ready GATK Bash Workflow with Pixi](https://riverxdata.github.io/river-docs/blog/gatk-variant-calling-bash-workflow-pixi-part1)
- **GATK documentation**: https://gatk.broadinstitute.org/
- **Repository root README**: `../../README.md`

## Support

For issues or questions:
- Open an issue on [GitHub](https://github.com/nttg8100/variant-calling-gatk-pipeline-best-practice-from-scratch/issues)
- Check the troubleshooting guide: `../../docs/TROUBLESHOOTING.md`
- Read the blog post: https://riverxdata.github.io/river-docs/

---

**Status**: ✅ Complete (16/16 steps implemented)  
**Last Updated**: February 2026
