# GATK Variant Calling Pipeline - Best Practice from Scratch

A complete GATK variant calling workflow demonstrating the transformation from bash scripts to enterprise-grade Nextflow pipelines. This repository serves as both a learning resource and production baseline for genomic variant calling.

## Project Goals

1. **Demonstrate GATK Best Practices**: Implement the standard germline variant calling workflow following GATK guidelines
2. **Compare Approaches**: Show the evolution from bash scripts to Nextflow workflows
3. **Enable Reproducibility**: Use modern tools (Pixi, containers) for environment management
4. **Production-Ready**: Build towards scalable, testable, and maintainable pipelines

## Repository Structure

```
variant-calling-gatk-pipeline-best-practice-from-scratch/
├── pixi.toml                    # Environment specification
├── download_data.sh             # Script to download test data
├── gatk_pipeline.sh             # Bash workflow implementation
├── README.md                    # This file
├── .gitignore                   # Git ignore patterns
├── data/                        # Input FASTQ files (not tracked)
├── reference/                   # Reference genome and indices (not tracked)
└── results/                     # Pipeline outputs (not tracked)
```

## Quick Start

### Prerequisites

- **Linux** system (tested on Ubuntu/Debian)
- **Pixi** package manager ([installation guide](https://pixi.sh))
- **8GB RAM** minimum (16GB recommended)
- **10GB disk space** for test data and results

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/variant-calling-gatk-pipeline-best-practice-from-scratch.git
cd variant-calling-gatk-pipeline-best-practice-from-scratch

# Install dependencies with Pixi
pixi install

# Verify installation
pixi run bwa 2>&1 | head -3
pixi run samtools --version
pixi run gatk --version
pixi run fastqc --version
```

### Download Test Data

We provide a convenient script to download all required test data:

```bash
# Run the download script (uses Pixi environment for BWA indexing)
pixi run bash download_data.sh
```

This script will:
- Create `data/` and `reference/` directories
- Download reference genome (genome.fasta) and indices
- Download known sites VCF files (dbSNP and Mills indels)
- Create BWA index for the reference genome
- Download test FASTQ files (test_1.fastq.gz, test_2.fastq.gz)

**Estimated download time**: 1-2 minutes  
**Total download size**: ~35 MB

### Run the Pipeline

```bash
# Make the script executable
chmod +x gatk_pipeline.sh

# Execute the workflow
pixi run bash gatk_pipeline.sh 2>&1 | tee pipeline.log
```

**Expected runtime**: ~2 minutes on modern hardware

## Workflow Overview

The pipeline implements the GATK germline variant calling best practices:

### Pipeline Steps

1. **Quality Control (FastQC)**
   - Assess read quality metrics
   - Identify potential sequencing issues
   - Generate HTML reports

2. **Read Alignment (BWA-MEM)**
   - Align paired-end reads to reference genome
   - Add read group information
   - Output in BAM format

3. **BAM Sorting (SAMtools)**
   - Sort alignments by genomic coordinate
   - Required for downstream GATK tools

4. **BAM Indexing (SAMtools)**
   - Create index for efficient BAM access

5. **Mark Duplicates (GATK MarkDuplicates)**
   - Identify PCR and optical duplicates
   - Generate duplication metrics

6. **Base Quality Score Recalibration - Table (GATK BaseRecalibrator)**
   - Model systematic errors in base quality scores
   - Use known variants to train model

7. **Apply BQSR (GATK ApplyBQSR)**
   - Apply recalibration to base qualities
   - Produce analysis-ready BAM

8. **Variant Calling (GATK HaplotypeCaller)**
   - Call SNPs and indels via local assembly
   - Generate VCF output

9. **VCF Indexing (GATK IndexFeatureFile)**
   - Create index for efficient VCF access

### Performance Metrics (Test Data)

| Step             | Duration  | Percentage |
| ---------------- | --------- | ---------- |
| FastQC           | 8s        | 6.6%       |
| Alignment (BWA)  | 23s       | 18.9%      |
| Sorting          | 4s        | 3.3%       |
| MarkDuplicates   | 15s       | 12.3%      |
| BaseRecalibrator | 21s       | 17.2%      |
| ApplyBQSR        | 24s       | 19.7%      |
| HaplotypeCaller  | 13s       | 10.7%      |
| VCF Indexing     | 13s       | 10.7%      |
| **Total**        | **~122s** | **100%**   |

## Tool Versions

Installed via Pixi (locked in `pixi.lock`):

- **BWA**: 0.7.19-r1273
- **SAMtools**: 1.23
- **GATK**: 4.6.2.0
- **FastQC**: v0.12.1

## Results and Validation

### Expected Outputs

After successful execution, you'll find:

```
results/
├── fastqc/
│   ├── test_1_fastqc.html       # Read 1 QC report
│   ├── test_1_fastqc.zip
│   ├── test_2_fastqc.html       # Read 2 QC report
│   └── test_2_fastqc.zip
├── test.bam                     # Raw aligned BAM
├── test.sorted.bam              # Sorted BAM
├── test.sorted.bam.bai          # BAM index
├── test.dedup.bam               # Deduplicated BAM
├── test.dedup.bai               # Dedup BAM index
├── test.metrics.txt             # Duplication metrics
├── test.recal.table             # BQSR calibration table
├── test.recal.bam               # Final analysis-ready BAM
├── test.recal.bai               # Final BAM index
├── test.vcf                     # Variant calls
└── test.vcf.idx                 # VCF index
```

### Validation

```bash
# Check BAM statistics
pixi run samtools flagstat results/test.recal.bam

# Expected output:
# 533,554 total reads
# 3,876 duplicates (0.73%)
# 10,125 mapped (1.90%)

# View duplication metrics
cat results/test.metrics.txt | grep -A 2 "LIBRARY"

# Check VCF structure
pixi run bash -c "grep -v '^##' results/test.vcf | head"
```

**Note**: Low mapping rate (1.90%) is expected—the test reference is only a 40kb region of chr22.

## Known Limitations

### Current Bash Implementation

1. **No Parallelization**: Processes one sample at a time
   - 1 sample: 2 minutes
   - 10 samples: 20 minutes (serial)
   - 100 samples: 3.3 hours (serial)

2. **Hard-coded Resources**: Fixed thread count (THREADS=2)
   - No dynamic allocation based on system resources
   - No memory limits or monitoring

3. **Limited Error Handling**: Basic `set -euo pipefail`
   - No resume capability
   - Must restart entire pipeline on failure
   - Limited error context

4. **Manual Scaling**: Requires custom loops for multiple samples
   - No batch processing support
   - Difficult to track multiple samples
   - No progress reporting

5. **No Validation**: Manual output inspection required
   - No automated QC checks
   - No regression testing
   - No comparison to baseline results

6. **Portability Issues**: Hard-coded paths and assumptions
   - Environment-dependent execution
   - No container support
   - "Works on my machine" problems

## Roadmap: Nextflow Transformation

This bash implementation serves as a baseline for transformation to enterprise-grade Nextflow:

### Phase 1: Core Nextflow Migration (Upcoming)
- [ ] Convert bash script to Nextflow DSL2
- [ ] Create modular processes for each step
- [ ] Implement channel-based data flow
- [ ] Add configuration management
- [ ] Enable parallel sample processing

### Phase 2: Testing & Validation
- [ ] Implement nf-test framework
- [ ] Create unit tests for each process
- [ ] Add integration tests
- [ ] Validate against bash baseline (byte-for-byte)
- [ ] Establish regression testing

### Phase 3: Quality & Standards
- [ ] Apply nf-lint for code quality
- [ ] Follow nf-core best practices
- [ ] Add comprehensive documentation
- [ ] Create reusable modules
- [ ] Implement nf-core schema

### Phase 4: CI/CD Integration
- [ ] Set up GitHub Actions workflows
- [ ] Automate testing on multiple platforms
- [ ] Add Docker/Singularity containers
- [ ] Implement release management
- [ ] Add performance benchmarking

### Phase 5: Production Deployment
- [ ] Configure for HPC (SLURM/PBS/SGE)
- [ ] Add cloud execution profiles (AWS, GCP, Azure)
- [ ] Optimize resource allocation
- [ ] Implement monitoring and logging
- [ ] Add multi-sample processing

### Success Criteria

The Nextflow version must:

✅ Produce **identical results** to bash version (validated via checksums)  
✅ Support **multiple samples** in parallel  
✅ **Resume** from any failure point  
✅ Run on **local, HPC, and cloud** environments  
✅ Process **10 samples in <5 minutes** (vs. 20 min with bash)  
✅ Pass **nf-lint** standards  
✅ Achieve **100% test coverage** with nf-test  
✅ Include **automated CI/CD** validation  

## Related Resources

### Blog Posts
- [Building a Production-Ready GATK Variant Calling Bash Workflow with Pixi (Part 1)](https://riverxdata.com/blog/gatk-variant-calling-bash-workflow-pixi-part1)
- [Pixi - New conda era](https://riverxdata.com/blog/pixi-is-new-conda-based-era)
- [Containers on HPC: From Docker to Singularity and Apptainer](https://riverxdata.com/blog/containers-hpc-docker-singularity-apptainer)

### External References
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- [NGS101 GATK Tutorial](https://ngs101.com/how-to-analyze-whole-genome-sequencing-data-for-absolute-beginners-part-1-from-raw-reads-to-high-quality-variants-using-gatk/)
- [nf-core Test Datasets](https://github.com/nf-core/test-datasets)
- [Pixi Documentation](https://pixi.sh)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

### Development Workflow

```bash
# Create a feature branch
git checkout -b feat/my-feature

# Make changes and test
pixi run bash gatk_pipeline.sh 2>&1 | tee test.log

# Commit with conventional commits
git commit -m "feat: add new feature"

# Push and create PR
git push origin feat/my-feature
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- **GATK Team** at Broad Institute for best practices guidelines
- **nf-core Community** for test datasets and standards
- **NGS101** for excellent beginner tutorials
- **Pixi Team** for modern package management

## Contact

For questions or feedback:
- Open an issue on GitHub
- Visit [RiverXData Blog](https://riverxdata.com/blog)

---

**Status**: Active Development | Bash Implementation Complete | Nextflow Migration Coming Soon
