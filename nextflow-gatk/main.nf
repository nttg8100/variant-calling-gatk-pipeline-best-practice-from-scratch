#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    GATK Variant Calling Pipeline - Complete 16-Step Nextflow Version
========================================================================================
    Based on the bash workflow - Full migration from Part 1
========================================================================================
*/

// Include modules
include { FASTQC } from './modules/fastqc'
include { TRIM_GALORE } from './modules/trim_galore'
include { BWA_MEM } from './modules/bwa_mem'
include { SAMTOOLS_SORT } from './modules/samtools_sort'
include { GATK_MARKDUPLICATES } from './modules/gatk_markduplicates'
include { GATK_BASERECALIBRATOR } from './modules/gatk_baserecalibrator'
include { GATK_APPLYBQSR } from './modules/gatk_applybqsr'
include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller'
include { GATK_GENOTYPEGVCFS } from './modules/gatk_genotypegvcfs'
include { GATK_SELECTVARIANTS_SNP } from './modules/gatk_selectvariants_snp'
include { GATK_VARIANTFILTRATION_SNP } from './modules/gatk_variantfiltration_snp'
include { GATK_SELECTVARIANTS_INDEL } from './modules/gatk_selectvariants_indel'
include { GATK_VARIANTFILTRATION_INDEL } from './modules/gatk_variantfiltration_indel'
include { GATK_MERGEVCFS } from './modules/gatk_mergevcfs'

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow GATK_VARIANT_CALLING {
    
    take:
    reads_ch        // channel: [ val(meta), [ path(read1), path(read2) ] ]
    reference_ch    // channel: [ path(fasta), path(fai), path(dict) ]
    bwa_index_ch    // channel: [ path(amb), path(ann), path(bwt), path(pac), path(sa) ]
    dbsnp_ch        // channel: [ path(vcf), path(tbi) ]
    known_indels_ch // channel: [ path(vcf), path(tbi) ]
    
    main:
    ch_versions = Channel.empty()
    
    //
    // STEP 1: Quality Control with FastQC
    //
    FASTQC (
        reads_ch
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    
    //
    // STEP 2: Adapter Trimming and Quality Filtering
    //
    TRIM_GALORE (
        reads_ch
    )
    ch_versions = ch_versions.mix(TRIM_GALORE.out.versions)
    
    //
    // STEP 3: Read Alignment with BWA-MEM
    //
    BWA_MEM (
        TRIM_GALORE.out.reads,
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] },
        bwa_index_ch
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    
    //
    // STEP 4: Sort BAM file
    //
    SAMTOOLS_SORT (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    
    //
    // STEP 5: Mark Duplicates
    //
    GATK_MARKDUPLICATES (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(GATK_MARKDUPLICATES.out.versions)
    
    //
    // STEP 6: Base Quality Score Recalibration - Generate table
    //
    GATK_BASERECALIBRATOR (
        GATK_MARKDUPLICATES.out.bam.join(GATK_MARKDUPLICATES.out.bai),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] },
        dbsnp_ch.map { it[0] },
        dbsnp_ch.map { it[1] },
        known_indels_ch.map { it[0] },
        known_indels_ch.map { it[1] }
    )
    ch_versions = ch_versions.mix(GATK_BASERECALIBRATOR.out.versions)
    
    //
    // STEP 7: Apply BQSR
    //
    GATK_APPLYBQSR (
        GATK_MARKDUPLICATES.out.bam
            .join(GATK_MARKDUPLICATES.out.bai)
            .join(GATK_BASERECALIBRATOR.out.table),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] }
    )
    ch_versions = ch_versions.mix(GATK_APPLYBQSR.out.versions)
    
    //
    // STEP 9: Variant Calling with HaplotypeCaller (GVCF mode)
    //
    GATK_HAPLOTYPECALLER (
        GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] },
        dbsnp_ch.map { it[0] },
        dbsnp_ch.map { it[1] }
    )
    ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER.out.versions)
    
    //
    // STEP 10: Genotype GVCFs
    //
    GATK_GENOTYPEGVCFS (
        GATK_HAPLOTYPECALLER.out.gvcf.join(GATK_HAPLOTYPECALLER.out.tbi),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] }
    )
    ch_versions = ch_versions.mix(GATK_GENOTYPEGVCFS.out.versions)
    
    //
    // STEP 11: Select and Filter SNPs
    //
    GATK_SELECTVARIANTS_SNP (
        GATK_GENOTYPEGVCFS.out.vcf.join(GATK_GENOTYPEGVCFS.out.tbi),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] }
    )
    ch_versions = ch_versions.mix(GATK_SELECTVARIANTS_SNP.out.versions)
    
    GATK_VARIANTFILTRATION_SNP (
        GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.tbi),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] }
    )
    ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION_SNP.out.versions)
    
    //
    // STEP 12: Select and Filter Indels
    //
    GATK_SELECTVARIANTS_INDEL (
        GATK_GENOTYPEGVCFS.out.vcf.join(GATK_GENOTYPEGVCFS.out.tbi),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] }
    )
    ch_versions = ch_versions.mix(GATK_SELECTVARIANTS_INDEL.out.versions)
    
    GATK_VARIANTFILTRATION_INDEL (
        GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.tbi),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] }
    )
    ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION_INDEL.out.versions)
    
    //
    // STEP 13: Merge filtered SNPs and Indels
    //
    GATK_MERGEVCFS (
        GATK_VARIANTFILTRATION_SNP.out.vcf
            .join(GATK_VARIANTFILTRATION_SNP.out.tbi)
            .join(GATK_VARIANTFILTRATION_INDEL.out.vcf)
            .join(GATK_VARIANTFILTRATION_INDEL.out.tbi),
        reference_ch.map { it[0] },
        reference_ch.map { it[1] },
        reference_ch.map { it[2] }
    )
    ch_versions = ch_versions.mix(GATK_MERGEVCFS.out.versions)
    
    emit:
    fastqc_html = FASTQC.out.html
    fastqc_zip  = FASTQC.out.zip
    trimmed_reads = TRIM_GALORE.out.reads
    final_bam = GATK_APPLYBQSR.out.bam
    gvcf = GATK_HAPLOTYPECALLER.out.gvcf
    raw_vcf = GATK_GENOTYPEGVCFS.out.vcf
    final_vcf = GATK_MERGEVCFS.out.vcf
    versions = ch_versions
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    //
    // Create input channel from samplesheet or input parameters
    //
    if (params.input) {
        // Read from samplesheet CSV
        Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                def meta = [:]
                meta.id = row.sample
                def reads = []
                reads.add(file(row.fastq_1))
                if (row.fastq_2) {
                    reads.add(file(row.fastq_2))
                }
                return [ meta, reads ]
            }
            .set { ch_input }
    } else {
        // Direct parameters
        def meta = [:]
        meta.id = params.sample ?: 'sample1'
        
        def reads = []
        reads.add(file(params.fastq_r1))
        if (params.fastq_r2) {
            reads.add(file(params.fastq_r2))
        }
        
        ch_input = Channel.of([ meta, reads ])
    }
    
    //
    // Prepare reference genome channel
    //
    reference_ch = Channel.of([
        file(params.reference),
        file("${params.reference}.fai"),
        file(params.reference.toString().replace('.fasta', '.dict').replace('.fa', '.dict'))
    ])
    
    //
    // Prepare BWA index files channel
    //
    bwa_index_ch = Channel.of([
        file("${params.reference}.amb"),
        file("${params.reference}.ann"),
        file("${params.reference}.bwt"),
        file("${params.reference}.pac"),
        file("${params.reference}.sa")
    ]).collect()
    
    //
    // Prepare known sites channels
    //
    dbsnp_ch = Channel.of([
        file(params.dbsnp),
        file("${params.dbsnp}.tbi")
    ])
    
    known_indels_ch = Channel.of([
        file(params.known_indels),
        file("${params.known_indels}.tbi")
    ])
    
    //
    // RUN WORKFLOW
    //
    GATK_VARIANT_CALLING (
        ch_input,
        reference_ch,
        bwa_index_ch,
        dbsnp_ch,
        known_indels_ch
    )
}

/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/

workflow.onComplete {
    log.info ( workflow.success ? """
        ==========================================
        Pipeline Completed Successfully!
        ==========================================
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        WorkDir      : ${workflow.workDir}
        Exit status  : ${workflow.exitStatus}
        Results      : ${params.outdir}
        ==========================================
        """ : """
        ==========================================
        Pipeline Failed
        ==========================================
        Failed: ${workflow.errorReport}
        Exit status : ${workflow.exitStatus}
        ==========================================
        """
    )
}
