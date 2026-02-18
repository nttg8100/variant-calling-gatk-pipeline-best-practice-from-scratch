#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    GATK Variant Calling Pipeline - Nextflow Version (Part 2)
========================================================================================
    Based on the bash workflow from Part 1
    Github: https://github.com/your-repo/gatk-variant-calling
========================================================================================
*/

// Include modules
include { FASTQC } from './modules/fastqc'

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow GATK_VARIANT_CALLING {
    
    take:
    reads_ch    // channel: [ val(meta), [ path(read1), path(read2) ] ]
    
    main:
    ch_versions = Channel.empty()
    
    //
    // STEP 1: Quality Control with FastQC
    //
    FASTQC (
        reads_ch
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    
    emit:
    fastqc_html = FASTQC.out.html
    fastqc_zip  = FASTQC.out.zip
    versions    = ch_versions
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
    // RUN WORKFLOW
    //
    GATK_VARIANT_CALLING (
        ch_input
    )
}

/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/

workflow.onComplete {
    log.info ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        WorkDir      : ${workflow.workDir}
        Exit status  : ${workflow.exitStatus}
        Results      : ${params.outdir}
        """ : """
        Failed: ${workflow.errorReport}
        Exit status : ${workflow.exitStatus}
        """
    )
}
