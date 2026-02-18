process GATK_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/bqsr/${meta.id}", mode: 'copy'
    
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict
    path dbsnp
    path dbsnp_tbi
    path known_indels
    path known_indels_tbi
    
    output:
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Base Quality Score Recalibration - Generate table
    gatk BaseRecalibrator \\
        -I $bam \\
        -R $reference \\
        --known-sites $dbsnp \\
        --known-sites $known_indels \\
        $args \\
        -O ${prefix}_recal.table
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+' || echo "4.4.0.0")
    END_VERSIONS
    """
}
