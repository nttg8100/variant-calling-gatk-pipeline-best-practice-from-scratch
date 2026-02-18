process GATK_APPLYBQSR {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/aligned/${meta.id}", mode: 'copy'
    
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
    
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path reference
    path reference_fai
    path reference_dict
    
    output:
    tuple val(meta), path("*_recal.bam"), emit: bam
    tuple val(meta), path("*_recal.bai"), emit: bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Apply Base Quality Score Recalibration
    gatk ApplyBQSR \\
        -I $bam \\
        -R $reference \\
        --bqsr-recal-file $recal_table \\
        $args \\
        --create-output-bam-index true \\
        -O ${prefix}_recal.bam
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E 'GATK v[0-9.]+' | sed 's/.*v\\([0-9.]\\+\\).*/\\1/' || echo "4.4.0.0")
    END_VERSIONS
    """
}
