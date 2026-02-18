process GATK_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_high'
    
    publishDir "${params.outdir}/aligned/${meta.id}", mode: 'copy'
    
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*_dedup.bam"), emit: bam
    tuple val(meta), path("*_dedup.bai"), emit: bai
    tuple val(meta), path("*_metrics.txt"), emit: metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Mark Duplicates with GATK
    gatk MarkDuplicates \\
        -I $bam \\
        -O ${prefix}_dedup.bam \\
        -M ${prefix}_metrics.txt \\
        $args \\
        --CREATE_INDEX true
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E 'GATK v[0-9.]+' | sed 's/.*v\\([0-9.]\\+\\).*/\\1/' || echo "4.4.0.0")
    END_VERSIONS
    """
}
