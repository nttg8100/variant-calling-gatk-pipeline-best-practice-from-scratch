process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    
    publishDir "${params.outdir}/aligned/${meta.id}", mode: 'copy', pattern: "*.bam"
    
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'
    
    input:
    tuple val(meta), path(reads)
    path reference
    path reference_fai
    path reference_dict
    path bwa_index
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: '-M'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}_lib"
    
    """
    # Read Alignment with BWA-MEM
    bwa mem \\
        -t $task.cpus \\
        $args \\
        -R "${read_group}" \\
        $reference \\
        ${reads[0]} ${reads[1]} \\
        | samtools view -Sb - > ${prefix}_aligned.bam
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -oP 'Version: \\K[0-9.]+' || echo "unknown")
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
