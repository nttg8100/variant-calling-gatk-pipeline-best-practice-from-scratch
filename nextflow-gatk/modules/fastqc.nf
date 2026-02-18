process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/qc/${meta.id}", mode: 'copy'
    
    // Use BioContainers FastQC image
    // Works with both Docker and Singularity (no root required)
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path "versions.yml"            , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Determine if single-end or paired-end
    def input_files = reads instanceof List ? reads.join(' ') : reads
    
    """
    # Quality Control with FastQC
    # This generates HTML and ZIP files with quality metrics
    fastqc \\
        $args \\
        --threads $task.cpus \\
        $input_files
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
