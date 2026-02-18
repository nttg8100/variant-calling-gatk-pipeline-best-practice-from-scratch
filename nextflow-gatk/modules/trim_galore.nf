process TRIM_GALORE {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/trimmed/${meta.id}", mode: 'copy'
    
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_val_{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("*_trimming_report.txt"), emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: '--quality 20 --length 50 --fastqc'
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Adapter Trimming and Quality Filtering
    trim_galore \\
        --paired \\
        $args \\
        --cores $task.cpus \\
        ${reads[0]} ${reads[1]}
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trim_galore: \$(trim_galore --version | grep -oP 'version \\K[0-9.]+')
    END_VERSIONS
    """
}
