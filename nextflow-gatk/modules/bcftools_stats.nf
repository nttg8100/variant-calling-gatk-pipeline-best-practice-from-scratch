process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/bcftools:1.17--haef29d1_0'
    
    publishDir "${params.outdir}/variants/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), path(raw_vcf), path(raw_tbi), path(filtered_vcf), path(filtered_tbi)
    
    output:
    tuple val(meta), path("*_variant_stats_raw.txt"), emit: raw_stats
    tuple val(meta), path("*_variant_stats_filtered.txt"), emit: filtered_stats
    tuple val(meta), path("*_raw_snp_count.txt"), emit: raw_snp_count
    tuple val(meta), path("*_raw_indel_count.txt"), emit: raw_indel_count
    tuple val(meta), path("*_filtered_snp_count.txt"), emit: filtered_snp_count
    tuple val(meta), path("*_filtered_indel_count.txt"), emit: filtered_indel_count
    path("versions.yml"), emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    # Generate bcftools stats for raw variants
    bcftools stats ${raw_vcf} > ${prefix}_variant_stats_raw.txt
    
    # Generate bcftools stats for filtered variants
    bcftools stats ${filtered_vcf} > ${prefix}_variant_stats_filtered.txt
    
    # Count raw variants by type
    bcftools view -v snps ${raw_vcf} | bcftools query -f '.\\n' | wc -l > ${prefix}_raw_snp_count.txt
    bcftools view -v indels ${raw_vcf} | bcftools query -f '.\\n' | wc -l > ${prefix}_raw_indel_count.txt
    
    # Count filtered variants by type (PASS only)
    bcftools view -f PASS -v snps ${filtered_vcf} | bcftools query -f '.\\n' | wc -l > ${prefix}_filtered_snp_count.txt
    bcftools view -f PASS -v indels ${filtered_vcf} | bcftools query -f '.\\n' | wc -l > ${prefix}_filtered_indel_count.txt
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
