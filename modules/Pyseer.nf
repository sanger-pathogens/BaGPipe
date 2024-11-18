process PyseerUnitig {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path unitig
    path pheno
    path k_matrix

    output:
    path "*", emit: pyseer_out

    script:
    """
    pyseer \
        --lmm \
        --phenotypes "${pheno}" \
        --kmers "${unitig}" \
        --similarity "${k_matrix}" \
        --phenotype-column "${params.antibiotic}" \
        --output-patterns "kmer_patterns_${params.antibiotic}.txt" \
        --min-af ${params.pyseer_min_af} \
        --max-af ${params.pyseer_max_af} \
        --cpu ${task.cpus} > "gwas_${params.antibiotic}_kmers.txt"
    count_patterns.py "kmer_patterns_${params.antibiotic}.txt" > "kmer_pattern_count_${params.antibiotic}.txt"
    qq_plot.py "gwas_${params.antibiotic}_kmers.txt"
    """
}

process PyseerVariants {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path variants
    path pheno
    path k_matrix

    output:
    path "*", emit: pyseer_out

    script:
    """
    pyseer \
        --lmm \
        --phenotypes "${pheno}" \
        --vcf "${variants}" \
        --similarity "${k_matrix}" \
        --phenotype-column "${params.antibiotic}" \
        --output-patterns "var_patterns_${params.antibiotic}.txt" \
        --min-af ${params.pyseer_min_af} \
        --max-af ${params.pyseer_max_af} \
        --cpu ${task.cpus} > "gwas_${params.antibiotic}_var.txt"
    count_patterns.py "var_patterns_${params.antibiotic}.txt" > "var_pattern_count_${params.antibiotic}.txt"
    qq_plot.py "gwas_${params.antibiotic}_var.txt"
    """
}

process PyseerPreAbs {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path pre_abs
    path pheno
    path k_matrix

    output:
    path "*", emit: pyseer_out

    script:
    """
    pyseer \
        --lmm \
        --phenotypes "${pheno}" \
        --pres "${pre_abs}" \
        --similarity "${k_matrix}" \
        --phenotype-column ${params.antibiotic} \
        --output-patterns "gene_patterns_${params.antibiotic}.txt" \
        --min-af ${params.pyseer_min_af} \
        --max-af ${params.pyseer_max_af} \
        --cpu ${task.cpus} > "gwas_${params.antibiotic}_preabs.txt"
    count_patterns.py "gene_patterns_${params.antibiotic}.txt" > "gene_pattern_count_${params.antibiotic}.txt"
    qq_plot.py "gwas_${params.antibiotic}_preabs.txt"
    """
}
