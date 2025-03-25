process PyseerUnitig {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/sangerpathogens/pyseer:1.3.11"

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
        --phenotype-column "${params.chosen_phenotype}" \
        --output-patterns "kmer_patterns_${params.chosen_phenotype}.txt" \
        --min-af ${params.pyseer_min_af} \
        --max-af ${params.pyseer_max_af} \
        --cpu ${task.cpus} > "gwas_${params.chosen_phenotype}_kmers.txt"
    count_patterns.py "kmer_patterns_${params.chosen_phenotype}.txt" > "kmer_pattern_count_${params.chosen_phenotype}.txt"
    qq_plot.py "gwas_${params.chosen_phenotype}_kmers.txt"
    """
}

process PyseerVariants {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/sangerpathogens/pyseer:1.3.11"

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
        --phenotype-column "${params.chosen_phenotype}" \
        --output-patterns "var_patterns_${params.chosen_phenotype}.txt" \
        --min-af ${params.pyseer_min_af} \
        --max-af ${params.pyseer_max_af} \
        --cpu ${task.cpus} > "gwas_${params.chosen_phenotype}_var.txt"
    count_patterns.py "var_patterns_${params.chosen_phenotype}.txt" > "var_pattern_count_${params.chosen_phenotype}.txt"
    qq_plot.py "gwas_${params.chosen_phenotype}_var.txt"
    """
}

process PyseerPreAbs {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/sangerpathogens/pyseer:1.3.11"

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
        --phenotype-column ${params.chosen_phenotype} \
        --output-patterns "gene_patterns_${params.chosen_phenotype}.txt" \
        --min-af ${params.pyseer_min_af} \
        --max-af ${params.pyseer_max_af} \
        --cpu ${task.cpus} > "gwas_${params.chosen_phenotype}_preabs.txt"
    count_patterns.py "gene_patterns_${params.chosen_phenotype}.txt" > "gene_pattern_count_${params.chosen_phenotype}.txt"
    qq_plot.py "gwas_${params.chosen_phenotype}_preabs.txt"
    """
}
