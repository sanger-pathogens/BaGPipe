process SignificantKmers {
    publishDir "${params.outdir}/significant_unitigs", mode: 'copy', overwrite: true

    input:
    path pyseer_result

    output:
    path "*", emit: sig_kmer_out

    script:
    """
    threshold=\$(grep 'Threshold:' kmer_pattern_count_${params.chosen_phenotype}.txt | cut -f2)
    cat <(head -1 gwas_${params.chosen_phenotype}_kmers.txt) \
        <(awk -v thresh=\$threshold '\$4<thresh {print \$0}' \
        gwas_${params.chosen_phenotype}_kmers.txt) > significant_kmers.txt
    """
}

process KmerMap {
    tag "${fasta}"
    publishDir "${params.outdir}/significant_unitigs/Manhattan", mode: 'copy', overwrite: true
    container "quay.io/sangerpathogens/pyseer:1.3.11"

    input:
    tuple path(fasta), path(gff), path(sig_kmer)

    output:
    path "*.plot", emit: kmer_map_out

    script:
    """
    phandango_mapper "${sig_kmer}" "${fasta}" "mapped_kmers_${fasta}.plot"
    """
}

process WriteReferenceText {
    publishDir "${params.outdir}/significant_unitigs", mode: 'copy', overwrite: true

    input:
    path manifest_ch 
    path ref_manifest_ch

    output:
    path "*", emit: write_ref_text_out

    script:
    output_file = "references.txt"
    """
    while IFS=\$'\\t' read -r col1 col2
    do
        echo -e "\${col1}\\t\${col2}\\tref"
    done < "${ref_manifest_ch}" > "${output_file}"

    while IFS=, read -r sample_id assembly_path
    do
        if [[ \${sample_id} != "sample_id" ]]; then
            echo -e "\${assembly_path}\\t\${sample_id}.gff3\\tdraft"
        fi
    done < "${manifest_ch}" >> "${output_file}"
    """
}

process AnnotateKmers {
    publishDir "${params.outdir}/annotated_unitigs", mode: 'copy', overwrite: true
    container "quay.io/sangerpathogens/pyseer:1.3.11"

    input:
    path sig_kmer
    path reftxt
    path gff_files

    output:
    path "*.tsv", emit: annotated_kmers_out

    script:
    """
    annotate_hits_pyseer "${sig_kmer}" "${reftxt}" annotated_kmers.tsv
    summarise_annotations.py annotated_kmers.tsv > gene_hits.tsv
    """
}

process GeneHitPlot {
    publishDir "${params.outdir}/annotated_unitigs", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/r-ggrepel:0.6.5--r3.3.2_0"

    input:
    path genehit
    path plot_script

    output:
    path "*.pdf", emit: genehit_plot_out

    script:
    """
    Rscript ${plot_script} ${genehit} 
    """
}
