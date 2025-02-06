// Pyseer kinship matrix (genotype/design matrix) from variant presence absence
process PyseerGenotypeMatrix {
    publishDir "${params.outdir}/kinship_matrix", mode: 'copy', overwrite: true
    container "quay.io/sangerpathogens/pyseer:1.3.11"

    input:
    path variants 
    path manifest_ch 

    output:
    path "*K.tsv", emit: kinship_matrix

    script:
    """
    similarity_pyseer --version > version.txt
    tail -n +2 "${manifest_ch}" | cut -d',' -f1 > sample_list.txt
    similarity_pyseer --vcf "${variants}" sample_list.txt > genotype_K.tsv
    """
}