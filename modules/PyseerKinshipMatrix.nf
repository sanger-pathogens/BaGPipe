// Pyseer kinship matrix
process PyseerKinshipMatrix {
    publishDir "${params.outdir}/kinship_matrix", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path tree 

    output:
    path "*K.tsv", emit: kinship_matrix

    script:
    """
    phylogeny_distance.py \
        --lmm "${tree}" > phylogeny_K.tsv
    """
}