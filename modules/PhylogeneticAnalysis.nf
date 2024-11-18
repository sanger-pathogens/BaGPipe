// Build a phylogeny from the core gene alignment using IQ-TREE
process PhylogeneticAnalysis {
    publishDir "${params.outdir}/iqtree", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    input:
    path alignment

    output:
    path "*", emit: iqtree_out  //TODO Do we really want to output every file in work dir into this channel? Convenient for publishing perhaps, but might be better to be specific.
    path "*.treefile", emit: phylo_tree

    script:
    """
    iqtree \
        -s "${alignment}" \
        -pre core_tree \
        -fast \
        -m ${params.iqtree_model} \
        -nt ${task.cpus}
    """
}