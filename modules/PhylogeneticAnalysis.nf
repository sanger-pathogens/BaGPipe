// Extract SNPs from core gene alignment
process ExtractSNPs {
    publishDir "${params.outdir}/snp-sites", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/snp-sites:2.5.1--he4a0461_5"

    input:
    path(panaroo_core_aln)

    output:
    path(core_alignment_snps), emit: core_alignment_snps

    script:
    core_alignment_snps = "${panaroo_core_aln.baseName}_snps.aln"
    """
    snp-sites \\
        -c ${panaroo_core_aln} \\
        -o ${core_alignment_snps}
    """
}

// Extract SNPs from core gene alignment
process ConstantSitesFreq {
    publishDir "${params.outdir}/snp-sites", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/snp-sites:2.5.1--he4a0461_5"

    input:
    path(panaroo_core_aln)

    output:
    path(core_alignment_constant_sites), emit: constant_sites_freq

    script:
    core_alignment_constant_sites = "${panaroo_core_aln.baseName}.conscount"
    """
    snp-sites \\
        -C ${panaroo_core_aln} \\
    > ${core_alignment_constant_sites}
    """
}

// Build a phylogeny from the core gene alignment using IQ-TREE
process IQTree {
    publishDir "${params.outdir}/iqtree", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    input:
    tuple path(snp_aln), path(constant_sites_freq)

    output:
    path "*", emit: iqtree_out  //TODO Do we really want to output every file in work dir into this channel? Convenient for publishing perhaps, but might be better to be specific.
    path "*.treefile", emit: phylo_tree

    script:
    """
    iqtree \\
        -pre core_tree \\
        -fconst \$(head -n 1 ${constant_sites_freq}) \\
        -s ${snp_aln} \\
        -nt AUTO \\
        -ntmax ${task.cpus} \\
        -bb ${params.iqtree_bootstrap_trees} \\
        -m ${params.iqtree_model} \\
        -seed ${params.iqtree_seed} \\
        ${params.iqtree_args}
    """
}

workflow PhylogeneticAnalysis {
    take:
    alignment

    main:
    ExtractSNPs(alignment)
    ConstantSitesFreq(alignment)

    ExtractSNPs.out.core_alignment_snps
    | combine(ConstantSitesFreq.out.constant_sites_freq)
    | set { phylo_input }

    IQTree(phylo_input)

    emit:
    snp_aln = ExtractSNPs.out.core_alignment_snps
    phylo_tree = IQTree.out.phylo_tree
}