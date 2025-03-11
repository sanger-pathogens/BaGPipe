// Build a pangenome, including a multiple sequence alignment of core genes (MAFFT), using Panaroo 
process PanarooAnalysis {
    publishDir "${params.outdir}/panaroo", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/panaroo:1.5.2--pyhdfd78af_0"

    input:
    path gff_files, stageAs: 'gffs/*'

    output:
    path "panaroo_output", emit: panaroo_out
    path "panaroo_output/core_gene_alignment.aln", emit: panaroo_output_core_aln
    path "panaroo_output/gene_presence_absence.Rtab", emit: panaroo_output_pre_abs

    script:
    """
    panaroo \
        -i gffs/* \
        -o panaroo_output \
        --clean-mode ${params.panaroo_clean_mode} \
        --remove-invalid-genes \
        -a ${params.panaroo_alignment} \
        --aligner ${params.panaroo_aligner} \
        --core_threshold ${params.panaroo_core_threshold} \
        -t ${task.cpus}
    """
}