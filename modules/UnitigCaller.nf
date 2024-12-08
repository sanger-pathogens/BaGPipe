process UnitigCaller {
    publishDir "${params.outdir}/unitigs", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/unitig-counter:1.1.0--h56fc30b_0"

    input:
    path manifest

    output:
    path "*", emit: unitig_all
    path "output/*.txt.gz", emit: unitig_out

    script:
    """
    cat ${manifest} | tr ',' '\t' > assemblies_list.txt
    unitig-counter \\
        -strains assemblies_list.txt \\
        -output output \\
        -nb-cores ${task.cpus}
    gzip output/unitigs.txt
    """
}