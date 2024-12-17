// Annotate fasta genome assemblies using Prokka
process BaktaAnnotate {
    tag "${sample_id}"
    
    container "quay.io/biocontainers/bakta:1.9.4--pyhdfd78af_0"

    publishDir "${params.outdir}/bakta", mode:'copy', overwrite: true

    input:
    tuple val(sample_id), path(assembly_path), path(bakta_db)

    output:
    tuple val(sample_id), path("${sample_id}"), emit: bakta_output
    tuple val(sample_id), path("${sample_id}/*.gff3"), emit: bakta_output_gff

    script:
    amended_id = "${sample_id}".replaceAll(/[^\w.-]/, '_')
    gff = "${amended_id}.gff3"
    """
    bakta \\
        ${assembly_path} \\
        ${params.bakta_args} \\
        --threads ${task.cpus} \\
        --output ${sample_id} \\
        --prefix ${amended_id} \\
        --locus-tag ${amended_id} \\
        --db ${bakta_db} \\
        --keep-contig-headers

    # Remove non-ASCII characters from GFF
    sed -i "s/[‘’]/'/g" "${sample_id}/${gff}"
    """
}
