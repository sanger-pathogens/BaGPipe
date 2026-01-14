process AMRFinderPlus {
    tag "${sample_id}"
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.23--hf69ffd2_0'

    publishDir "${params.outdir}/amr_finder", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(assembly_path), path(amrfinder_db)

    output:
    tuple val(sample_id), path("${sample_id}_amrfinder.tsv"), emit: amrfinder_output

    script:
    amended_id = "${sample_id}".replaceAll(/[^\w.-]/, '_')
    out = "${amended_id}_amrfinder.tsv"
    // Append /latest if the directory contains versioned subdirectories
    db_path = (amrfinder_db.toString() =~ /latest$/) ? amrfinder_db : amrfinder_db + "/latest"
    """
    amrfinder -n ${assembly_path} -o ${out} --threads ${task.cpus} --database ${db_path} ${params.amrfinder_args ?: ''}
    """
}

process AMRFinderSummary {
    publishDir "${params.outdir}/amr_finder", mode: 'copy', overwrite: true

    input:
    path amr_files

    output:
    path "amr_summary.tsv", emit: amrfinder_summary

    script:
    """
    # Concatenate all AMRFinder TSVs (skip errors if none)
    cat ${amr_files} > combined_amrfinder.tsv || true

    python3 - ${amr_files} <<'PY'
import sys
from collections import Counter

files = sys.argv[1:]
rows = []
header = None
idx = None

for f in files:
    try:
        with open(f) as fh:
            h = fh.readline().rstrip("\\n")
            if header is None:
                header = h
                cols = header.split("\\t")
                # Look for the Element symbol column
                if "Element symbol" in cols:
                    idx = cols.index("Element symbol")
                elif "gene_symbol" in cols:
                    idx = cols.index("gene_symbol")
                else:
                    # Default to column 5 (0-indexed)
                    idx = 5
            
            for line in fh:
                if line.strip():
                    rows.append(line.rstrip("\\n").split("\\t"))
    except Exception as e:
        pass

cnt = Counter()
for r in rows:
    if len(r) > idx:
        gene = r[idx].strip()
        if gene and gene != "NA":
            cnt[gene] += 1

with open("amr_summary.tsv","w") as out:
    out.write("Element symbol\\tcount\\n")
    for gene, c in cnt.most_common():
        out.write(f"{gene}\\t{c}\\n")
PY
    """
}
