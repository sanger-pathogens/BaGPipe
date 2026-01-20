/*
 * ScreenUnitigs.nf
 *
 * External screening workflow:
 *  1) Build FASTA of significant sequences from significant_kmers.txt (col: 'variant')
 *  2) Call presence/absence across new assemblies using unitig-caller --simple
 *  3) Summarise per-sample / per-unitig, and optionally per-gene
 *
 * Inputs:
 *  - screening_manifest (CSV): sample_id,absolute_path_to_fasta.fa
 *  - significant_kmers (TSV): header with 'variant'
 *  - screening_gene_hits (optional TSV)
 *  - outdir: output directory
 */

workflow screen_unitigs {

    take:
    screening_manifest
    significant_kmers
    genehit
    outdir

    main:

    screening_outdir = "${outdir}/screening"

    /*
     * Input channels
     */
    ch_manifest = screening_manifest
    ch_sigkmers = significant_kmers
    ch_annotated_kmers = genehit


    /*
     * 1) Build FASTA of significant unitigs
     */
    buildSigFasta(
        ch_sigkmers,
        screening_outdir
    )

    ch_sigfa = buildSigFasta.out.sigfa

    /*
     * 2) Presence/absence calling (unitig-caller only)
     */
    unitigCallerSimple(
        ch_manifest,
        ch_sigfa,
        screening_outdir
    )

    ch_calls = unitigCallerSimple.out.calls_rtab

    /*
     * 3) Summaries (+ optional gene roll-up)
     */
    // make the screening_summaries script available inside the task
    script_ch = Channel.value(file("${projectDir}/scripts/screening_summaries.py"))

    screeningSummaries(
        ch_calls,
        ch_annotated_kmers,
        script_ch,
        screening_outdir
    )

    emit:
    per_sample = screeningSummaries.out.per_sample
    per_unitig = screeningSummaries.out.per_unitig
    per_gene   = screeningSummaries.out.per_gene
}


/*
 * Build significant_unitigs.fa from significant_kmers.txt
 */
process buildSigFasta {

    tag "buildSigFasta"

    input:
    path sigkmers_tsv
    val  outdir

    output:
    path "${outdir}/significant_unitigs.fa", emit: sigfa

    container 'python:3.11'

    script:
    """
    mkdir -p ${outdir}

    python - <<'PY'
import sys, csv
from pathlib import Path

inp = Path("${sigkmers_tsv}")
out = Path("${outdir}") / "significant_unitigs.fa"

with inp.open() as f, out.open("w") as g:
    r = csv.DictReader(f, delimiter='\\t')
    if 'variant' not in r.fieldnames:
        sys.exit("ERROR: Input TSV missing 'variant' header")

    seen = set()
    for row in r:
        seq = row['variant'].strip().upper()
        if not seq:
            continue
        if any(c not in 'ACGTN' for c in seq):
            continue
        if seq in seen:
            continue
        seen.add(seq)
        g.write(f">{seq}\\n{seq}\\n")
PY
    """
}


/*
 * unitig-caller --simple
 */
process unitigCallerSimple {

    tag "unitigCallerSimple"

    input:
    path manifest_csv
    path sigfa
    val  outdir

    output:
    path "${outdir}/screening_calls.unitigs.Rtab", emit: calls_rtab
    path "${outdir}/screening_calls.log"

    container 'quay.io/biocontainers/unitig-caller:1.3.1--py311heec5c76_1'

    script:
    """
    set -euo pipefail
    mkdir -p ${outdir}

    # Extract absolute FASTA paths from CSV (column 2), skip header line
    tail -n +2 ${manifest_csv} | cut -d, -f2 > ${outdir}/new_refs.txt

    unitig-caller --simple \
        --refs ${outdir}/new_refs.txt \
        --unitigs ${sigfa} \
        --out ${outdir}/screening_calls \
        2> ${outdir}/screening_calls.log
    
    # Convert pyseer-style output to rtab matrix expected by screening_summaries.py
    python - <<'PY'
from pathlib import Path
import os

py = Path("${outdir}") / "screening_calls.pyseer"
rtab = Path("${outdir}") / "screening_calls.unitigs.Rtab"
if not py.exists():
    raise SystemExit(0)

samples = []
data = {}
with py.open() as f:
    for line in f:
        line = line.strip()
        if not line or '|' not in line:
            continue
        seq, rest = line.split('|', maxsplit=1)
        seq = seq.strip()
        parts = rest.strip().split()
        if seq not in data:
            data[seq] = {}
        for p in parts:
            if ':' not in p:
                continue
            sample, count = p.split(':', 1)
            if sample not in samples:
                samples.append(sample)
            try:
                cnt = int(count)
            except ValueError:
                try:
                    cnt = int(float(count))
                except:
                    cnt = 0
            data[seq][sample] = cnt

if not samples:
    # no samples found; still create empty rtab
    rtab.write_text(os.linesep)
else:
    with rtab.open('w') as g:
        g.write('\t' + '\t'.join(samples) + os.linesep)
        for seq, row in data.items():
            vals = [str(row.get(s, 0)) for s in samples]
            g.write(seq + '\t' + '\t'.join(vals) + os.linesep)
PY
    """
}


/*
 * Summaries + optional gene roll-up
 */
process screeningSummaries {
    tag "screeningSummaries"
    publishDir "${params.outdir}/screening", mode: 'copy', overwrite: true

    input:
    path calls_rtab
    path genehit
    path script
    val  outdir

    output:
    path "${outdir}/screening_per_sample.tsv", emit: per_sample
    path "${outdir}/screening_per_unitig.tsv", emit: per_unitig
    path "${outdir}/screening_per_gene.tsv", emit: per_gene

    container 'python:3.11'

    script:
    """
    set -euo pipefail


    # if AnnotateKmers produced annotated_kmers.tsv it will be staged in the work dir
    python screening_summaries.py \
        --rtab ${calls_rtab} \
        --annotated_kmers annotated_kmers.tsv \
        --outdir ${outdir}
    """
}
