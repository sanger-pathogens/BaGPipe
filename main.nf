#!/usr/bin/env nextflow

log.info """\
    BAGPIPE PIPELINE PARAMETERS
    ==========================
    genus: ${params.genus}
    phenotypes: ${params.phenotypes}
    chosen_phenotype: ${params.chosen_phenotype}
    genotype_method: ${params.genotype_method}
    """
    .stripIndent()

/*
========================================================================================
    HELP
========================================================================================
*/

def printHelp() {
    log.info """
    Usage:
    nextflow run main.nf

    Options:

      --manifest                   Manifest containing paths to FASTA files (mandatory)
      --phenotypes                 A tab file containing phenotypes for all samples (mandatory)
      --chosen_phenotype           Phenotype to use for the GWAS analysis (mandatory)
      --genus                      Genus name for samples (mandatory if starting from FASTA files)
      --genotype_method            Genotype method to run GWAS, from a choice of three (unitig|pa|snp) (mandatory)
                                   Note: unitig is recommended.
      --annotation_method          Tool to use for annotation (prokka|bakta). Default: bakta. (optional)
      --reference                  Manifest containing paths to reference FASTA and GFF files (mandatory for significant k-mer/unitig analysis)
      --mygff                      Input a manifest containing paths to already annotated GFF files; must match sample_ids in manifest (optional)
      --mytree                     Input user preferred phylogenetic tree (optional)
      --mvcf                       Input already mergerd vcf.gz file (optional)
      --fe                         Run GWAS using fixed model (SEER) (optional)
      --help                       Print this help message (optional)

    Alternative Options for Some Processes:

    [PanarooAnalysis]
      --panaroo_clean_mode       Default: "strict"
      --panaroo_alignment        Default: "core"
      --panaroo_aligner          Default: "mafft"
      --panaroo_core_threshold   Default: 0.95
    [PhylogeneticAnalysis]
      --iqtree_model             Evolutionary model that will be used by IQTree. Default: "GTR"
      --iqtree_bootstrap_trees   Number of bootstrap trees to use. Default: 1000
      --iqtree_seed              Random seed to ensure pipeline runs are deterministic. Default: 1234
      --iqtree_args              String supplying additional IQtree options. Incompatible with the following (reserved for use by this pipeline): -pre, -fconst, -s, -nt, -ntmax, -mem, -bb, -m, -seed. Default: ""
    [Pyseer]
      --pyseer_min_af            Default: 0.01
      --pyseer_max_af            Default: 0.99

    """.stripIndent()
}


/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { validate_parameters } from './modules/helper_functions'
include { ProkkaAnnotate } from './modules/ProkkaAnnotate'
include { BaktaAnnotate } from './modules/BaktaAnnotate'
include { PanarooAnalysis } from './modules/PanarooAnalysis'
include { PhylogeneticAnalysis } from './modules/PhylogeneticAnalysis'
include { PyseerKinshipMatrix } from './modules/PyseerKinshipMatrix'
include { UnitigCaller } from './modules/UnitigCaller'
include { MergeVCF } from './modules/MergeVCF'
include { PyseerGenotypeMatrix } from './modules/PyseerGenotypeMatrix'
include {
    PyseerUnitig;
    PyseerPreAbs;
    PyseerVariants
} from './modules/Pyseer'
include {
    SignificantKmers;
    KmerMap;
    ManhattanPlotR;
    WriteReferenceText;
    AnnotateKmers;
    GeneHitPlot
} from './modules/SignificantKmerAnalysis'
include { 
    screen_unitigs;
    buildSigFasta;
    unitigCallerSimple;
    screeningSummaries
} from './modules/ScreenUnitigs'
include {
    AMRFinderPlus;
    AMRFinderSummary
} from './modules/AMRFinderPlus'

/*
========================================================================================
    PARAMETER DEFAULTS
========================================================================================
*/

if (!params.containsKey('amrfinder_args')) {
    params.amrfinder_args = ''
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    // Print help message if needed
    if (params.help) {
        printHelp()
        exit(0)
    }

    // Validate inputs
    validate_parameters()

    // Parse assembly manifest
    manifest_ch = Channel.fromPath(params.manifest)

    genomes_ch = manifest_ch.splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_id, row.assembly_path) }

    // Get phenotypes file
    pheno = Channel.fromPath(params.phenotypes)

    // Generate annotations if necessary
    if (params.mygff) {
        gff_files = Channel.fromPath(params.mygff)
            .splitCsv(header: true, sep: ',')
            .map { row -> tuple(row.sample_id, row.ann_genome_path) }
    } else {
        if (params.annotation_method == "prokka") {
            ProkkaAnnotate(genomes_ch)
            gff_files = ProkkaAnnotate.out.prokka_output_gff
        } else if (params.annotation_method == "bakta") {
            bakta_db = Channel.fromPath(params.bakta_db, checkIfExists: true)
            bakta_input = genomes_ch.combine(bakta_db)
            BaktaAnnotate(bakta_input)
            gff_files = BaktaAnnotate.out.bakta_output_gff
        } else {
            throw new Exception("Unknown annotator passed to `--annotation_method`. Please choose from: bakta, prokka")
        }
    }
    gff_files = gff_files
        .map { it -> it[1] }
        .collect()

    // Optional: run AMRFinderPlus on assemblies and summarise results
    if (params.run_amrfinder) {
        if (params.amrfinder_db) {
            amr_db = Channel.fromPath(params.amrfinder_db, checkIfExists: true)
        } else {
            amr_db = Channel.value(null)
        }

        amr_input = genomes_ch.combine(amr_db)
        AMRFinderPlus(amr_input)
        amr_files = AMRFinderPlus.out.amrfinder_output
            .map { it -> it[1] }
            .collect()
        AMRFinderSummary(amr_files)
    }

    // Generate tree if necessary
    if (params.mytree) {
        tree = Channel.fromPath(params.mytree)
    } else {
        PanarooAnalysis(gff_files)
        alignment = PanarooAnalysis.out.panaroo_output_core_aln

        PhylogeneticAnalysis(alignment)
        tree = PhylogeneticAnalysis.out.phylo_tree
    }

    // Presence/Absence approach
    if (params.genotype_method == "pa") {
        PyseerKinshipMatrix(tree)
        k_matrix = PyseerKinshipMatrix.out.kinship_matrix
        pre_abs = PanarooAnalysis.out.panaroo_output_pre_abs

        PyseerPreAbs(pre_abs, pheno, k_matrix)
    }
    
    // SNP-based approach
    if (params.genotype_method == "snp") {
        if (params.mvcf) {
            // But currently this has a problem: Must use distance matrix with fixed effects (SEER)
            // This is the classifical way of using MDS
            // Distance matrix needs mash (1); or we could extract patristic distances from a phylogeny (2).
                // 1>> mash sketch -s 10000 -o samples *.fa \ mash dist samples.msh samples.msh | square_mash > mash.tsv
                // 2>> python scripts/phylogeny_distance.py core_genome.tree > phylogeny_distances.tsv
            // But, we can use core_gene VCF to generate kinship matrix and use FaST-LMM model?
                // similarity_pyseer --vcf core_gene_snps.vcf sample_list.txt > genotype_kinship.tsv
            if (params.fe) {
                // fixed effect way of doing GWAS
            } else {
                variants = Channel.fromPath(params.mvcf)
                PyseerGenotypeMatrix(variants, manifest_ch)
                k_matrix = PyseerGenotypeMatrix.out.kinship_matrix
                PyseerVariants(variants, pheno, k_matrix)
            }
        } else {
            // Either use snippy to call variant from fasta, given a reference, then use Process: MergeVCF
            // Or ask user to input another manifest containing directory of all vcf, then use Process: MergeVCF
        }
    }
    
    // Unitig approach
    if (params.genotype_method == "unitig") {
        PyseerKinshipMatrix(tree)
        k_matrix = PyseerKinshipMatrix.out.kinship_matrix

        UnitigCaller(manifest_ch)
        unitig = UnitigCaller.out.unitig_out

        PyseerUnitig(unitig, pheno, k_matrix)
        pyseer_result = PyseerUnitig.out.pyseer_out

        SignificantKmers(pyseer_result)
        sig_kmer = SignificantKmers.out.sig_kmer_out
    }

    // Annotate reference with significant kmers
    if (params.reference) {
        ref_manifest_ch = Channel.fromPath(params.reference)

        ref_ch = ref_manifest_ch.splitCsv(header: false, sep:"\t")
            .map { row -> tuple(row[0], row[1]) }
            .combine(sig_kmer)

        // Extract reference GFF path for Manhattan plot annotation
        ref_gff_ch = ref_manifest_ch.splitCsv(header: false, sep:"\t")
            .map { row -> file(row[1]) }
            .first()

        // Preserve Phandango-compatible plots for interactive visualisation
        KmerMap(ref_ch)

        WriteReferenceText(manifest_ch, ref_manifest_ch)
        reftxt = WriteReferenceText.out.write_ref_text_out

        AnnotateKmers(sig_kmer, reftxt, gff_files)
        genehit = AnnotateKmers.out.annotated_kmers_out

        // Flatten potential list emissions and keep only annotated_kmers.tsv for Manhattan plots
        annot_kmers = genehit.flatMap { item ->
            (item instanceof List) ? item : [item]
        }.filter { f -> f.name == 'annotated_kmers.tsv' || f.toString().endsWith('annotated_kmers.tsv') }

        // Create Manhattan plots from annotated kmers (which have genomic position data)
        manhattan_script_ch = Channel.value(file("${projectDir}/scripts/manhattan_kmers.R"))

        threshold_file_ch = pyseer_result
            .flatMap { item -> (item instanceof List) ? item : [item] }
            .filter { f -> f.name == "kmer_pattern_count_${params.chosen_phenotype}.txt" || f.toString().endsWith("kmer_pattern_count_${params.chosen_phenotype}.txt") }
            .view { f -> "Using threshold file: ${f}" }

        annot_with_threshold = annot_kmers.combine(threshold_file_ch)
        
        ManhattanPlotR(annot_with_threshold, manhattan_script_ch, ref_gff_ch)

        plot_script = Channel.value(file("${projectDir}/scripts/gene_hit_summary_plot.R"))
        GeneHitPlot(genehit, plot_script)

        // Screen unitigs if requested
        if (params.screen_unitigs) {
            screening_manifest_ch = Channel.fromPath(params.screening_manifest)

            screen_unitigs(
                screening_manifest_ch,
                sig_kmer,
                genehit,
                params.outdir
            )
        }
    }

}
