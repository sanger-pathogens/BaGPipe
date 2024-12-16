# BaGPipe

## Overview
Running a bacterial GWAS isn't that hard, but the pre- and post-processing of data can be boring. In particular, a user must contend with data conversions, estimate computational resource requirements and wait for one process to finish before another can start. These tasks take precious time away from you, and sometimes bring in small errors that impact the reliability of results.

**BaGPipe** is a nextflow pipeline that integrates a series of bioinformatic tools into a standard, tested workflow for performing bacterial GWAS on large datasets (see Figure below). 

![Flowchart describing the pipeline](images/bgwas_pipeline_0422.drawio.png)

## What does it do?

The most comprehensive way of running BaGPipe starts with only genome assemblies. In this mode, it generates all other data required as input for bacterial GWAS by:
- Creating k-mers (short DNA sequences of length k) or unitigs (longer, non-overlapping sequences assembled from k-mers)
- Annotating assemblies
- Performing a pangenome analysis and creating a core phylogenetic tree
- Generating a pairwise distance matrix 

Users have the flexibility to enter the workflow at alternative starting points in this pre-processing stage, e.g. supplying their own tree.

For the association analysis, [Pyseer](https://github.com/mgalardini/pyseer) (Lees et al, 2018) was used for its speed and design in addressing the common problems of bacterial GWAS like population strucutre, recombination rate, and multiple testing. By default, BaGPipe uses of a linear mixed model algorithm and unitigs as the input genotype (options recommended by Pyseer authors).

In the post-processing stage, BaGPipe automatically performs significant unitig analysis. If the user provides reference files, it can conveniently produce Manhattan plots, annotate the significant unitigs, and eventually produce a gene-hit plot. 

BaGPipe facilitates GWAS analysis on large datasets by easing the computational burden. It optimises the use of requested memory and CPU, which can be customised if necessary. Additionally, an automated resource escalation strategy ensures that a process will be re-run with higher resource requests if the process failed due to lack of memory or runtime on an HPC cluster node.

## Requirements

### Software

- A Unix-like operating system environment (e.g. Linux, macOS, Windows with [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)) with Bash 3.2 or later.
- Java 11 or later ([OpenJDK](https://openjdk.org/)) as a nextflow dependency
- [nextflow](https://www.nextflow.io/) as workflow management system
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) for pipeline dependency management

> NOTE: Java and Nextflow can be installed with instructions on https://www.nextflow.io/docs/latest/install.html.

### Hardware (recommended)

- \>= 16GB RAM
- \>= 2 CPUs/cores
- Enough space for images (< 10GB) and intermediate data

## Getting started

### Example command
Running BaGPipe is easy! You just need to specify a few input files (see [Inputs](#inputs) below). Some example input is provided in the [example_input](./example_input) folder. These example inputs will work with the below commands, provided assemblies available in the [Pyseer tutorial](https://pyseer.readthedocs.io/en/master/tutorial.html) are downloaded and paths in the inputs updated accordingly.

Different configuration profiles can be specified depending on the compute environment in which BaGPipe is run. In most cases, BaGPipe should be run in an HPC environment, but these environments are specific to host institutions. To run BaGPipe with an nf-core profile appropriate for your institution, find a config file [here](https://github.com/nf-core/configs/tree/master/conf). If a config file exists for your institution, run it using `-profile <institution>`. Different profiles can be combined in a list, with later profiles overriding previous ones.

You can also enable different methods for dependency management using profiles. Currently supported profiles for this purpose include: `singularity` and `docker`.

For instance to run on the Sanger HPC, using singularity to handle pipeline dependencies:
```
bsub -q oversubscribed -M 4000 -R "rusage[mem=4000] select[mem>4000]" -o test.o -e test.e \
  nextflow run sanger-pathogens/BaGPipe \
    -profile sanger,singularity \
    --manifest genome_manifest.csv \
    --genotype_method unitig \
    --reference reference_manifest.tsv \
    --genus Streptococcus \
    --phenotypes pheno.tab \
    --antibiotic penicillin \
    -resume
```

If you would like to make code changes it may be better to clone this repository and run the pipeline using `nextflow run main.nf`, instead of pulling directly from GitHub (as in the example command above).

For more options, please explore in the help message (using `--help`).

### Execution and Caching

Users can supply parameters to the pipeline in two ways, via:
- Command Line Interface (CLI), as shown in the [example command](#example-command).
- Nextflow configuration file (similar to [nextflow.config](./nextflow.config))

If a configuration file is specified, the user must establish the paths for all input files and define the specific analysis along with its related parameters. 

In case of any error during a run, the user can rectify the error and restart the pipeline using the `-resume` nextflow option. Doing so utilises cached files from the last run, without the need to re-run everything again from the top.

## Inputs

The user must specify, through the option `--genotype_method`, one of the three variant genotype methods:
1. k-mers/unitigs (`unitig`)
2. gene presence or absence (`pa`)
3. SNPs and INDELs (Insertions and Deletions) (`snp`)

If unsure, I recommend the `unitig` approach, as it is the approach recommended by the Pyseer authors. Currently, the `unitig` option is fully tested but the others are incomplete.

| Input Type        | Format                                                                | Use Case                                                                                        |
|-------------------|-----------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| Genomes           | A CSV file listing all sample IDs (first column) and the paths to their FASTA files (second column) | Compulsory                                                                                      |
| Phenotypes        | A TSV file listing all sample IDs (first column) and whether the isolate belongs to each phenotype (0 or 1, each of the other columns) | Compulsory                                                                                      |
| References        | A TSV file listing the paths to reference assemblies (first column) and the paths to their corresponding GFF files (second column) | Significant unitig analysis                                                                     |
| Annotated genomes | A CSV file listing all sample IDs (first column) and the paths to their GFF files (second column) | If the user prefers their own GFFs in “unitig” mode                                             |
| Phylogenetic tree | A phylogenetic tree file for the pangenome in NEWICK format           | If the user prefers their own tree in “unitig” mode                                             |
| Variant file      | A CSV file listing all sample IDs (first column) and the paths to their VCF files (second column) | If the user prefers their own VCFs in “snp” mode                                                |
| Merged variant file | A merged VCF file for all samples                                     | If the user prefers their own merged VCFs in “snp” mode                                         |


## How it works?
### 1. Annotation
[Prokka](https://github.com/tseemann/prokka) (v.1.14.6) (Seemann, 2014) annotates assembled genomic DNA by employing various tools to predict genomic coordinates and gene functions. Most tools, except Prodigal, provide both coordinates and descriptions; Prodigal solely identifies coding sequence coordinates. Prokka then refines the annotations by comparing the predicted gene coordinates to multiple databases, escalating from smaller, reliable ones to larger, domain-specific collections, ultimately utilising curated protein family models. For each sample, Prokka outputs 10 files with a user-defined "sample_id" prefix. BaGPipe collects all the GFF files into one Nextflow channel for subsequent processing.

### 2. Building multiple sequence alignment of core genes from pangenome
From the annotated assemblies from Prokka or user input, BaGPipe runs [Panaroo](https://github.com/gtonkinhill/panaroo) (v.1.3.4) (Tonkin-Hill et al, 2020) with the option to build a multiple sequence alignment of the core genes using Multiple Alignment using Fast Fourier Transform (MAFFT). Panaroo uses pangenome evolution models instead of gene accumulation curves, and it constructs a graphical pangenome representation, where genes from the nodes and edges link genes that are adjacent on at least one contig. This initial graph is then refined through several cleaning steps to address common issues in genome annotation.

Panaroo outputs several files, including the gene presence/absence format and core and accessory genome alignments, created using the alignment tool MAFFT. The user can select from other alignment tools supported by Panaroo such as Prank or Clustal Omega. The annotated pangenome graph produced by Panaroo is in GML format, compatible with Cytoscape for convenient visualisation. The user can access these outputs and use them in alternative onward analyses, such as the representation of phylogeny and clustering.

### 3. Building phylogeny from the core gene alignment
To control for population structure in the downstream association study, BaGPipe uses [IQ-TREE](https://github.com/iqtree/iqtree2) (v.2.2.6) (Nguyen et al, 2015) to build a phylogeny from the core gene alignment. IQ-TREE is a fast stochastic algorithm which constructs phylogenetic trees using a maximum likelihood method. It employs a fast hill-climbing nearest neighbour interchange (NNI) algorithm to generate initially optimal trees, retaining the best topologies based on likelihood. Following this, a stochastic NNI step perturbs these selected trees, allowing for potential escape from local optima. The algorithm re-applies hill-climbing NNI to perturbed trees, iteratively refining them. This recursive process continues until the algorithm identifies the best tree, which remains unchanged through 100 random perturbations. BaGPipe selectively stores the “.treefile” file into a channel for downstream processes.

IQ-TREE generates multiple output files, including: 
1.	`*.iqtree`: The primary report file which provides detailed computational results and a text-based representation of the final tree. The user can examine this file to understand the outcomes of the run.
2.	`*.treefile`: The maximum likelihood tree in NEWICK format, which can be viewed in a tree viewer program that supports the format, such as FigTree or iTOL. 
3.	`*.log`: A log file which records the entire process of the run

### 4. Association study
BaGPipe uses [Pyseer](https://github.com/mgalardini/pyseer) (v.1.3.11), a non-phylogenetic method and a Python implementation of SEER, to perform GWAS using a linear mixed model. BaGPipe adapts the flexibility given by Pyseer so that the user can choose their preferred methods for doing the association study.

A typical bacterial GWAS requires three types of input: genotype, phenotype, and interaction. One of the options is to study the association between a binary phenotype (i.e. antibiotic resistance) and a binary genotype (i.e. presence/absence) which is created upstream, and the interaction can be the kinship matrix which is produced automatically from the upstream phylogeny in the pipeline. Pyseer can also calculate a pairwise distance matrix that can be used as an interaction from provided genome assemblies using Mash.

From the interaction information it has, Pyseer applies MDS to manage population structure. An alternative, faster method of GWAS supported by Pyseer uses k-mers/unitigs of various lengths to represent genotypic variations is supported by Pyseer. This method has become best-practice and is recommended by the creators of Pyseer.

To control for multiple testing, BaGPipe also counts the number of patterns using scripts from Pyseer and outputs it in a TXT file starting with “pattern_count”. The user can check if this count is reasonably lower than the total number of tested genotypes (k-mers/unitigs/SNPs, etc.), which can aid the decision whether to perform multiple testing corrections. However, this is rarely needed if the k-mer/unitg approach is used.

As a standard protocol, a Q-Q plot is produced after GWAS to assess whether the observed distribution of p-values deviates from the expected distribution under the null hypothesis, i.e., no association between genetic markers and the trait or condition of interest. If the observed p-values deviate systematically from the expected line (which represents the null hypothesis), this could indicate a problem with the study design, population structure, or other confounding factors that may have inflated or deflated the p-values.

If the user selects the k-mer/unitig genotype method (“unitig”), BaGPipe takes steps further in processing those k-mer/unitig outputs, making them more intuitive for the user to analyse the results. First, only the unitigs which exceeded the significance threshold are filtered through. Then, the significant k-mers/unitigs are mapped iteratively to each of the references provided by the user to produce materials for Manhattan plots.

Additionally, BaGPipe automatically annotates the significant k-mers/unitigs. It utilises all the references and the input assemblies as drafts for annotation, then finds annotations for the significant k-mers/unitigs, outputting “gene_hits.tsv”. Annotations labelled “ref” permit partial matches between the k-mers/unitigs and reference sequence, while those marked as “draft” necessitate an exact match. For each significant k-mer/unitig, BaGPipe enriches the output by including the coding genes it resides in, alongside both the closest upstream and downstream coding genes.

While the user can inspect this file directly, there is an example R script provided that can be customised for visualisation. BaGPipe automatically plots these in three different sizes and the user can have a quick look to the overall distribution of p-value and effect size of the annotated genes.

## Output and visualisation
Output files are systematically organised into a hierarchical folder arrangement. The outputs of each process are arranged in independent folders. Additionally, Nextflow generates a separate “work” folder to store intermediate output files, which acts as a specific log for each run.

## Docker and Nextflow implementation of the pipeline
A Docker container has been developed to facilitate the installation of all essential tools and packages, as well as to ensure access to particular versions of these resources. With the Nextflow implementation, the user can run the multi-modular pipeline in one line, specifying which type of input files and GWAS methods they want. BaGPipe can also conveniently be run on an HPC environment supporting the Load Sharing Facility (LSF) executor for swifter processing, as this allows each process to be dispatched as an independent job.

## Customisation
### Configuration file
A single default configuration file is supplied, encompassing all options and parameters for BaGPipe. The user can tailor this file for conducting specific analyses. To simplify the experience for users at all levels of programming knowledge, the configuration file limits modifications to only key parameters. However, experienced users have the option to alter any parameter by changing the source code. Additionally, the user can edit the profile configuration file to run BaGPipe matching the resource settings of their institution’s HPC platform (https://nf-co.re/configs).

## Credits
This pipeline is developed and co-maintained by Charles Wei as part of a Project for Systems Biology at the University of Cambridge. Ewan Harrison, Dinesh Aggarwal, and William Roberts-Sengier supervised the project and provided ideas on pipeline design, implementation and testing.

[Pyseer](https://github.com/mgalardini/pyseer), developed and maintained by [Marco Galardini](https://github.com/mgalardini) and [John Lees](https://github.com/johnlees), along with it's accompanying documentation proved invaluable to the development of the pipeline. We thank them for making this amazing tool available!

Thanks must also go to Beth Blane, who curated and provided the *Staphylococcus aureus* dataset used in testing the pipeline.
