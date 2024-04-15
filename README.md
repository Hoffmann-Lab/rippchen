# rippchen

...are tasty! Acquire a taste for peak calling from *IP-seq experiments or for expression, gene-fusion, methylation and gene set analyses.

Rippchen leverages on bashbone, which is a bash/biobash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses. Rippchen makes use of bashbones best-practice parameterized and run-time tweaked software wrappers and compiles them into a multi-threaded pipeline to analyze model and non-model organisms.

# Outline

- [Rippchen](#rippchen)
  - [Features](#features)
  - [Covered Tasks](#covered-tasks)
- [License](#license)
- [Download](#download)
- [Installation](#installation)
  - [Full installation of all third party tools used in bashbones biobash library](#full-installation-of-all-third-party-tools-used-in-bashbones-biobash-library)
  - [Upgrade to a newer release or extend an existing bashbone installation](#upgrade-to-a-newer-release-or-extend-an-existing-bashbone-installation)
  - [Update tools](#update-tools)
- [Usage](#usage)
  - [Retrieve SRA datasets](#retrieve-sra-datasets)
  - [Retrieve genomes](#retrieve-genomes)
  - [Merge/collate features](#mergecollate-features)
  - [Generate pseudo replicates](#generate-pseudo-replicates)
  - [Sample info file](#sample-info-file)
  - [Adapter sequences](#adapter-sequences)
  - [Examples](#examples)
    - [Index genomes](#index-genomes)
    - [FastQ data processing](#fastq-data-processing)
    - [Alignment data processing](#alignment-data-processing)
    - [Multiple inputs](#multiple-inputs)
      - [Differential analysis](#differential-analysis)
      - [Peak calling](#peak-calling)
    - [Start, redo or resume](#start-redo-or-resume)
- [Outputs](#outputs)
  - [counts](#counts)
  - [coexpressed](#coexpressed)
    - [Cluster enrichment test directories](#cluster-enrichment-test-directories)
  - [deseq](#deseq)
    - [Enrichment test directories](#enrichment-test-directories)
  - [diego](#diego)
  - [metilene](#metilene)
  - [stats](#stats)
  - [mapped](#mapped)
  - [mecall](#mecall)
  - [peaks](#peaks)
  - [qualities](#qualities)
  - [Pre-processing directories](#pre-processing-directories)
- [Third-party software](#third-party-software)
- [Supplementary information](#supplementary-information)
- [Closing remarks](#closing-remarks)

## Features
[&#x25B2; back to top](#rippchen)

- Most software related parameters will be inferred directly from your data so that all functions require just a minimalistic set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)
- Extensive logging and different verbosity levels
- Start, redo or resume from any point within your data analysis pipeline

## Covered Tasks
[&#x25B2; back to top](#rippchen)

- For paired-end and single-end derived raw sequencing or prior mapped read data
  - RNA-Seq protocols
  - DNA-Seq protocols
  - Bisulfite converted DNA-Seq protocols
- Data preprocessing (quality check, adapter clipping, quality trimming, error correction, artificial rRNA depletion)
- Read alignment and post-processing
  - knapsack problem based slicing of alignment files for parallel task execution
  - sorting, filtering, unique alignment extraction, removal of optical duplicates
  - generation of pools and pseudo-replicates
- Gene fusion detection
- Methyl-C calling and prediction of differentially methylated regions
- Expression analysis
  - Read quantification, TPM and Z-score normalization and heatmap plotting
  - Inference of strand specific library preparation methods
  - Inference of differential expression as well as co-expression clusters
  - Detection of differential splice junctions and differential exon usage
  - Gene ontology (GO) gene set enrichment and over representation analysis plus semantic similarity based clustering
- Implementation of Encode3 best-practice ChIP-Seq Peak calling 
  - Peak calling from RIP-Seq, MeRIP-Seq, m6A-Seq and other related *IP-Seq data
  - Inference of effective genome sizes

# License
[&#x25B2; back to top](#rippchen)

The whole project is licensed under the GPL v3 (see LICENSE file for details). <br>
**except** the the third-party tools set-upped during installation. Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download
[&#x25B2; back to top](#rippchen)

This will download you a copy which includes the latest developments

```bash
git clone --recursive https://github.com/Hoffmann-Lab/rippchen
```

To check out the latest release (**irregularly compiled**) do

```bash
cd rippchen
git checkout --recurse-submodules $(git describe --tags)
```

# Installation
[&#x25B2; back to top](#rippchen)

```bash
scripts/setup.sh -h
```

## Full installation of all third party tools used in bashbones biobash library
[&#x25B2; back to top](#rippchen)

When using the `-g` switch (**recommended**), the setup routine will create conda environments or setups software from source according to enclosed configuration files, URLs respectively. Without `-g` switch, software is installed in latest available version, which may lead to unexpected behavior and errors. During setup, current configuration files will be written to `<path/of/installation/config>`.

```bash
scripts/setup.sh -g -i all -d <path/to/installation>
source <path/of/installation>/latest/rippchen/activate.sh
bashbone -h
```

## Upgrade to a newer release or extend an existing bashbone installation
[&#x25B2; back to top](#rippchen)

Use the `-g` switch, in order to also upgrade conda environments that fail the comparison with the supplied configuration files. **Attention**: This switch will downgrade tools, if the initial installation was done for cutting edge tools i.e. without `-g`.

```bash
scripts/setup.sh -g -i upgrade -d <path/of/installation>
```

## Update tools
[&#x25B2; back to top](#rippchen)

Trimmomatic, segemehl, STAR-Fusion, GEM, mdless and gztool will be installed next to the conda environments. Their latest versions and download URLs will be automatically inferred.

```bash
scripts/setup.sh -i trimmomatic,segemehl,starfusion,gem,mdless,gztool -d <path/of/installation>
```

# Usage
[&#x25B2; back to top](#rippchen)

To load rippchen, bashbone respectively, execute

```
source <path/of/installation/latest/rippchen/activate.sh> -c true
rippchen.sh -h

BASIC OPTIONS
-h       | --help                     : prints this message
-dev     | --devel                    : prints list of keywords in processing order for advanced pipeline control
-e       | --env                      : list tools and versions in setupped environment
-v       | --verbosity [value]        : set level of verbosity. default: 0
                                        0 - get simple status updates
                                        1 - get status updates and commands
                                        2 - get full output
-o       | --out [path]               : output directory. default: $OUTDIR
-l       | --log [path]               : output directory. default: $OUTDIR/run.log
-tmp     | --tmp                      : temporary directory. default: ${TMPDIR-/tmp}/rippchen.XXXXXXXXXX
                                        NOTE: respects TMPDIR environment variable
-k       | --keep                     : keep temporary and unnecessary files
-t       | --threads [value]          : number of threads. default: $THREADS
-xmem    | --max-memory [value]       : fraction or total amount of allocatable memory in MB. default: $MAXMEMORY MB i.e. currently available memory
-mem     | --memory [value]           : allocatable memory per instance of memory greedy tools in MB. defines internal number of parallel instances
                                        default: $MEMORY which allows for $MTHREADS instances and $MTHREADS SAM/BAM slices according to -xmem
                                        NOTE: needs to be raised in case of GCThreads, HeapSize or OutOfMemory errors
-resume  | --resume-from [string]     : resume from a specific pipeline step (see -dev)
                                        NOTE: define entry point before skip and redo
-skip    | --skip [string,..]         : skip specific pipeline step(s). comma separated (see -dev)
-redo    | --redo [string,..]         : redo specific pipeline step(s). comma separated (see -dev)


INDEXING OPTIONS
-x       | --index                    : triggers creation of all requiered genome and annotation indices plus md5 sums
-b       | --bisulfite [string]       : triggers indices for for methylation analyses. use keyword WGBS
-g       | --genome [path]            : genome fasta input
-gtf     | --gtf [path]               : annotation gtf input
                                        NOTE: no gtf file implies star index creation without splice junctions database
-it      | --is-transcriptome         : fasta and gtf input is actually a transcriptome converted. uses -no-dsj
                                        NOTE: requires transcript_id as sequence ids and tags in gtf to sum up fractional counts. see genome2transcriptome.pl
-go      | --go [path]                : annotation go input
                                        NOTE: no go file implies disabled creation of org.db from semantically clustered gene ontology terms
-no-sege | --no-segemehl              : disables indexing for segemehl
-no-star | --no-star                  : disables indexing for STAR
-no-bwa  | --no-bwa                   : disables indexing for BWA
-no-dsj  | --no-diffsplicejunctions   : disables indexing for splice junction analysis
-salm    | --salmon                   : enable and run indexing for salmon
                                        NOTE: unless -it option, converts genome into transcriptome implicitly, thus requires transcript_id tag in gtf


DIFFERENTIAL EXPRESSION ANALYSIS OPTIONS
-c       | --comparisons [path,..]    : triggers differential expression analysis. tabular descriptor file(s) for pairwise comparisons
                                        NOTE: no file implies -no-dsj -no-dea -no-go -no-clust. see below for format information
-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed
                                        NOTE: no fasta file implies -no-map
-it      | --is-transcriptome         : fasta and gtf input is actually a transcriptome converted. uses -no-split, -no-uniq and -no-dsj
                                        NOTE: requires transcript_id as sequence ids and tags in gtf to sum up fractional counts. see genome2transcriptome.pl
-gtf     | --gtf [path]               : annotation gtf input
                                        NOTE: no gtf file implies -no-quant
-go      | --go [path]                : annotation go input
                                        NOTE: no go file implies -no-go
-1       | --fq1 [path,..]            : fastq input. single or first mate. comma separated or a file with all paths
-2       | --fq2 [path,..]            : fastq input. mate pair. comma separated or a file with all paths
-3       | --fq3 [path,..]            : fastq input. UMI sequences. comma separated or a file with all paths
-no-qual | --no-qualityanalysis       : disables intermediate quality analyses and thus adapter inference
                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' quality trimming
-no-clip | --no-clipping              : disables clipping off leading and trailing N's as well as adapter sequences when used with -a
                                      : NOTE: clipping also includes simple 3' quality trimming
-no-pclip| --no-polyntclipping        : disables removal of trailing mono-nucleotide sequences i.e. poly-(A|C|G|T)
-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma separated. default: automatically inferred unless -no-qual option
-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma separated. default: automatically inferred unless -no-qual option
-no-cor  | --no-correction            : disables majority based raw read error correction. recommended for bisulfite sequencing data
-no-rrm  | --no-rrnafilter            : disables rRNA filter
-no-map  | --no-mapping               : disables read alignment and downstream analyses
-d       | --distance                 : maximum read alignment edit distance in %. default: 5
-i       | --insertsize               : maximum allowed insert size for aligning mate pairs. default: 200000 (1000 for -is-transcriptome and -no-split)
                                      : NOTE: does not affect bwa. for segemehl, only multiple alignments are filtered
-no-split| --no-split                 : disables split read mapping. sets default insertsize to 1000
-no-sege | --no-segemehl              : disables mapping by segemehl
-no-star | --no-star                  : disables mapping by STAR
-no-bwa  | --no-bwa                   : disables mapping by BWA given -no-split option
-m       | --mapped [path,..]         : SAM/BAM input. comma separated or a file with all paths (replaces fastq input and processing)
-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
-no-sort | --no-sort                  : disables sorting alignments
-bl      | --blacklist [path|string]  : bedfile of regions or reference/chromosome name to filter alignments for
-sf      | --sizefilter [value:value] : fragment size filtering for alignments to be kept by a given range
                                        NOTE: recommendation for DNA-seq derivatives is 0:1000
-rmd     | --removeduplicates         : enables removing duplicates
-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+)\s*.*
                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
-cmo     | --clipmateoverlaps         : enables clipping of read mate overlaps
-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
-no-quant| --no-quantification        : disables per feature read quantification and TPM calculation plus downstream analyses
-salm    | --salmon                   : enables quantification from quasi-mappings by Salmon
                                        NOTE: unless -it option, converts genome into transcriptome implicitly, thus requires transcript_id tag in gtf
-qf      | --quantifyfeature [string] : switch to other feature with [string]_id tag in gtf for quantification. default: gene, with tag gene_id
-ql      | --quantifylevel [string]   : switch to other feature level for quantification. default: exon
-s       | --strandness [value]       : defines library strandness for all inputs. default: automatically inferred
                                        0 - unstranded
                                        1 - stranded (fr second strand)
                                        2 - reversely stranded (fr first strand)
-no-dsj  | --no-diffsplicejunctions   : disables differential splice junction analysis
-no-dea  | --no-diffexanalysis        : disables differential feature expression analysis plus downstream analyses
-fb      | --featurebiotype [string]  : regex of features with matching [-qf]_(bio)type tag in gtf for clustering and enrichments. default: protein_coding
-no-go   | --no-geneontology          : disables gene ontology enrichment analyses for differentially expressed features and co-expression clusters
-no-clust| --no-clustering            : disables feature co-expression clustering
-cf      | --clusterfilter [value(s)] : decide for a set of features by to be clustered for co-expression. default: 04
                                        NOTE: filters can be turned off by the value null or combined. e.g 01 (equals 10) or 023 or ..
                                        0 - padj <= 0.05 in at least one comparison defined in experiment summary file (see -c)
                                          to take effect, this filter requires upstream performed differential expression analysis
                                        1 - log2foldchange difference >= 0.5 in at least one comparison defined in experiment summary file (see -c)
                                          to take effect, this filter requires upstream performed differential expression analysis
                                        2 - basemean >=5 in at least one comparison defined in experiment summary file (see -c)
                                        3 - discard features within the lower 30% percentile of expression values
                                        4 - TPM >=5 in at least one sample
                                        5 - biotype (see -fb) of matching [-qf]_(bio)type tag in gtf


DIFFERENTIAL METHYLATION ANALYSIS OPTIONS
-b       | --bisulfite [string|value] : triggers methylation analysis. use keyword WGBS or length of RRBS diversity adapters (0 if none)
-c       | --comparisons [path,..]    : triggers differential methylation analysis. tabular descriptor file(s) for pairwise comparisons
                                        NOTE: no file implies -no-dma. see below for format information
-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed
                                        NOTE: no fasta file implies -no-map
-1       | --fq1 [path,..]            : fastq input. single or first mate. comma separated or a file with all paths
-2       | --fq2 [path,..]            : fastq input. mate pair. comma separated or a file with all paths
-3       | --fq3 [path,..]            : fastq input. UMI sequences. comma separated or a file with all paths
-no-qual | --no-qualityanalysis       : disables intermediate quality analyses and thus adapter inference
                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
-no-mspi | --no-mspiselection         : in case of RRBS, disables selection of MspI digested reads (use in case of multi-digestion enzymes)
-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach
                                      : NOTE: in addition, simple 5' clipping is conducted for WGBS data
-no-clip | --no-clipping              : disables clipping off leading and trailing N's as well as adapter sequences when used with -a
                                      : NOTE: clipping also includes simple 3' quality trimming
-no-pclip| --no-polyntclipping        : disables removal of trailing mono-nucleotide sequences i.e. poly-(A|C|G|T)
-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma separated. default: automatically inferred unless -no-qual option
                                      : NOTE: in case of RRBS internally prefixed by NN to address MspI cutting site end-repair bias
-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma separated. default: automatically inferred unless -no-qual option
                                      : NOTE: in case of RRBS, R2 is also 5' clipped by two nucleotides to address MspI cutting site end-repair bias
-no-map  | --no-mapping               : disables read alignment and downstream analyses
-d       | --distance                 : maximum read alignment edit distance in %. default: 5
-i       | --insertsize               : maximum allowed insert for aligning mate pairs. default: 1000
                                      : NOTE: does not affect bwa. for segemehl, only multiple alignments are filtered
-no-sege | --no-segemehl              : disables mapping by segemehl
-no-bwa  | --no-bwa                   : disables mapping by BWA
-m       | --mapped [path,..]         : SAM/BAM input. comma separated or a file with all paths (replaces fastq input and processing)
-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
-no-sort | --no-sort                  : disables sorting alignments
-bl      | --blacklist [path|string]  : bedfile of regions or reference/chromosome name to filter alignments for
-sf      | --sizefilter [value:value] : fragment size filtering for alignments to be kept by a given range
                                        NOTE: recommendation is 0:1000
-rmd     | --removeduplicates         : in case of RRBS, enables removing duplicates
-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+)\s*.*
                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
-no-rmd  | --no-removeduplicates      : in case of WGBS, disables removing duplicates
-no-cmo  | --no-clipmateoverlaps      : disables clipping of read mate overlaps
-no-qual | --no-qualityanalysis       : disables intermediate read and alignment quality analyses and thus adapter inference
                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
-no-mec  | --no-mecall                : disables calling of methylated CpGs plus downstream analyses
-cx      | --context [string]         : regex of context to call methylation in. default: CG
-no-medl | --no-methyldackel          : disables calling of methylation by methyldackel
-no-haarz| --no-haarz                 : disables calling of methylation by haarz
-no-dma  | --no-diffmeanalysis        : disables differential methylation analysis from minimum 10x covered nucleotides
-md      | --min-data [value]         : require at least %/100 or an absolute value of CpG methylation rates per condition (see -c). default: 0.8
-md-cap  | --min-data-cap [value]     : caps/upper bounds required CpG methylation rates per condition (see -md). default: no capping


FUSION DETECTION OPTIONS
-f       | --fusiondetection [string] : triggers gene fusion detection. configure blacklist filter by keyword [null|hg19|GRCh38|GRCm38|GRCm39]
-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed
-1       | --fq1 [path,..]            : fastq input. single or first mate. comma separated or a file with all paths
-2       | --fq2 [path,..]            : fastq input. mate pair. comma separated or a file with all paths
-no-qual | --no-qualityanalysis       : disables intermediate quality analyses and thus adapter inference
                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' quality trimming
-no-clip | --no-clipping              : disables clipping off leading and trailing N's as well as adapter sequences when used with -a
                                      : NOTE: clipping also includes simple 3' quality trimming
-no-pclip| --no-polyntclipping        : disables removal of trailing mono-nucleotide sequences i.e. poly-(A|C|G|T)
-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma separated. default: automatically inferred unless -no-qual option
-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma separated. default: automatically inferred unless -no-qual option
-no-cor  | --no-correction            : disables majority based raw read error correction
-no-rrm  | --no-rrnafilter            : disables rRNA filter
-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
-no-arr  | --no-arriba                : disables fusion detection by Arriba which requieres hg19|hg38|mm10 genome/gtf input
-no-sfus | --no-starfusion            : disables fusion detection by STAR-Fusion which requires CTAT resource as genome/gtf input


PEAK CALLING OPTIONS
-p       | --peakcalling [string]     : triggers peak calling by keyword
                                        CHIP  - account for biomodal asymmetry between sense/antisense mapped reads (transcription factors, PolII,..)
                                        CHIPB - CHIP analysis for broad peaks (histones). affects macs and gopeaks
                                        CUT   - account for more narrow biomodal asymmetry of sparsely mapped reads (CUT&TAG, ChIP-exo)
                                        CUTB  - CUT analysis for broad peaks (histones). affects macs and gopeaks
                                        ATAC  - call up reads instead of fragments
                                        RIP   - call up split-aligned reads from RNA *IP-Seq experiments (CLIP/m6A/meRIP)
-no-idr  | --no-idr                   : disables pseudo-replicates/pool generation to filter loosely called peaks by irreproducible discovery rates
                                        NOTE: -nr* and -tr* options ignored. genrich and peakachu will use all given files as replicates
-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed
                                        NOTE: no fasta file implies -no-map
-gtf     | --gtf [path]               : annotation gtf input
                                        NOTE: required by gem for RNA based *IP-Seq experiments unless given by -s option
-n1      | --normal-fq1 [path,..]     : optional, normal fastq input. single or first mate. comma separated or a file with all paths
-n2      | --normal-fq2 [path,..]     : optional, normal fastq input. mate pair. comma separated or a file with all paths
-n3      | --normal-fq3 [path,..]     : optional, normal fastq input. UMI sequences. comma separated or a file with all paths
-nr1     | --normal-repfq1 [path,..]  : optional, normal replicate fastq input. single or first mate. comma separated or a file with all paths
-nr2     | --normal-repfq2 [path,..]  : optional, normal replicate fastq input. mate pair. comma separated or a file with all paths
-nr3     | --normal-repfq3 [path,..]  : optional, normal replicate fastq input. UMI sequences. comma separated or a file with all paths
-t1      | --treat-fq1 [path,..]      : *IP-Seq fastq input. single or first mate. comma separated or a file with all paths
-t2      | --treat-fq2 [path,..]      : *IP-Seq fastq input. mate pair. comma separated or a file with all paths
-t3      | --treat-fq3 [path,..]      : *IP-Seq fastq input. UMI sequences. comma separated or a file with all paths
-tr1     | --treat-repfq1 [path,..]   : *IP-Seq replicate fastq input. single or first mate. comma separated or a file with all paths
-tr2     | --treat-repfq2 [path,..]   : *IP-Seq replicate fastq input. mate pair. comma separated or a file with all paths
-tr3     | --treat-repfq3 [path,..]   : *IP-Seq replicate fastq input. UMI sequences. comma separated or a file with all paths
-no-qual | --no-qualityanalysis       : disables intermediate quality analyses and thus adapter inference
                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' quality trimming
-no-clip | --no-clipping              : disables clipping off leading and trailing N's as well as adapter sequences when used with -a
                                      : NOTE: clipping also includes simple 3' quality trimming
-no-pclip| --no-polyntclipping        : disables removal of trailing mono-nucleotide sequences i.e. poly-(A|C|G|T)
-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma separated. default: automatically inferred unless -no-qual option
-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma separated. default: automatically inferred unless -no-qual option
-no-cor  | --no-correction            : disables majority based raw read error correction
-no-rrm  | --no-rrnafilter            : disables rRNA filter
-no-map  | --no-mapping               : disables read alignment and downstream analyses
-d       | --distance                 : maximum read alignment edit distance in %. default: 5
-i       | --insertsize               : maximum allowed insert for aligning mate pairs. default: 1000 (200000 for RIP)
                                        NOTE: does not affect bwa. for segemehl, only multiple alignments are filtered
-no-sege | --no-segemehl              : disables mapping by segemehl
-no-star | --no-star                  : disables mapping by STAR
-no-bwa  | --no-bwa                   : disables mapping by BWA
-nm      | --normal-map [path,..]     : normal SAM/BAM input. comma separated or a file with all paths (replaces fastq input and processing)
-nrm     | --normal-repmap [path,..]  : normal replicate SAM/BAM input. comma separated or a file with all paths (replaces fastq input and processing)
-tm      | --treat-map [path,..]      : *IP-Seq SAM/BAM input. comma separated or a file with all paths (replaces fastq input and processing)
-trm     | --treat-repmap [path,..]   : *IP-Seq replicate SAM/BAM input. comma separated or a file with all paths (replaces fastq input and processing)
-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
-no-sort | --no-sort                  : disables sorting alignments
-bl      | --blacklist [path|string]  : bedfile of regions or reference/chromosome name to filter alignments for
-sf      | --sizefilter [value:value] : fragment size filtering for alignments to be kept by a given range
                                        NOTE: relaxed or conservative recommendations
                                        0:120   or 0:100   - ATAC nucleosome free regions
                                        151:280 or 181:250 - ATAC mononucleosomes
                                        131:200            - MNase mononucleosomes
                                        311:480            - ATAC dinucleosomes
                                        501:650 or 551:620 - ATAC trinucleosomes
                                        0:120              - CUT&TAG or CUT&RUN transcription factors, PolII,..
                                        121:1000           - CUT&TAG histones (<120 fragments may be true signal from tagmentation of linkers or surface)
                                        151:1000           - CUT&RUN histones (<150 fragments may be true signal from tagmentation of linkers or surface)
                                        0:1000             - CHIP histones, transcription factors, PolII,..
-rmd     | --removeduplicates         : in case of CUT/CUTB, enables removing duplicates
-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+).*
                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
-no-rmd  | --no-removeduplicates      : disables removing duplicates. not recommended unless reads were mapped on a transcriptome.
-ct      | --cliptn5                  : enables +4/-5 soft-clipping to address Tn5 cutting site end-repair in e.g. ATAC and CUT&TAG experiments
-cmo     | --clipmateoverlaps         : enables clipping of read mate overlaps
-fs      | --fragmentsize [value]     : estimated size of sequenced fragments. default: 150
-s       | --strandness [value]       : defines library strandness for all inputs. default: 0 (automatically inferred for RIP)
                                        0 - unstranded
                                        1 - stranded (fr second strand)
                                        2 - reversely stranded (fr first strand)
-no-call | --no-call                  : disables peak calling and downstream analyses
-no-macs | --no-macs                  : disables peak calling by macs
-no-gem  | --no-gem                   : disables peak calling by gem
-no-homer| --no-homer                 : disables peak calling by homer
-no-rich | --no-genrich               : disables peak calling by genrich
-no-seacr| --no-seacr                 : disables peak calling by seacr
-no-gops | --no-gopeaks               : disables peak calling by gopeaks
-no-peaka| --no-peakachu              : disables peak calling by peakachu
-matk    | --matk                     : enables m6a peak calling by matk/deeprip
-m6a     | --m6aviewer                : enables m6a peak calling by m6aviewer - requieres user interaction
```

## Retrieve SRA datasets
[&#x25B2; back to top](#rippchen)

Use the enclosed script to fetch sequencing data from SRA

```bash
sra-dump.sh -h
```

## Retrieve genomes
[&#x25B2; back to top](#rippchen)

Use the enclosed script to fetch human hg19/hg38 or mouse mm10/mm11 genomes, gene and ontology annotations plus dbSNP and MSigDB. The Plug-n-play CTAT genome resource, made for gene fusion detection and shipped with STAR index, can be selected optionally.

```bash
dlgenome.sh -h
```

The genome, using annotation information, can be converted into a transcriptome or transcript-genome.

```bash
genome2transcriptome.pl
```

## Merge/collate features
[&#x25B2; back to top](#rippchen)

Use the enclosed script to merge e.g. exons of multiple transcripts of a gene with optional offsets for elongation or shrinkage.

```bash
mergexons.sh -h
```

## Generate pseudo replicates
[&#x25B2; back to top](#rippchen)

To shuffle and split NGS raw data in fastq format into two pseudo-replicates, use

```bash
shufnsplitfq.sh -h
```

## Sample info file
[&#x25B2; back to top](#rippchen)

In order to perform desired comparative tasks, some functions require a sample info file.

Assume this input:
<br>

| Treatment   | Replicate 1       | Replicate 2       |
| :---        | :---:             | :---:             |
| wild-type   | path/to/wt1.fq    | path/to/wt2.fq    |
| treatment A | path/to/trA_1.fq  | path/to/trA_2.fq  |
| treatment B | path/to/trB.n1.fq | path/to/trB.n2.fq |

And this desired output (N=2 vs N=2 each):

- wt_vs_A
- wt_vs_B
- A_vs_B

Then the info file should consist of:

- At least 4 tab-separated columns (`<name>`, `<main-factor>`, `NA`, `<replicate>`)
- Optionally, additional factors
- First column needs to consist of unique prefixes of input fastq basenames which can be expand to full file names

|        |     |     |     |        |
| ---    | --- | --- | --- | ---    |
| wt1    | wt  | NA  | N1  | female |
| wt2    | wt  | NA  | N2  | male   |
| trA_1  | A   | NA  | N1  | female |
| trA_2  | A   | NA  | N2  | male   |
| trB.n1 | B   | NA  | N1  | female |
| trB.n2 | B   | NA  | N2  | male   |


## Adapter sequences
[&#x25B2; back to top](#rippchen)

Adapter sequences listed below will be tested by FastQC, extracted from the reports and stored as arrays. In case of paired-end data, unknown adapter sequences will be extracted from mate overlaps utilizing BBMap.

```bash
source <path/of/installation>/latest/bashbone/activate.sh
declare -a fastq_R1=(<path/to/file> [..<path/to/file>])
declare -a fastq_R2=(<path/to/file> [..<path/to/file>])
declare -a adapter_R1 adapter_R2
preprocess::fastqc -t <threads> -o <outdir> -1 fastq_R1 [-2 fastq_R2] -a1 adapter_R1 [-a2 adapter_R2]

```

Further adapter sequences can be found in the Illumina Adapter Sequences Document (<https://www.illumina.com/search.html?q=Illumina Adapter Sequences Document>) or Illumina Adapter Sequences HTML (<https://support-docs.illumina.com/SHARE/adapter-sequences.htm>) and the resource of Trimmomatic (<https://github.com/usadellab/Trimmomatic/tree/main/adapters>), FastQC respectively (<https://github.com/s-andrews/FastQC/blob/master/Configuration>).

The following excerpt is independent of the indexing type, i.e. single, unique dual (UD) or combinatorial dual (CD).

Nextera (Transposase Sequence), TruSight, AmpliSeq, stranded total/mRNA Prep, Ribo-Zero Plus: CTGTCTCTTATACACATCT

TruSeq (Universal) Adapter with A prefix due to 3' primer A-tailing : AGATCGGAAGAGC

TruSeq full length DNA & RNA R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

TruSeq full length DNA MethC R1: AGATCGGAAGAGCACACGTCTGAAC R2: AGATCGGAAGAGCGTCGTGTAGGGA

TruSeq Small RNA 3': TGGAATTCTCGGGTGCCAAGG

TruSeq Small RNA 5': GTTCAGAGTTCTACAGTCCGACGATC

Ovation Methyl-Seq R1: AGATCGGAAGAGC R2: AAATCAAAAAAAC

## Examples
[&#x25B2; back to top](#rippchen)

This section showcases some usages without explaining each parameter in a broader detail. Check out the rippchen help page for more configuration options. Most of them will be opt-out settings.

### Index genomes
[&#x25B2; back to top](#rippchen)

Unknown genomes will be indexed for segemehl, STAR and BWA on the fly. But in order to execute rippchen in multiple instances on a new genome, it should be indexed beforehand. 

```bash
rippchen.sh -x -g <path/to/genome.fa> [-gtf <path/to/genome.gtf>] [-go <path/to/genome.go>]
```

Transcriptomes can be index for Salmon this way

```bash
rippchen.sh -x -g <path/to/transcriptome.fa> [-gtf <path/to/transcriptome.gtf>] --salmon
```

The `-b` option is a switch to index genomes for Methyl-C data analyses.

```bash
rippchen.sh -x -g <path/to/transcriptome.fa> [-gtf <path/to/transcriptome.gtf>] -b WGBS
```

### FastQ data processing
[&#x25B2; back to top](#rippchen)

RNA-seq data pre-processing only with automated adapter identification and clipping, but without Rcorrector for sequencing error correction and without SortMeRNA.

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] --no-correction --no-rrnafilter --no-mapping
```

RNA-seq data pre-processing, mapping by segemehl and alignment post-processing (i.e. unique alignment extraction, sorting, indexing).

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] --no-star --no-quantification
```

RNA-seq data pre-processing, mapping by STAR, but without alignment post-processing (i.e. unique alignment extraction, sorting, indexing).

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] --no-segemehl --no-uniqify --no-sort --no-quantification
```

RNA-seq data pre-processing, mapping by segemehl and alignment post-processing including UMI based de-duplication

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -3 <fastq> --no-star --no-quantification
```

RNA-seq data pre-processing with Illumina universal adapter clipping, mapping by segemehl and STAR, alignment post-processing and quantification

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -a1 AGATCGGAAGAGC [-a2 AGATCGGAAGAGC]
```

RNA-seq data pre-processing and **gene fusion detection**.

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -f <hg19|GRCh38|GRCm38|GRCm39>
```

DNA-seq data pre-processing mapping by BWA and alignment post-processing including removal of PCR duplicates

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] --no-split --no-segemehl --no-star --removeduplicates
```

DNA-seq data pre-processing, segemehl based mapping with adjusted accuracy, enabled duplicate removal and clipping of read mate overlaps (paired-end data only)

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> -2 <fastq> -d 10 --no-split --no-bwa --no-star --removeduplicates --clipmateoverlaps
```

RRBS (without diversity adapter) data pre-processing, mapping by BWA-meth and CpG methylation rate calling utilizing methyldackel (no de-duplication applied)

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -b 0 --no-segemehl --no-haarz
```

WGBS data pre-processing, mapping by segemehl and CpG methylation rate calling utilizing haarz without prior de-duplication

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> [-gtf <gtf>] -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -b WGBS --no-bwa --no-removeduplicates
```

### Alignment data processing
[&#x25B2; back to top](#rippchen)

Alignment files can be used as alternative starting point. Depending on their processing state, alignment post-processing steps have to be disabled (i.e. unique alignment extraction, sorting, de-duplication, mate overlap clipping, indexing).

```bash
find <bam> > <bams>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-m <bams> [--no-uniqify --no-sort ..]
```

### Multiple inputs
[&#x25B2; back to top](#rippchen)

Multiple inputs can be submitted as comma separated list

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>]
```

Or as a file which contains input file paths (one per line)

```bash
find <fastq> > <fastqs>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastqs> [-2 <fastqs>]
```

#### Differential analysis
[&#x25B2; back to top](#rippchen)

Given multiple inputs and a sample info file allows for differential gene expression, co-expression or methylation analysis. The following example conducts RNA-seq data pre-processing, mapping by segemehl and STAR, alignment post-processing including UMI based de-duplication, differential expression, splice junction and co-expression analyses.
<br>

For how to define the sample info file, see the [detailed usage](#usage)

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>] -3 <fastq,fastq,..> -c <sample-info>
```

#### Peak calling
[&#x25B2; back to top](#rippchen)

Infer peaks from paired CHIP-seq datasets of normal/control/input and IP treatment from alignments with mate distance smaller 1000nt (paired-end only) using Macs2, GEM and Peakachu with more stringent cutoffs instead of filtering by irreproducible discovery rate. 

```bash
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq> -n2 <fastq> -t1 <fastq> -t2 <fastq> -p CHIP --no-idr --sizefilter 0:1000 --no-genrich --no-seacr --no-gopeaks
```

Infer peaks from paired CHIP-seq datasets of normal/control and IP treatment with replicates.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq> [-n2 <fastq>] -t1 <fastq> [-t2 <fastq>] -tr1 <fastq> [-tr2 <fastq>]
-p CHIP --sizefilter 0:1000
```

### Start, redo or resume
[&#x25B2; back to top](#rippchen)

List all possible break points and keywords to control rippchen.

```bash
rippchen.sh -dev

md5     : check for md5sums and if necessary trigger genome indexing
fqual   : input quality metrics
mspi    : mspi cutting site selection
trim    : trimming
clip    : adapter clipping (& simple trimming)
pclip   : poly- mono-and di-nucleotide clipping
cor     : raw read correction
rrm     : rRNA filtering
arr     : Arriba gene fusion detection
sfus    : STAR-Fusion detection
sege    : segemehl mapping
star    : STAR mapping
bwa     : BWA mapping
mqual   : mapping/input quality metrics
uniq    : extraction of properly paired and uniquely mapped reads
sort    : sorting and indexing of sam/bam files
blist   : blacklist based alignment filtering
fsel    : fragment size based alignment selection
slice   : better dont touch! slicing of bams for parallelization, needs -prevtmp | --previoustmp [path to rippchen.XXXXXXXXXX]
rmd     : removing duplicates
ctn5    : clipping tn5 end-repaired sites
cmo     : clipping mate overlaps
stats   : fastq preprocessing and mapping statistics
rep     : pooling/generating replicates
macs    : peak calling by macs
gem     : peak calling by gem
homer   : peak calling by homer
rich    : peak calling by genrich
seacr   : peak calling by seacr
gopeaks : peak calling by gopeaks
peaka   : peak calling by peakachu
matk    : peak calling by matk
m6a     : peak calling by m6aViewer
medl    : methylation calling by methyldackel
haarz   : methylation calling by haarz
dma     : differentially methylation analysis
quant   : read quantification
dsj     : differential splice junction analysis
salm    : salmon mapping/quantification
tpm     : TPM calculation
dea     : pca and differential expression analysis
join    : counts joining
clust   : coexpression clustering
go      : go enrichment
```

Use comma separated lists to e.g. skip md5 check and fastq quality analysis.

```bash
rippchen.sh [...] -skip md5,fqual
```

Example how to resume from the segemehl mapping break point after previous data pre-processing.

```bash
rippchen.sh [...] -resume sege
```

Single tasks can be re-computed with the `redo` parameter and a comma separated list of arguments.

```bash
rippchen.sh [...] -redo quant,tpm
```

# Outputs
[&#x25B2; back to top](#rippchen)

Depending on the executed task one ore more of the following directories can be found in the output directory.

## counts
[&#x25B2; back to top](#rippchen)

This directory contains read counts per feature derived from Salmon or featureCounts as well as TPM transformed read counts, DESeq variance stabilization transformed read counts (VSC) and z-score normalized values of all samples.

## coexpressed
[&#x25B2; back to top](#rippchen)

This directory contains WGCNA based co-expression clustering of genes, optionally filtered by biotype, basemean readcount, fold-change and/or TPM cutoff. Clustering is performed on TPM transformed read counts and DESeq variance stabilization transformed read counts (VSC). Sub-directories `vsc` and `tpm` contain clusters plus eigen-gene correlation based or Z-Score based heatmaps as well as expression trajectory footprints. If computed, within each cluster or module directory, functional enrichment tables and plots can be found for GO and MSigDB categories.

- wgcna.cluster.\*.tsv : correlation/p-value/FDR per gene expression level with the trend of a cluster (eigen-gene)
- wgcna.cluster.heatmap.pdf : correlation/FDR per sample with the trend of a cluster (eigen-gene)
- wgcna.cluster.mean.\* : meaned values over replicates
- wgcna.tom.heatmap.pdf : topological overlap matrix from which cluster were inferred, annotated by underlying dendrogram from gene expression correlation values
- cluster.\*.heatmap.pdf : normalized unique alignment counts
- cluster.\*.trajectories.pdf : normalized unique alignment counts as lines per cluster
- cluster.\*.zscores.\* : row normalized values
- cluster.mean.\* : meaned values over replicates

### Cluster enrichment test directories
[&#x25B2; back to top](#rippchen)

Theses directories contain hyper-geometric/over-representation gene set test results using TPM>=1 filtered features in at least one sample as background.

- goenrichment.tsv : enriched GO/MSigDB terms with fdr <= 0.05
- dotplot.\*.pdf,ridgeplot.\*.pdf,barplot.\*.pdf : up to top 30 or all enriched GO/MSigDB terms
- treemap.pdf/semantic_space.pdf : semantically clustered GO terms
- barplot_reduced.\*.pdf : up to top 30 or all enriched GO/MSigDB terms after semantic clustering

## deseq
[&#x25B2; back to top](#rippchen)

This directory contains DESeq2 results of pairwise tests for differentially expressed features between conditions, joined heatmaps across all tests and PCA plots for the first three principle components from all input samples. Sub-directories hold DESeq result tables, volcano plots, MA-plots, tests for impact of additional cofactors/covariates (interaction terms), differences to other conditions (extra effects) and if computed, functional enrichment tables and plots for GO and MSigDB categories.

- experiments.tpm : tpm normalized unique alignment counts
- experiments.vsc : normalized (variance stabilized) unique alignment counts
- experiments.\*.zscores : row normalized values i.e. (value_sample-mean_row)/stdv_sample
- experiments.mean.\* : meaned values over replicates
- deseq.full.annotated.tsv : raw list of differentially expressed genes
- deseq.annotated.tsv : list of differentially expressed genes filtered by fdr <= 0.05
- deseq.\*fcshrunk : shrunked imprecise/over-estimated fold-changes of lowly expressed genes. used to select genes for heatmap plotting.
- heatmap.\*.localclust.pdf : up to top 50 differentially expressed genes by fc-shrunkage. hirarchically clustered columns within each condition
- heatmap.\*.globalclust.pdf : up to top 50 differentially expressed genes by fc-shrunkage. hirarchically clustered columns over all condition
- pca_\*_top\*.pdf : pca (PC1/PC2/PC3) from rlog normalized unique alignment counts of top X most variable genes

### Enrichment test directories
[&#x25B2; back to top](#rippchen)

Theses directories contain gene set enrichment test results and GSEA plots from log2FC shrunked data.

- goenrichment.tsv : enriched GO/MSigDB terms with fdr <= 0.05
- dotplot.\*.pdf,ridgeplot.\*.pdf,barplot.\*.pdf : up to top 30 or all enriched GO/MSigDB terms
- treemap.pdf/semantic_space.pdf : semantically clustered GO terms
- barplot_reduced.\*.pdf : up to top 30 or all enriched GO/MSigDB terms after semantic clustering

## diego
[&#x25B2; back to top](#rippchen)

This directory contains DIEGO results of pairwise test for differential splice junction usage.

## metilene
[&#x25B2; back to top](#rippchen)

This directory contains metilene results of pairwise test for differentially methylated regions.

## stats
[&#x25B2; back to top](#rippchen)

This directory contains simple read count statistics of raw read quality assessing steps, multi- and uniquely mapped reads. Bars are shown in an overlayed fashion. If a mate pair insert size filter was applied, sub-directories contains plot for fragment size distribution per mapper.

## mapped
[&#x25B2; back to top](#rippchen)

This directory contains raw and post-processed alignment files for each mapper applied.

## mecall
[&#x25B2; back to top](#rippchen)

This directory contains full methylation calling reports as well as context and 10-fold coverage filtered methylation rates in bed format plus joined and Z-Score normalized matrices across all input samples.

## peaks
[&#x25B2; back to top](#rippchen)

This directory contains peaks in narrowPeak format, homogenized for all peak callers applied. Unless natively reported by the tool, Macs2 lambda track derived summits and fold-changes are inserted. Sub-directories contain the un-processed results.

## qualities
[&#x25B2; back to top](#rippchen)

This directory contains FastQC quality report and inferred adapter sequences per pre-processing step.

## Pre-processing directories
[&#x25B2; back to top](#rippchen)

Theses directories contain fastq files after different pre-processing steps. Listed in processing order:

- trimmed : Quality trimmed data
- adapterclipped : Adapter clipped data
- polyntclipped : Data after removal of trailing mono- and di-nucleotide sequences i.e. poly-(A|C|G|T) and poly-(AC|AG..|GT). Latter not applied on Methyl-C sequencing data.
- corrected : Data after sequencing error correction
- rrnafiltered : Data after artificial rRNA depletion

# Third-party software
[&#x25B2; back to top](#rippchen)

| Tool | Source | DOI |
| ---  | ---    | --- |
| Arriba        | <https://github.com/suhrig/arriba/>                                 | NA |
| BamUtil       | <https://genome.sph.umich.edu/wiki/BamUtil>                         | 10.1101/gr.176552.114 |
| BBTools       | <https://jgi.doe.gov/data-and-tools/software-tools/bbtools>         | 10.1371/journal.pone.0185056 |
| BWA           | <https://github.com/lh3/bwa>                                        | 10.1093/bioinformatics/btp324 |
| BWA-mem2      | <https://github.com/bwa-mem2/bwa-mem2>                              | 10.1109/IPDPS.2019.00041 |
| BWA-meth      | <https://github.com/brentp/bwa-meth>                                | arXiv:1401.1129 |
| BCFtools      | <http://www.htslib.org/doc/bcftools.html>                           | 10.1093/bioinformatics/btr509 |
| BEDTools      | <https://bedtools.readthedocs.io>                                   | 10.1093/bioinformatics/btq033 |
| gztool        | <https://github.com/circulosmeos/gztool>                            | NA |
| cgpBigWig     | <https://github.com/cancerit/cgpBigWig>                             | NA |
| clusterProfiler | <https://guangchuangyu.github.io/software/clusterProfiler>        | 10.1089/omi.2011.0118 |
| Cutadapt      | <https://cutadapt.readthedocs.io/en/stable>                         | 10.14806/ej.17.1.200 |
| DANPOS3       | <https://github.com/sklasfeld/DANPOS3>                              | 10.1101/gr.142067.112 |
| deepTools2    | <https://deeptools.readthedocs.io/en/latest/index.html>             | 10.1093/nar/gkw257 |
| DESeq2        | <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>   | 10.1186/s13059-014-0550-8 <br> 10.1093/biostatistics/kxw041|
| DEXSeq        | <https://bioconductor.org/packages/release/bioc/html/DEXSeq.html>   | 10.1101/gr.133744.111 |
| DIEGO         | <http://www.bioinf.uni-leipzig.de/Software/DIEGO>                   | 10.1093/bioinformatics/btx690 |
| DGCA          | <https://github.com/andymckenzie/DGCA>                              | 10.1186/s12918-016-0349-1 |
| dupsifter     | <https://github.com/huishenlab/dupsifter>                           | 10.1093/bioinformatics/btad729 |
| fastqc        | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>         | NA |
| featureCounts | <http://subread.sourceforge.net>                                    | 10.1093/bioinformatics/btt656 |
| fgbio         | <http://fulcrumgenomics.github.io/fgbio/>                           | NA |
| freebayes     | <https://github.com/ekg/freebayes>                                  | arXiv:1207.3907 |
| GATK4         | <https://github.com/broadinstitute/gatk>                            | 10.1101/gr.107524.110 <br> 10.1038/ng.806 |
| GEM           | <https://groups.csail.mit.edu/cgs/gem>                              | 10.1371/journal.pcbi.1002638 |
| GNU Parallel  | <https://www.gnu.org/software/parallel/>                            | 10.5281/zenodo.1146014 |
| GoPeaks       | <https://github.com/maxsonBraunLab/gopeaks>                         | 10.1186/s13059-022-02707-w |
| GoSemSim      | <http://bioconductor.org/packages/release/bioc/html/GOSemSim.html>  | 10.1093/bioinformatics/btq064 |
| GSEABase      | <https://bioconductor.org/packages/release/bioc/html/GSEABase.html> | NA |
| HTSeq         | <https://htseq.readthedocs.io>                                      | 10.1093/bioinformatics/btu638 |
| IDR           | <https://github.com/nboley/idr>                                     | 10.1214/11-AOAS466 |
| IGV           | <http://software.broadinstitute.org/software/igv>                   | 10.1038/nbt.1754 |
| Intervene     | <https://github.com/asntech/intervene>                              | 10.1186/s12859-017-1708-7 |
| kent/UCSC utilities | <https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads> | 10.1093/bioinformatics/btq351 |
| khmer         | <https://khmer.readthedocs.io>                                      | 10.12688/f1000research.6924.1 |
| m6aViewer     | <http://dna2.leeds.ac.uk/m6a/>                                      | 10.1261/rna.058206.116 |
| Macs2         | <https://github.com/macs3-project/MACS>                             | 10.1186/gb-2008-9-9-r137 |
| MethylDackel  | <https://github.com/dpryan79/MethylDackel>                          | NA |
| metilene      | <https://www.bioinf.uni-leipzig.de/Software/metilene/>              | 10.1101/gr.196394.115 |
| moose2        | <http://grabherr.github.io/moose2/>                                 | 10.1186/s13040-017-0150-8 |
| PEAKachu      | <https://github.com/tbischler/PEAKachu>                             | NA |
| Picard        | <http://broadinstitute.github.io/picard>                            | NA |
| Platypus      | <https://rahmanteamdevelopment.github.io/Platypus>                  | 10.1038/ng.3036 |
| pugz          | <https://github.com/Piezoid/pugz>                                   | 10.1109/IPDPSW.2019.00042 |
| rapidgzip     | <https://github.com/mxmlnkn/rapidgzip>                              | 10.1145/3588195.3592992 |
| RAxML         | <https://cme.h-its.org/exelixis/web/software/raxml/index.html>      | 10.1093/bioinformatics/btl446 |
| Rcorrector    | <https://github.com/mourisl/Rcorrector>                             | 10.1186/s13742-015-0089-y |
| RSeQC         | <http://rseqc.sourceforge.net>                                      | 10.1093/bioinformatics/bts356 |
| REVIGO        | <https://code.google.com/archive/p/revigo-standalone>               | 10.1371/journal.pone.0021800 |
| RRVGO         | <https://ssayols.github.io/rrvgo>                                   | 10.17912/micropub.biology.000811 |
| Salmon        | <https://combine-lab.github.io/salmon/>                             | 10.1038/nmeth.4197 |
| SalmonTE      | <https://github.com/hyunhwan-jeong/SalmonTE>                        | 10.1142/9789813235533_0016 |
| SAMtools      | <http://www.htslib.org/doc/samtools.html>                           | 10.1093/bioinformatics/btp352 |
| SEACR         | <https://github.com/FredHutch/SEACR>                                | 10.1186/s13072-019-0287-4 |
| segemehl      | <http://www.bioinf.uni-leipzig.de/Software/segemehl>                | 10.1186/gb-2014-15-2-r34 <br> 10.1371/journal.pcbi.1000502 |
| SnpEff        | <https://pcingola.github.io/SnpEff>                                 | 10.4161/fly.19695 |
| SortMeRNA     | <https://bioinfo.lifl.fr/RNA/sortmerna>                             | 10.1093/bioinformatics/bts611 |
| STAR          | <https://github.com/alexdobin/STAR>                                 | 10.1093/bioinformatics/bts635 |
| STAR-Fusion   | <https://github.com/STAR-Fusion/STAR-Fusion/wiki>                   | 10.1101/120295 |
| Trimmomatic   | <http://www.usadellab.org/cms/?page=trimmomatic>                    | 10.1093/bioinformatics/btu170 |
| UMI-tools     | <https://github.com/CGATOxford/UMI-tools>                           | 10.1101/gr.209601.116 |
| VarDict       | <https://github.com/AstraZeneca-NGS/VarDict>                        | 10.1093/nar/gkw227 |
| VarScan       | <http://dkoboldt.github.io/varscan>                                 | 10.1101/gr.129684.111 |
| vcflib        | <https://github.com/vcflib/vcflib>                                  | 10.1371/journal.pcbi.1009123 |
| Vt            | <https://genome.sph.umich.edu/wiki/Vt>                              | 10.1093/bioinformatics/btv112 |
| WGCNA         | <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA> | 10.1186/1471-2105-9-559 |

# Supplementary information
[&#x25B2; back to top](#rippchen)

Rippchen can be executed in parallel instances and thus are able to be submitted as jobs into a queuing system like a Sun Grid Engine (SGE). This could be easily done by utilizing `commander::qsubcmd`. This function makes use of array jobs, which further allows to wait for completion of all jobs, handle single exit codes and alter used resources via `commander::qalter.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
declare -a cmds=()
for i in *R1.fastq.gz; do
	j=${i/R1/R2}
	commander::makecmd -a cmds -v i -v j -c <<-'EOF'
		rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> -1 $i -2 $j
	EOF
  # or simply cmds+=("rippchen.sh [...]")
done
commander::qsubcmd -r -p <env> -t <threads> -i <instances> -n <jobname> -o <logdir> -a cmds
commander::qstat
commander::qalter -p <jobname|jobid> -i <instances>
```

In some rare cases a glibc pthreads bug (<https://sourceware.org/bugzilla/show_bug.cgi?id=23275>) may cause pigz failures (`internal threads error`) and premature termination of tools leveraging on it e.g. Cutadapt and pigz. One can circumvent this by e.g. making use of an alternative pthreads library e.g. compiled without lock elision via `LD_PRELOAD`

```bash
source <path/of/installation/latest/bashbone/activate.sh>
LD_PRELOAD=</path/to/no-elision/libpthread.so.0> <command>
```

# Closing remarks
[&#x25B2; back to top](#rippchen)

Rippchen is a continuously developed pipeline which leverages on bashbone biobash library and which is actively used in my daily work. As a single developer it may take me a while to fix errors and issues. Feature requests cannot be handled so far, but I am happy to receive pull request.
