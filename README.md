# Rippchen

...are tasty! Acquire a taste for peak calling from *IP-Seq experiments or for differential expression- and ontology analysis from RNA-Seq data

Rippchen leverages on bashbone, which is a bash/biobash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses. Rippchen makes use of bashbones best-practice parameterized and run-time tweaked software wrappers and compiles them into a multi-threaded pipeline for analyses of model AND non-model organisms.

## Features

- Most software related parameters will be inferred directly from your data so that all functions require just a minimalistic set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)
- Extensive logging and different verbosity levels
- Start, redo or resume from any point within your data analysis pipeline

## Covered Tasks

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

The whole project is licensed under the GPL v3 (see LICENSE file for details). <br>
**except** the the third-party tools set-upped during installation. Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download

This will download you a copy which includes the latest developments

```bash
git clone --recursive https://github.com/Hoffmann-Lab/rippchen
```

To check out the latest release (irregularly compiled) do

```bash
cd rippchen
git checkout --recurse-submodules $(git describe --tags)
```

# Quick start (without installation)
Please see installation section to get all third-party tools set-upped and subsequently all functions to work properly.

Load the library and list available functions. Each function comes with a usage. Or check out the Rippchen usage.

```bash
source ./activate.sh
bashbone -h
rippchen.sh -h
```

# Installation

## Full installation from scratch

```bash
./setup.sh -i all -d <path/to/installation>
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -h
```

### Upgrade to a new release or if bashbone was previously installed

```bash
./setup.sh -i upgrade -d <path/of/installation>
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -h
```

### Update tools

The setup routine will always install the latest software via conda, which can be updated by running the related setup functions again.

```bash
./setup.sh -i conda_tools -d <path/of/installation>
```

Trimmomatic, segemehl, STAR-Fusion, GEM, mdless and gztool will be installed next to the conda environments. If new releases are available, they will be automatically fetched and installed upon running the related setup functions again.

```bash
./setup.sh -i trimmomatic,segemehl,starfusion,gem,mdless,gztool -d <path/of/installation>
```

# Usage

To load rippchen, bashbone respectively, execute
```bash
source <path/of/installation/latest/rippchen/activate.sh>
bashbone -h
```

In order to get all function work properly, enable bashbone to use conda environments. Conda and bashbone it self can be disabled analogously.
```bash
bashbone -c
# disable conda
bashbone -s
# disable bashbone
bashbone -x
```

Shortcut:

```bash
source <path/of/installation/latest/rippchen/activate.sh> -c true
bashbone -s
```

## Retrieve SRA datasets

Use the enclosed script to fetch sequencing data from SRA

```bash
source <path/of/installation/latest/rippchen/activate.sh> -c true
sra-dump.sh -h
```

## Retrieve genomes

Use the enclosed script to fetch human hg19/hg38 or mouse mm9/mm10/mm11 genomes and annotations. Plug-n-play CTAT genome resource made for gene fusion detection and shipped with STAR index can be selected optionally.

```bash
source <path/of/installation/latest/rippchen/activate.sh> -c true
dlgenome.sh -h
```

## Index genomes

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -x -g <path/to/genome.fa> -gtf <path/to/genome.gtf>
```

## Sample info file

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

- At least 4 columns (`<name>`, `<main-factor>`, `NA`, `<replicate>`)
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

## Examples

This section showcases some usages without explaining each parameter in a broader detail. Check out the rippchen help page for more configuration options. Most of them will be opt-out settings.

### Data pre-processing and mapping plus gene fusion detection

Data pre-processing without mapping by segemehl or STAR and without Rcorrector for sequencing error correction and SortMeRNA for artificial rRNA depletion.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -no-cor -no-rrm -no-sege -no-star
```

Data pre-processing without mapping by segemehl or STAR but enabled fusion detection.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -fusion -no-sege -no-star
```

Data pre-processing, mapping by segemehl and STAR and alignment post-processing (i.e. unique read extraction, sorting, indexing).

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>]
```

Data pre-processing with Illumina universal adapter removal, mapping by segemehl and STAR and alignment post-processing (i.e. unique read extraction, sorting, indexing). Sequences can be found in the Illumina Adapter Sequences Document (<https://www.illumina.com/search.html?q=Illumina Adapter Sequences Document>) or Illumina Adapter Sequences HTML (<https://support-docs.illumina.com/SHARE/adapter-sequences.htm>) and the resource of Trimmomatic (<https://github.com/usadellab/Trimmomatic/tree/main/adapters>), FastQC respectively (<https://github.com/s-andrews/FastQC/blob/master/Configuration>).

The following excerpt is independent of the indexing type, i.e. single, unique dual (UD) or combinatorial dual (CD).

Nextera (Transposase Sequence), TruSight, AmpliSeq, stranded total/mRNA Prep, Ribo-Zero Plus: CTGTCTCTTATACACATCT

TruSeq (Universal) Adapter with A prefix due to 3' primer A-tailing : AGATCGGAAGAGC

TruSeq full length DNA & RNA R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

TruSeq full length DNA MethC R1: AGATCGGAAGAGCACACGTCTGAAC R2: AGATCGGAAGAGCGTCGTGTAGGGA

TruSeq Small RNA: TGGAATTCTCGGGTGCCAAGG

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -a1 AGATCGGAAGAGC [-a2 AGATCGGAAGAGC]
```

Data pre-processing, mapping by segemehl and STAR and disabled post-processing (i.e. unique read extraction, sorting, indexing).

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -no-uniq -no-sort -no-idx
```

Multiple inputs can be submitted as comma separated list.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>]
```

Tweak the amount of allowed mismatches in % during mapping and enable clipping of read mate overlaps (paired-end data only).

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -d 10 -cmo
```

### Differential expression analysis

Pairwise differential gene expression analysis with disabled co-expression analysis.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>] -c <sample-info> -no-clust
```

### Peak calling

Infer peaks from paired datasets of normal/control and ChIP-Seq using macs2, GEM and Peakachu.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq,fastq,...> [-n2 <fastq,fastq,...>] -t1 <fastq,fastq,...> [-t2 <fastq,fastq,...>]
```

Infer peaks from paired datasets of normal/control and *IP-Seq with replicates and without utilizing macs2.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq,fastq,...> [-n2 <fastq,fastq,...>] -t1 <fastq,fastq,...> [-t2 <fastq,fastq,...>] \
-tr1 <fastq,fastq,...> [-tr2 <fastq,fastq,...>] -no-macs
```

### Fusion detection

Predict gene fusions utilizing STAR-Fusion and Arriba.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>] -f
```

### Differential methylation analysis

Pairwise differential methylation expression analysis from WGBS data.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>] -c <sample-info> -b WGBS
```

Pairwise differential methylation expression analysis from RRBS data without diversity adapters, thus disbaled MspI cutting site filtering and just using the in place simple quality trimming along with adapter removal to address end-repair bias.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>] -c <sample-info> -a1 AGATCGGAAGAGC [-a2 AGATCGGAAGAGC] -b 0 -no-mspi -no-trim
```

### Start, redo or resume

List all possible break points and keywords to control Rippchen.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh -dev
```

Use comma separated lists to e.g. skip md5 check and quality analysis.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh [...] -skip md5,fqual
```

Example how to resume from the segemehl mapping break point after previous data pre-processing.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh [...] -resume sege
```

Single tasks can be re-computed with the `redo` parameter and a comma separated list of arguments.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
rippchen.sh [...] -redo quant,tpm
```

# Outputs

## coexpressed

- WGCNA based coexpression clustering of genes, optionally filtered by biotype, which have a basemean readcount >=5 and abs(log2FC)>=0.5 in at least one deseq comparison
- Clustering is performed on TPM transformed read counts and DESeq variance stabilization transformed read counts (VSC)
- Sub-directories `vsc` and `tpm` contain clusters composed of modules plus z-score based heatmaps and log2FC trajectory footprints
- If computed, within each cluster or module directory, GO enrichment tables and plots can be found for the main three GO categories

## counts

- raw read counts per gene or exon as well as TPM transformed read counts, deseq variance stabilization transformed read counts (VSC) and z-score normalized values of all samples

## deseq

- Pairwise tests for differentially expressed genes between groups, contstrained by optionally input factors with heatmaps
- PCA plots for the first three principle components of all sample
- Sub-directories hold volcano plot like ma-plots, VSC based heatmaps of top 50 abs(log2FC) genes and DESeq tables
- If computed, GO enrichment tables and plots can be found for the main three GO categories

## diego

- Pairwise test results for differential splice junctions

## stats

- simple read count statistics of raw read quality assessing steps, multi- and uniquely mapped reads. Bars are shown in an overlayed fashion

# Third-party software

## In production

| Tool | Source | DOI |
| ---  | ---    | --- |
| Arriba        | <https://github.com/suhrig/arriba/>                                 | NA |
| BamUtil       | <https://genome.sph.umich.edu/wiki/BamUtil>                         | 10.1101/gr.176552.114 |
| BEDTools      | <https://bedtools.readthedocs.io>                                   | 10.1093/bioinformatics/btq033 |
| BWA           | <https://github.com/lh3/bwa>                                        | 10.1093/bioinformatics/btp324 |
| Cutadapt      | <https://cutadapt.readthedocs.io/en/stable>                         | 10.14806/ej.17.1.200 |
| DESeq2        | <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>   | 10.1186/s13059-014-0550-8 <br> 10.1093/biostatistics/kxw041|
| DEXSeq        | <https://bioconductor.org/packages/release/bioc/html/DEXSeq.html>   | 10.1101/gr.133744.111 |
| DIEGO         | <http://www.bioinf.uni-leipzig.de/Software/DIEGO>                   | 10.1093/bioinformatics/btx690 |
| DGCA          | <https://github.com/andymckenzie/DGCA>                              | 10.1186/s12918-016-0349-1 |
| fastqc        | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>         | NA |
| featureCounts | <http://subread.sourceforge.net>                                    | 10.1093/bioinformatics/btt656  |
| GEM           | <https://groups.csail.mit.edu/cgs/gem>                              | 10.1371/journal.pcbi.1002638 |
| GNU Parallel  | <https://www.gnu.org/software/parallel/>                            | 10.5281/zenodo.1146014 |
| GoSemSim      | <http://bioconductor.org/packages/release/bioc/html/GOSemSim.html>  | 10.1093/bioinformatics/btq064 |
| GSEABase      | <https://bioconductor.org/packages/release/bioc/html/GSEABase.html> | NA |
| HTSeq         | <https://htseq.readthedocs.io>                                      | 10.1093/bioinformatics/btu638 |
| IDR           | <https://github.com/kundajelab/idr>                                 | 10.1214/11-AOAS466 |
| khmer         | <https://khmer.readthedocs.io>                                      | 10.12688/f1000research.6924.1 |
| Macs2         | <https://github.com/macs3-project/MACS>                             | 10.1186/gb-2008-9-9-r137 |
| metilene      | <https://www.bioinf.uni-leipzig.de/Software/metilene/>              | 10.1101/gr.196394.115 |
| Picard        | <http://broadinstitute.github.io/picard>                            | NA |
| Rcorrector    | <https://github.com/mourisl/Rcorrector>                             | 10.1186/s13742-015-0089-y |
| ReSeqC        | <http://rseqc.sourceforge.net>                                      | 10.1093/bioinformatics/bts356 |
| REVIGO        | <https://code.google.com/archive/p/revigo-standalone>               | 10.1371/journal.pone.0021800 |
| SAMtools      | <http://www.htslib.org/doc/samtools.html>                           | 10.1093/bioinformatics/btp352 |
| segemehl      | <http://www.bioinf.uni-leipzig.de/Software/segemehl>                | 10.1186/gb-2014-15-2-r34 <br> 10.1371/journal.pcbi.1000502 |
| SortMeRNA     | <https://bioinfo.lifl.fr/RNA/sortmerna>                             | 10.1093/bioinformatics/bts611 |
| STAR          | <https://github.com/alexdobin/STAR>                                 | 10.1093/bioinformatics/bts635 |
| STAR-Fusion   | <https://github.com/STAR-Fusion/STAR-Fusion/wiki>                   | 10.1101/120295 |
| Trimmomatic   | <http://www.usadellab.org/cms/?page=trimmomatic>                    | 10.1093/bioinformatics/btu170 |
| WGCNA         | <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA> | 10.1186/1471-2105-9-559 |

## In preparation

| Tool | Source | DOI |
| ---  | ---    | --- |
| clusterProfiler | <https://guangchuangyu.github.io/software/clusterProfiler> | NA |
| HISAT2          | <https://daehwankimlab.github.io/hisat2>                   | 10.1038/nmeth.3317 |

# Supplementary information

Rippchen can be executed in parallel instances and thus are able to be submitted as jobs into a queuing system like a Sun Grid Engine (SGE). This could be easily done by utilizing `commander::qsubcmd`. This function makes use of array jobs, which further allows to wait for completion of all jobs, handle single exit codes and alter used resources via `commander::qalter.

```bash
source <path/of/installation/latest/rippchen/activate.sh>
declare -a cmds=()
for i in *R1.fastq.gz; do
	j=${i/R1/R2}
	sh=job_$(basename $i .R1.fastq.gz)
	commander::makecmd -a cmd1 -c {COMMANDER[0]}<<- CMD
		rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> -1 $i -2 $j
	CMD
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

Bashbone is a continuously developed library and actively used in my daily work. As a single developer it may take me a while to fix errors and issues. Feature requests cannot be handled so far, but I am happy to receive pull request.
