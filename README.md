# Rippchen
---

...are tasty! Acquire a taste for peak calling from *IP-Seq experiments or for differential expression- and ontology analysis from RNA-Seq data

Rippchen leverages on bashbone, which is a bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses. Rippchen makes use of bashbones best-practice parameterized and run-time tweaked software wrappers and compiles them into a multi-threaded pipeline for analyses of model AND non-model organisms. 

## Features

- Most software related parameters will be inferred directly from your data so that all functions require just a minimalistic set of input arguments
- Benefit from a non-root stand-alone installer without need for any prerequisites
- Get genomes, annotations from Ensembl, variants from GATK resource bundle and RAW sequencing data from NCBI Sequence Read Archive (SRA)
- Extensive logging and different verbosity levels
- Start, redo or resume from any point within your data analysis pipeline

## Covered Tasks

- For paired-end and single-end derived raw sequencing or prior mapped read data
- Data preprocessing (quality check, adapter clipping, quality trimming, error correction, artificial rRNA depletion)
- Read alignment and post-processing
  - knapsack problem based slicing of alignment files for parallel task execution
  - sorting, filtering, unique alignment extraction, removal of optical duplicates
  - generation of pools and pseudo-replicates
- Read quantification, TPM and Z-score normalization (automated inference of strand specific library preparation methods)
- Inference of differential expression as well as co-expression clusters
- Detection of differential splice junctions and differential exon usage
- Gene ontology (GO) gene set enrichment and over representation analysis plus semantic similarity based clustering
- Free implementation of Encode3 best-practice ChIP-Seq Peak calling (automated inference of effective genome sizes)
- Peak calling from RIP-Seq, MeRIP-Seq, m6A-Seq and other related *IP-Seq data

# License
---

The whole project is licensed under the GPL v3 (see LICENSE file for details). <br>
**except** the the third-party tools set-upped during installation. Please refer to the corresponding licenses

Copyleft (C) 2020, Konstantin Riege

# Download
---

```bash
git clone --recursive https://github.com/koriege/rippchen.git
git checkout $(git describe --tags)
```

# Quick start (without installation)
Please see installation section to get all third-party tools set-upped and subsequently all functions to work properly.

Load the library and list available functions. Each function comes with a usage. Or check out the Rippchen usage.

```bash
source activate.sh
bashbone -h
rippchen.sh -h
```

# Installation
---

```bash
setup -i all -d <path/to/installation>
source <path/of/installation/activate.sh>
rippchen.sh -h
```

## Update to a newer release

```bash
setup -i upgrade -d <path/of/installation>
source <path/of/installation/activate.sh>
bashbone -h
```

# Usage
---

## Retrieve SRA datasets

Use the enclosed script to fetch sequencing data from SRA

```bash
source <path/of/installation/activate.sh> -c true
sra-dump.sh -h
```

## Retrieve genomes

Use the enclosed script to fetch human hg19/hg38 or mouse mm9/mm10 genomes and annotations.

```bash
source <path/of/installation/activate.sh> -c true
dlgenome.sh -h
```

## Index genomes

```bash
source <path/of/installation/activate.sh>
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
- wt_vs_b
- A_vs_B 

Then the info file should consist of:

- At least 4 columns (`<name>`, `<main-factor>`, `single-end|paired-end`, `<replicate>`)
- Optionally, additional factors
- Unique prefixes of input fastq basenames in the first column which expand to the full file name

|        |     |            |     |        |
| ---    | --- | ---        | --- | ---    |
| wt1    | wt  | single-end | N1  | female |
| wt2    | wt  | single-end | N2  | male   |
| trA_1  | A   | single-end | N1  | female |
| trA_2  | A   | single-end | N2  | male   |
| trB.n1 | B   | single-end | N1  | female |
| trB.n2 | B   | single-end | N2  | male   |


## Examples

This section showcases some usages without explaining each parameter in a broader detail. Check out the Rippchen help page for more configuration options. Most of them will be opt-out settings.

### Data pre-processing and mapping

Data pre-processing without mapping by segemehl or STAR and without Rcorrector for sequencing error correction and SortMeRNA for artificial rRNA depletion.

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -no-cor -no-rrm -no-sege -no-star
```

Data pre-processing, mapping by segemehl and STAR and alignment post-processing (i.e. unique read extraction, sorting, indexing).

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>]
```

Data pre-processing with Illumina universal adapter removal, mapping by segemehl and STAR and alignment post-processing (i.e. unique read extraction, sorting, indexing). More sequences can be found via Illumina Adapter Sequences Document (<https://www.illumina.com/search.html?q=Illumina Adapter Sequences Document&filter=manuals&p=1>), the resource of Trimmomatic (<https://github.com/timflutre/trimmomatic/tree/master/adapters>), FastQC respectively (<https://github.com/s-andrews/FastQC/blob/master/Configuration/contaminant_list.txt>).

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -a1 AGATCGGAAGAG [-a2 AGATCGGAAGAG]
```

Data pre-processing, mapping by segemehl and STAR and disabled post-processing (i.e. unique read extraction, sorting, indexing).

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -no-uniq -no-sort -no-idx
```

Multiple inputs can be submitted as comma separated list.

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>]
```

Tweak the amount of allowed mismatches in % during mapping and enable clipping of read mate overlaps (paired-end data only).

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq> [-2 <fastq>] -d 10 -cmo
```

### Differential expression analysis

Pairwise differential gene expression analysis with disabled co-expression analysis.

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-1 <fastq,fastq,...> [-2 <fastq,fastq,...>] -c <sample-info> -no-clust
```

### Peak calling

Infer peaks from paired datasets of normal/control and ChIP-Seq using macs2 and GEM.

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq,fastq,...> [-n2 <fastq,fastq,...>] -t1 <fastq,fastq,...> [-t2 <fastq,fastq,...>]
```

Infer peaks from paired datasets of normal/control and *IP-Seq with replicates and without utilizing macs2.

```bash
source <path/of/installation/activate.sh>
rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> \
-n1 <fastq,fastq,...> [-n2 <fastq,fastq,...>] -t1 <fastq,fastq,...> [-t2 <fastq,fastq,...>] \
-tr1 <fastq,fastq,...> [-tr2 <fastq,fastq,...>] -no-macs
```

### Start, redo or resume

List all possible break points and keywords to control Rippchen.

```bash
source <path/of/installation/activate.sh>
rippchen.sh -dev
```

Use comma separated lists to e.g. skip md5 check and quality analysis.

```bash
source <path/of/installation/activate.sh>
rippchen.sh [...] -skip md5,qual 
```

Example how to resume from the segemehl mapping break point after previous data pre-processing.

```bash
source <path/of/installation/activate.sh>
rippchen.sh [...] -resume sege
```

Single tasks can be re-computed with the `redo` parameter and a comma separated list of arguments.

```bash
source <path/of/installation/activate.sh>
rippchen.sh [...] -redo quant,tpm
```

# Third-party software
---

## In production

| Tool | Source | DOI |
| ---  | ---    | --- |
| BamUtil       | <https://genome.sph.umich.edu/wiki/BamUtil>                         | 10.1101/gr.176552.114 |
| BEDTools      | <https://bedtools.readthedocs.io>                                   | 10.1093/bioinformatics/btq033 |
| Cutadapt      | <https://cutadapt.readthedocs.io/en/stable>                         | 10.14806/ej.17.1.200 |
| DESeq2        | <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>   | 10.1186/s13059-014-0550-8 |
| DEXSeq        | <https://bioconductor.org/packages/release/bioc/html/DEXSeq.html>   | 10.1101/gr.133744.111 |
| DIEGO         | <http://www.bioinf.uni-leipzig.de/Software/DIEGO>                   | 10.1093/bioinformatics/btx690 |
| DGCA          | <https://github.com/andymckenzie/DGCA>                              | 10.1186/s12918-016-0349-1 |
| fastqc        | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc>         | NA |
| featureCounts | <http://subread.sourceforge.net>                                    | 10.1093/bioinformatics/btt656  |
| GEM           | <https://groups.csail.mit.edu/cgs/gem>                              | 10.1371/journal.pcbi.1002638 |
| GSEABase      | <https://bioconductor.org/packages/release/bioc/html/GSEABase.html> | NA |
| HTSeq         | <https://htseq.readthedocs.io>                                      | 10.1093/bioinformatics/btu638 |
| IDR           | <https://github.com/kundajelab/idr>                                 | 10.1214/11-AOAS466 |
| khmer         | <https://khmer.readthedocs.io>                                      | 10.12688/f1000research.6924.1 |
| Macs2         | <https://github.com/macs3-project/MACS>                             | 10.1186/gb-2008-9-9-r137 |
| Picard        | <http://broadinstitute.github.io/picard>                            | NA |
| Rcorrector    | <https://github.com/mourisl/Rcorrector>                             | 10.1186/s13742-015-0089-y |
| ReSeqC        | <http://rseqc.sourceforge.net>                                      | 10.1093/bioinformatics/bts356 |
| REVIGO        | <https://code.google.com/archive/p/revigo-standalone>               | 10.1371/journal.pone.0021800 |
| segemehl      | <http://www.bioinf.uni-leipzig.de/Software/segemehl>                | 10.1186/gb-2014-15-2-r34 <br> 10.1371/journal.pcbi.1000502 |
| SortMeRNA     | <https://bioinfo.lifl.fr/RNA/sortmerna>                             | 10.1093/bioinformatics/bts611 |
| STAR          | <https://github.com/alexdobin/STAR>                                 | 10.1093/bioinformatics/bts635 |
| Trimmomatic   | <http://www.usadellab.org/cms/?page=trimmomatic>                    | 10.1093/bioinformatics/btu170 |
| WGCNA         | <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA> | 10.1186/1471-2105-9-559 |

## In preparation

| Tool | Source | DOI |
| ---  | ---    | --- |
| BWA             | <https://github.com/lh3/bwa>                               | 10.1093/bioinformatics/btp324 |
| clusterProfiler | <https://guangchuangyu.github.io/software/clusterProfiler> | NA |
| HISAT2          | <https://daehwankimlab.github.io/hisat2>                   | 10.1038/nmeth.3317 |
| STAR-Fusion     | <https://github.com/STAR-Fusion/STAR-Fusion/wiki>          | 10.1101/120295 |

# Supplementary information

Rippchen can be executed in parallel instances and thus is able to be submitted as a job into a queuing system like a Sun Grid Engine (SGE). This could be easily done by amending the following code snipped.

```bash
for i in *R1.fastq.gz; do
	j=${i/R1/R2}
	sh=job_$(basename $i .R1.fastq.gz)
	cat <<- EOF > $sh.sh
		#!/usr/bin/env bash
		source <path/of/installation/activate.sh>
		rippchen.sh -v 2 -t <threads> -g <fasta> -gtf <gtf> -o <outdir> -l <logfile> -tmp <tmpdir> -1 $i -2 $j
	EOF
	echo "rm -f $sh.+(log|err) && qsub -pe <env> <threads> -l 'h=<hostname>|<hostname>' -S /bin/bash -e $sh.err -o $sh.log -V -cwd $sh.sh"
done
```

In some cases a glibc pthreads bug (<https://sourceware.org/bugzilla/show_bug.cgi?id=23275>) may cause pigz failures (`internal threads error`) and premature termination of toola leveraging on it e.g. Cutadapt. One can circumvent this by upgrading the operating system or making use of an alternative pthreads library and `LD_PRELOAD`

```bash
source <path/of/installation/activate.sh>
LD_PRELOAD=/lib64/noelision/libpthread.so.0 rippchen.sh [...]
```