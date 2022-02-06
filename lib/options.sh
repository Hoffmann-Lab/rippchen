#! /usr/bin/env bash
# (c) Konstantin Riege

options::usage() {
	commander::print {COMMANDER[0]}<<- EOF
		DESCRIPTION
		RIPPCHEN are tasty!
		acquire a taste for peak calling from *IP-Seq experiments or for differential expression- and ontology analysis from RNA-Seq data


		VERSION
		$VERSION
		utilizing bashbone $BASHBONE_VERSION


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
		-tmp     | --tmp                      : temporary directory. default: $TMPDIR/rippchen.XXXXXXXXXX
		-r       | --remove                   : remove temporary and unnecessary files upon succesful termination
		-t       | --threads [value]          : number of threads. default: $THREADS
		-xmem    | --max-memory [value]       : total amount of allocatable memory in MB. default: $MAXMEMORY MB i.e. currently available memory
		-mem     | --memory [value]           : allocatable memory per instance of memory greedy tools in MB. defines internal number of parallel instances
		                                        default: $MEMORY which allows for $MTHREADS instances and $MTHREADS SAM/BAM slices according to -xmem
		                                        NOTE: needs to be raised in case of GCThreads, HeapSize or OutOfMemory errors
		-resume  | --resume-from [string]     : resume from a specific pipeline step (see -dev)
		-skip    | --skip [string,..]         : skip specific pipeline step(s). comma seperated (see -dev)
		-redo    | --redo [string,..]         : redo specific pipeline step(s). comma seperated (see -dev)


		INDEXING OPTIONS
		-x       | --index                    : triggers creation of all requiered genome and annotation indices plus md5 sums
		-b       | --bisulfite [string]       : triggers indices for for methylation analyses. use keyword WGBS
		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed (see dlgenome.sh)
		-gtf     | --gtf [path]               : annotation gtf input. default: [-g].gtf (see dlgenome.sh)
		-no-sege | --no-segemehl              : disables indexing for segemehl
		-no-star | --no-star                  : disables indexing for STAR. use this option for indexing of plug-n-play CTAT resource
		                                        NOTE: md5 sum of CTAT [-g].star.idx/SA file needs to be manually added to [-g].md5.sh file afterwards
		-no-bwa  | --no-bwa                   : disables indexing for BWA
		-no-dsj  | --no-diffsplicejunctions   : disables indexing for splice junction analysis



		DIFFERENTIAL EXPRESSION ANALYSIS OPTIONS
		-c       | --comparisons [path,..]    : triggers differential expression analysis. tabular descriptor file(s) for pairwise comparisons
		                                        NOTE: no file implies -no-dsj -no-dea -no-go -no-clust. see below for format information
		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed (see dlgenome.sh)
		                                        NOTE: no fasta file implies -no-map
		-gtf     | --gtf [path]               : annotation gtf input. default: [-g].gtf (see dlgenome.sh)
		                                        NOTE: no gtf file implies -no-quant
		-1       | --fq1 [path,..]            : fastq input. single or first mate. comma seperated or a file with all paths
		-2       | --fq2 [path,..]            : fastq input. mate pair. comma seperated or a file with all paths
		-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' clipping
		-no-clip | --no-clipping              : disables removal of poly N, mono- and di-nucleotide ends as well as adapter sequences when used with -a
		                                      : NOTE: clipping includes simple 3' quality trimming anyways
		-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma seperated
		-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma seperated (can be the same as [-a1]. no revere complement required)
		-no-cor  | --no-correction            : disables majority based raw read error correction. recommended for bisulfite sequencing data
		-no-rrm  | --no-rrnafilter            : disables rRNA filter
		-no-map  | --no-mapping               : disables read alignment and downstream analyses
		-d       | --distance                 : maximum read alignment edit distance in %. default: 5
		-i       | --insertsize               : maximum allowed insert for aligning mate pairs. default: 200000
		-no-split| --no-split                 : disables split read mapping
		-no-sege | --no-segemehl              : disables mapping by segemehl
		-no-star | --no-star                  : disables mapping by STAR
		-no-bwa  | --no-bwa                   : disables mapping by BWA given -no-split option
		-m       | --mapped [path,..]         : SAM/BAM input. comma seperated or a file with all paths (replaces fastq input and processing)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
		-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                  : disables sorting alignments
		-rmd     | --removeduplicates         : enables removing duplicates
		-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+)\s*.*
		                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
		-cmo     | --clipmateoverlaps         : enables clipping of read mate overlaps
		-no-idx  | --no-index                 : disables indexing alignments
		-no-qual | --no-qualityanalysis       : disables intermediate read and alignment quality analyses
		                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
		-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
		-no-quant| --no-quantification        : disables per feature read quantification and TPM calculation plus downstream analyses
		-ql      | --quantifylevel            : switch to other feature type for quantification. default: exon
		                                        NOTE: quantifying using a different feature will break differential expression analysis
		-qt      | --quantifytag              : switch to other feature tag for quantification. default: gene_id
		-s       | --strandness [value]       : defines library strandness for all inputs. default: automatically inferred
		                                        0 - unstranded
		                                        1 - stranded (fr second strand)
		                                        2 - reversely stranded (fr first strand)
		-no-dsj  | --no-diffsplicejunctions   : disables differential splice junction analysis
		-no-dea  | --no-diffexanalysis        : disables differential feature expression analysis plus downstream analyses
		-no-go   | --no-geneontology          : disables gene ontology enrichment analyses for differentially expressed features and co-expression clusters
		-no-clust| --no-clustering            : disables feature co-expression clustering
		-cf      | --clusterfilter [value]    : decide for a set of features by to be clustered for co-expression. default: 0
		                                        0 - padj <= 0.05 in at least one comparison defined in experiment summary file (see -c)
		                                        	to take effect, this filter requires upstream performed differential expression analysis
		                                        1 - log2foldchange difference >= 0.5 in at least one comparison defined in experiment summary file (see -c)
		                                        	to take effect, this filter requires upstream performed differential expression analysis
		                                        2 - basemean/TPM >=5 in at least one comparison/sample
		                                        	to take effect, this filter requires either
		                                        	- upstream performed differential expression analysis as defined in experiment summary file (see -c)
		                                        	- upstream performed quantification and TPM calculation
		                                        3 - discard features within the lower 30% percentile of expression values
		                                        NOTE: filter values can be combined. e.g 01 (equals 10) or 023 or ..
		-cb      | --clusterbiotype [string]  : regex of features with a gene_(bio)type tag (see -gtf) to be clustered. default: .


		DIFFERENTIAL METHYLATION ANALYSIS OPTIONS
		-b       | --bisulfite [string|value] : triggers methylation analysis. use keyword WGBS or length of RRBS diversity adapters (0 if none)
		-c       | --comparisons [path,..]    : triggers differential methylation analysis. tabular descriptor file(s) for pairwise comparisons
		                                        NOTE: no file implies -no-dma. see below for format information
		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed (see dlgenome.sh)
		                                        NOTE: no fasta file implies -no-map
		-1       | --fq1 [path,..]            : fastq input. single or first mate. comma seperated or a file with all paths
		-2       | --fq2 [path,..]            : fastq input. mate pair. comma seperated or a file with all paths
		-no-mspi | --no-mspiselection         : in case of RRBS, disables selection of MspI digested reads (use in case of multi-digestion enzymes)
		-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' clipping
		-no-clip | --no-clipping              : disables removal of poly N, mono- and di-nucleotide ends as well as adapter sequences when used with -a
		                                      : NOTE: clipping includes simple 3' quality trimming anyways
		-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma seperated
		-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma seperated (can be the same as [-a1]. no revere complement required)
		-no-map  | --no-mapping               : disables read alignment and downstream analyses
		-d       | --distance                 : maximum read alignment edit distance in %. default: 5
		-i       | --insertsize               : maximum allowed insert for aligning mate pairs. default: 200000
		-no-sege | --no-segemehl              : disables mapping by segemehl
		-no-bwa  | --no-bwa                   : disables mapping by BWA
		-m       | --mapped [path,..]         : SAM/BAM input. comma seperated or a file with all paths (replaces fastq input and processing)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
		-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                  : disables sorting alignments
		-rmd     | --removeduplicates         : in case of RRBS, enables removing duplicates
		-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+)\s*.*
		                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
		-no-rmd  | --no-removeduplicates      : in case of WGBS, disables removing duplicates
		-no-cmo  | --no-clipmateoverlaps      : disables clipping of read mate overlaps
		-no-idx  | --no-index                 : disables indexing alignments
		-no-qual | --no-qualityanalysis       : disables intermediate read and alignment quality analyses
		                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
		-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
		-no-mec  | --no-mecall                : disables calling of methylated Cs plus downstream analyses
		-no-dma  | --no-diffmeanalysis        : disables differential CpG methylation analysis from minimum 10x covered CpGs
		-md      | --min-data [value]         : require at least %/100 or an absolute value of CpG methylation rates per condition (see -c). default: 0.8
		-md-cap  | --min-data-cap [value]     : caps/upper bounds required CpG methylation rates per condition (see -md). default: no capping


		FUSION DETECTION OPTIONS
		-f       | --fusiondetection [string] : triggers gene fusion detection. use one of the accepted keywords [null|hg19|hg38|mm10]
		-g       | --genome [path]            : genome fasta input. without, only preprocessing is performed (see dlgenome.sh)

		-1       | --fq1 [path,..]            : fastq input. single or first mate. comma seperated or a file with all paths
		-2       | --fq2 [path,..]            : fastq input. mate pair. comma seperated or a file with all paths
		-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' clipping
		-no-clip | --no-clipping              : disables removal of poly N, mono- and di-nucleotide ends as well as adapter sequences when used with -a
		                                      : NOTE: clipping includes simple 3' quality trimming anyways
		-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma seperated
		-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma seperated (can be the same as [-a1]. no revere complement required)
		-no-cor  | --no-correction            : disables majority based raw read error correction. recommended for bisulfite sequencing data
		-no-rrm  | --no-rrnafilter            : disables rRNA filter
		-no-qual | --no-qualityanalysis       : disables intermediate read and alignment quality analyses
		                                        NOTE: given -no-qual and unless -no-stats option, intermediate per file analyses replaced by bulk analysis
		-no-stats| --no-statistics            : disables statistics from read and alignment quality analyses
		-no-arr  | --no-arriba                : disables fusion detection by Arriba which requieres hg19|hg38|mm10 genome/gtf input
		-no-sfus | --no-starfusion            : disables fusion detection by STAR-Fusion which requires CTAT resource as genome/gtf input


		PEAK CALLING OPTIONS
		-n1      | --normal-fq1 [path,..]     : normal fastq input. single or first mate. comma seperated or a file with all paths
		-n2      | --normal-fq2 [path,..]     : normal fastq input. mate pair. comma seperated or a file with all paths
		-nr1     | --normal-repfq1 [path,..]  : normal replicate fastq input. single or first mate. comma seperated or a file with all paths
		-nr2     | --normal-repfq2 [path,..]  : normal replicate fastq input. mate pair. comma seperated or a file with all paths
		-t1      | --treat-fq1 [path,..]      : *IP-Seq fastq input. single or first mate. comma seperated or a file with all paths
		-t2      | --treat-fq2 [path,..]      : *IP-Seq fastq input. mate pair. comma seperated or a file with all paths
		-tr1     | --treat-repfq1 [path,..]   : *IP-Seq replicate fastq input. single or first mate. comma seperated or a file with all paths
		-tr2     | --treat-repfq2 [path,..]   : *IP-Seq replicate fastq input. mate pair. comma seperated or a file with all paths
		-no-trim | --no-trimming              : disables quality trimming utilizing a conservative sliding window approach and simple 5' clipping
		-no-clip | --no-clipping              : disables removal of poly N, mono- and di-nucleotide ends as well as adapter sequences when used with -a
		                                      : NOTE: clipping includes simple 3' quality trimming anyways
		-a1      | --adapter1 [string,..]     : adapter sequence(s) of single or first mate. comma seperated
		-a2      | --adapter2 [string,..]     : adapter sequence(s) of mate pair. comma seperated (can be the same as [-a1]. no revere complement required)
		-no-cor  | --no-correction            : disables majority based raw read error correction. recommended for bisulfite sequencing data
		-no-rrm  | --no-rrnafilter            : disables rRNA filter
		-no-map  | --no-mapping               : disables read alignment and downstream analyses
		-d       | --distance                 : maximum read alignment edit distance in %. default: 5
		-i       | --insertsize               : maximum allowed insert for aligning mate pairs. default: 200000
		-rip     | --rna-ip                   : switch type of *IP-Seq experiment to RNA based *IP-Seq (e.g. meRIP, m6A, CLIP). default: ChIP
		                                        NOTE: implies split read mapping
		-no-map  | --no-mapping               : disables read alignment and downstream analyses
		-no-sege | --no-segemehl              : disables mapping by segemehl
		-no-star | --no-star                  : disables mapping by STAR
		-no-bwa  | --no-bwa                   : disables mapping by BWA given -no-split option
		-nm      | --normal-map [path,..]     : normal SAM/BAM input. comma seperated or a file with all paths (replaces fastq input and processing)
		-nrm     | --normal-repmap [path,..]  : normal replicate SAM/BAM input. comma seperated or a file with all paths (replaces fastq input and processing)
		-tm      | --treat-map [path,..]      : *IP-Seq SAM/BAM input. comma seperated or a file with all paths (replaces fastq input and processing)
		-trm     | --treat-repmap [path,..]   : *IP-Seq replicate SAM/BAM input. comma seperated or a file with all paths (replaces fastq input and processing)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input. default: custom
		-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                  : disables sorting alignments
		-no-rmd  | --no-removeduplicates      : disables removing duplicates - not recommended unless reads were mapped on a transcriptome
		-rx      | --regex [string]           : regex of read name identifier with grouped tile information. default: \S+:(\d+):(\d+):(\d+).*
		                                        NOTE: necessary for successful optical deduplication. to disable or if unavailable, set to null
		-cmo     | --clipmateoverlaps         : enables clipping of read mate overlaps
		-no-idx  | --no-index                 : disables indexing alignments
		-fs      | --fragmentsize [value]     : estimated size of sequenced fragments. default: 200
		-no-macs | --no-macs                  : disables peak calling by macs
		-no-gem  | --no-gem                   : disables peak calling by gem
		-no-peaka| --no-peakachu              : disables peak calling by peakachu
		-m6a     | --m6aviewer                : enables peak calling by m6aviewer - not fully tested yet and requieres user interaction
		-sp      | --strict-peaks             : use a more strict macs and gem parameterization
		-pp      | --pointy-peaks             : enables gem to report only very pointy narrow peaks. use e.g. for CLIP (see -rip) or ChIP-exo
		-no-idr  | --no-idr                   : disregards replicates and irreproducible discovery rates. allows unpaired input (see -t1, -t2)
		                                        NOTE: peakachu will use all files given by -n1/-n2 and -t1/-t2 as replicates


		DIFFERENTIAL ANALYSES DESCRIPTOR FILE FOR OPTION
		-c       | --comparisons [path,..]    : tabular descriptor file(s) for pairwise comparisons

		this file requires 4 or more tab seperated columns without header: sample condition info replicate [factor ..]
		sample names must match unique prefixes of input fastq or SAM/BAM basenames
		pairwise analyses will be performed according to condition column, the primary factor

        assume the folling input: wt1.R1.fq wt1.R2.fq wt2.fq trA_1.fq trA_2.fq trB.n1.fq trB.n2_1.fq trB.n2_2.fq

        example1 to get the following output: wt_vs_A wt_vs_b A_vs_B (N=2 vs N=2 each)
		wt1      wt   NA   N1   PE   female
		wt2      wt   NA   N2   SE   male
		trA_1    A    NA   N1   PE   female
		trA_2    A    NA   N2   PE   male
		trB.n1   B    NA   N1   SE   female
		trB.n2   B    NA   N2   PE   male

        example1 to get the following output: wt_vs_tr (N=2 vs N=4)
		wt1      wt   NA   N1   PE   wt   female
		wt2      wt   NA   N2   SE   wt   male
		trA_1    tr   NA   N1   PE   A    female
		trA_2    tr   NA   N2   PE   A    male
		trB.n1   tr   NA   N3   SE   B    female
		trB.n2   tr   NA   N4   PE   B    male


		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 1
}

options::developer() {
	cat <<- EOF
		DESCRIPTION
		In case of restarting or to resume an analysis use the identifiers below, listed in processing order

		DEVELOPER OPTIONS
		md5   : check for md5sums and if necessary trigger genome indexing
		fqual : quality analysis for input, trim, clip, cor, rrm
		mspi  : mspi cutting site selection
		trim  : trimming
		clip  : adapter clipping (& simple trimming)
		cor   : raw read correction
		rrm   : rRNA filtering

		arr   : Arriba gene fusion detection
		sfus  : STAR-Fusion detection

		sege  : segemehl mapping
		star  : STAR mapping
		bwa   : BWA mapping
		mqual : mapping, uniq, rmd, cmo

		uniq  : extraction of properly paired and uniquely mapped reads
		sort  : sorting and indexing of sam/bam files
		rep   : pooling/generating replicates
		slice : better dont touch! slicing of bams for parallelization, needs -prevtmp | --previoustmp [path]
		rmd   : removing duplicates
		cmo   : clipping mate overlaps
		idx   : intermediate and final bam indexing
		stats : fastq preprocessing and mapping statistics

		macs  : peak calling by macs
		gem   : peak calling by gem
		peaka : peak calling by peakachu
		m6a   : peak calling by m6aViewer

		mec   : methylation calling
		dma   : differentially methylation analysis

		quant : read quantification
		tpm   : TPM calculation
		dsj   : differential splice junction analysis
		dea   : pca and differential expression analysis
		join  : counts joining
		clust : coexpression clustering
		go    : go enrichment
	EOF
	exit 0
}

options::checkopt (){
	local e arg=false
	declare -a mapdata

	case $1 in
		-h        | --help) (options::usage); exit 0;;
		-e        | --env) bashbone -t; exit 0;;
		-dev      | --devel) options::developer;;
		-prevtmp  | --previoustmp) arg=true; PREVIOUSTMPDIR="$2";;
		-resume   | --resume-from) arg=true; options::resume "$2";;
		-skip     | --skip) arg=true; options::skip "$2";;
		-redo     | --redo) arg=true; options::redo "$2";;

		-tmp      | --tmp) arg=true; TMPDIR="$2";;
		-r        | --remove) CLEANUP=true;;
		-v        | --verbosity) arg=true; VERBOSITY=$2;;
		-t        | --threads) arg=true; THREADS=$2;;
		-mem      | --memory) arg=true; MEMORY=$2;;
		-xmem     | --max-memory) arg=true; MAXMEMORY=$2;;

		-x        | --index) INDEX=true;;
		-g        | --genome) arg=true; GENOME="$2";;
		-gtf      | --gtf) arg=true; GTF="$2";;
		-o        | --out) arg=true; OUTDIR="$2";;
		-l        | --log) arg=true; LOG="$2";;

		-1        | --fq1 | -n1 | --normal-fq1) arg=true; nfq1="$2";;
		-2        | --fq2 | -n2 | --normal-fq2) arg=true; nfq2="$2";;
		-nr1      | --normal-repfq1) arg=true; nrfq1="$2";;
		-nr2      | --normal-repfq2) arg=true; nrfq2="$2";;
		-t1       | --treat-fq1) arg=true; tfq1="$2";;
		-t2       | --treat-fq2) arg=true; tfq2="$2";;
		-tr1      | --treat-repfq1) arg=true; rfq1="$2";;
		-tr2      | --treat-repfq2) arg=true; rfq2="$2";;

		-f        | --fusiondetection) arg=true; FUSIONS=$2; [[ $FUSIONS == "null" ]] && noarr=true;;
		-no-arr   | --no-arriba) noarr=true;;
		-no-sfus  | --no-starfusion) nosfus=true;;

		-a1       | --adapter1) arg=true; mapfile -t -d ',' ADAPTER1 < <(printf '%s' "$2");;
		-a2       | --adapter2) arg=true; mapfile -t -d ',' ADAPTER2 < <(printf '%s' "$2");;
		-d        | --distance) arg=true; DISTANCE=$2;;
		-i        | --insertsize) arg=true; INSERTSIZE=$2;;
		-no-split | --no-split) nosplitreads=true;;
		-no-qual  | --no-qualityanalysis) noqual=true;;
		-no-trim  | --no-trimming) notrim=true;;
		-no-clip  | --no-clipping) noclip=true;;
		-no-cor   | --no-correction) nocor=true;;
		-no-rrm   | --no-rrnafilter) norrm=true;;
		-no-map   | --no-mapping) nosege=true; nostar=true; nobwa=true;;
		-no-sege  | --no-segemehl) nosege=true;;
		-no-star  | --no-star) nostar=true;;
		-no-bwa   | --no-bwa) nobwa=true;;

		-m        | --mapped | -nm | --normal-mapped) arg=true; nmap="$2";;
		-nrm      | --normal-repmapped) arg=true; nrmap="$2";;
		-tm       | --treat-mapped) arg=true; tmap="$2";;
		-trm      | --treat-repmapped) arg=true; rmap="$2";;
		-mn       | --mapper-name) arg=true; MAPNAME="$2";;

		-no-uniq  | --no-uniqify) nouniq=true;;
		-no-sort  | --no-sort) nosort=true;;
		-rx       | --regex) arg=true; REGEX="$2";;
		-rmd      | --removeduplicates) normd=false;;
		-no-rmd   | --no-removeduplicates) normd=true;;
		-cmo      | --clipmateoverlaps) nocmo=false;;
		-no-cmo   | --no-clipmateoverlaps) nocmo=true;;
		-no-idx   | --no-index) noidx=true;;
		-no-stats | --no-statistics) nostats=true;;

		-rip      | --rna-ip) RIPSEQ=true;;
		-fs       | --fragmentsize) arg=true; FRAGMENTSIZE=$2;;
		-sp       | --strict-peaks) STRICTPEAKS=true;;
		-pp       | --pointy-peaks) POINTYPEAKS=true;;
		-no-macs  | --no-macs) nomacs=true;;
		-no-gem   | --no-gem) nogem=true;;
		-no-peaka | --no-peakachu) nopeaka=true;;
		-no-idr   | --no-idr) noidr=true;;
		-m6a      | --m6aviewer) nom6a=false;;

		-c        | --comparisons) arg=true; mapfile -t -d ',' COMPARISONS < <(printf '%s' "$2");;

		-b        | --bisulfite) arg=true; DIVERSITY="$2"; BISULFITE=true; RRBS=false; nocor=true; norrm=true; [[ "$DIVERSITY" == "WGBS" ]] && { RRBS=false; normd=${normd:-false}; } || { RRBS=true; normd=${normd:-true}; };;
		-no-mspi  | --no-mspiselection) nomspi=true;;
		-no-mec   | --no-mecall) nomec=true;;
		-no-dma   | --no-diffmeanalysis) nodma=true;;
		-md       | --min-data) arg=true; MINDATA=$2;;
		-md-cap   | --min-data-cap) arg=true; MINDATACAP=$2;;

		-s        | --strandness) arg=true; STRANDNESS=$2;;
		-ql       | --quantifylevel) arg=true; QUANTIFYFLEVEL="$2";;
		-qt       | --quantifytag) arg=true; QUANTIFYTAG="$2";;
		-cf       | --clusterfilter) arg=true; CLUSTERFILTER=$2;;
		-cb       | --clusterbiotype) arg=true; CLUSTERBIOTYPE="$2";;
		-no-quant | --no-quantification) noquant=true; nodsj=true; nodea=true; noclust=true; nogo=true;;
		-no-dsj   | --no-diffsplicejunctions) nodsj=true;;
		-no-dea   | --no-diffexanalysis) nodea=true;;
		-no-clust | --no-clustering) noclust=true;;
		-no-go    | --no-geneontology) nogo=true;;

		-*) commander::printerr "illegal option $1"; return 1;;
		*) commander::printerr "illegal option $2"; return 1;;
	esac

	$arg && {
		[[ ! $2 ]] && commander::printerr "argument missing for option $1" && return 1
		[[ "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
		return 0
	} || {
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
		return 0
	}
}

options::resume(){
	local s enable=false
	# don't Smd5, Sslice !
	for s in fqual mspi trim clip cor rrm arr sfus sege star bwa mqual uniq sort rep rmd cmo idx stats macs gem peaka m6a mec dma quant tpm dsj dea join clust go; do
		eval "\${S$s:=true}" # unless S$s already set to false by -redo, do skip
		$enable || [[ "$1" == "$s" ]] && {
			enable=true
			eval "S$s=false"
		}
	done
}

options::skip(){
	local x s
	declare -a mapdata
	mapfile -t -d ',' mapdata < <(printf '%s' "$1")
	for x in "${mapdata[@]}"; do
		for s in md5 fqual mspi trim clip cor rrm arr sfus sege star bwa mqual uniq sort rep slice rmd cmo idx stats macs gem peaka m6a mec dma quant tpm dsj dea join clust go; do
			[[ "$x" == "$s" ]] && eval "S$s=true"
		done
	done
}

options::redo(){
	local x s
	declare -a mapdata
	mapfile -t -d ',' mapdata < <(printf '%s' "$1")
	for s in fqual mspi trim clip cor rrm arr sfus sege star bwa mqual uniq sort rep rmd cmo idx stats macs gem peaka m6a mec dma quant tpm dsj dea join clust go; do
		eval "\${S$s:=true}" # unless (no|S)$s alredy set to false by -resume, do skip
	done
	for x in "${mapdata[@]}"; do
		for s in fqual mspi trim clip cor rrm arr sfus sege star bwa mqual uniq sort rep rmd cmo idx stats macs gem peaka m6a mec dma quant tpm dsj dea join clust go; do
			[[ "$x" == "$s" ]] && eval "S$s=false"
		done
	done
}
