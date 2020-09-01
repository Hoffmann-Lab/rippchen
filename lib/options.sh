#! /usr/bin/env bash
# (c) Konstantin Riege

options::usage() {
	commander::print {COMMANDER[0]}<<- EOF
		DESCRIPTION
		RIPPCHEN are tasty!
		acquire a taste for peak calling from *IP-Seq experiments or for differential expression- and ontology analysis from RNA-Seq data

		VERSION
		$VERSION
		utilizing bashbone $BASHBONEVERSION

		SYNOPSIS INDEXING
		rippchen.sh -x -g genome.fa -gtf genome.gtf

		SYNOPSIS PREPROCESSING
		rippchen.sh -1 ctr1.fq,ctr2.fq,treat1.fq,treat2.fq -g genome.fa -no-quant

		SYNOPSIS DIFFERENTIAL EXPRESSION ANALYSIS
		rippchen.sh -1 ctr1.fq,ctr2.fq,treat1.fq,treat2.fq -g genome.fa -gtf genome.gtf -c cmp.txt

		SYNOPSIS PEAK CALLING
		rippchen.sh -1 ctrA.fq,ctrB.fq -t1 treatA1.fq,treatB1.fq -r1 treatA2.fq,treatB2.fq -g genome.fa

		SYNOPSIS PEAK CALLING AND DIFFERENTIAL EXPRESSION ANALYSIS
		rippchen.sh -1 ctrA.fq,ctrB.fq -t1 treatA1.fq,treatB1.fq -r1 treatA2.fq,treatB2.fq -g genome.fa -gt genome.gtf -c cmp.txt

		BASIC OPTIONS
		-h       | --help                     : prints this message
		-dev     | --devel                    : prints extended pipeline options
		-r       | --remove                   : remove temporary and unnecessary files upon succesful termination
		-v       | --verbosity [value]        : set level of verbosity - default: 0
		                                        0 - get simple status updates
		                                        1 - get status updates and commands
		                                        2 - get full output
		-x       | --index                    : create all requiered genome indices if necessary and exit
		-g       | --genome [path]            : genome fasta input, without only preprocessing is performed
		-gtf     | --gtf [path]               : annotation gtf input - optional, default: [-g].gtf
		-a1      | --adapter1 [string,..]     : adapter sequence(s) - optional. single or first pair, comma seperated
		-a2      | --adapter2 [string,..]     : adapter sequence(s) - optional. second pair, comma seperated
		-o       | --out [path]               : output directory - default: $OUTDIR
		-l       | --log [path]               : output directory - default: $OUTDIR/run.log
		-tmp     | --tmp                      : temporary directory - default: $TMPDIR/rippchen.XXXXXXXXXX
		-t       | --threads [value]          : threads - predicted default: $THREADS
		-mem     | --memory [value]           : amout of memory for creating bam slices and processing them in parallel instances
		                                        available: $MAXMEMORY
		                                        default: 30000 (allows for $MTHREADS instances)
		                                        NOTE: needs to be raised in case of GCThreads, HeapSize or OutOfMemory errors
		-resume  | --resume-from [string]     : resume from a specific pipeline step - see -dev|--devel
		-skip    | --skip [string,..]         : skip specific pipeline step(s) - see -dev|--devel, comma seperated
		-redo    | --redo [string,..]         : just rerun specific pipeline step(s) - see -dev|--devel, comma seperated
		-no-qual | --no-qualityanalysis       : disables read quality analysis
		-no-trim | --no-trimming              : disables quality trimming
		-no-clip | --no-clipping              : disables removal of adapter sequences if -a|--adapter is used
		-no-cor  | --no-correction            : disables majority based raw read error correction
		-no-rrm  | --no-rrnafilter            : disables rRNA filter
		-no-sege | --no-segemehl              : disables mapping by Segemehl
		-no-star | --no-star                  : disables mapping by STAR
		-no-stats| --no-statistics            : disables fastq preprocessing and mapping statistics

		ALIGNMENT OPTIONS
		-d       | --distance                 : maximum read alignment edit distance in % - default: 5
		-i       | --insertsize               : maximum allowed insert for aligning mate pairs - default: 200000
		-no-split| --no-split                 : disable split read mapping
		-no-uniq | --no-uniqify               : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                  : disables sorting alignments
		-no-idx  | --no-index                 : disables indexing alignments
		-cmo     | --clipmateoverlaps         : enable clipping of read mate overlaps

		QUANTIFICATION OPTIONS
		-ql      | --quantifylevel            : switch to other feature type for quantification - default: exon
		                                        NOTE: quantifying using a different feature will break differential expression analysis
		-qt      | --quantifytag              : switch to other feature tag for quantification - default: gene_id
		-no-quant| --no-quantification        : disables per feature read quantification plus downstream analyses

		PEAK CALLING OPTIONS
		-rip     | --rna-ip                   : switch type of *IP-Seq experiment to RNA based *IP-Seq (e.g. meRIP, m6A, CLIP)
		-n1      | --normal-fq1 [path,..]     : normal fastq input - single or first pair, comma seperated or a file with all paths
		-n2      | --normal-fq2 [path,..]     : normal fastq input - optional. second pair, comma seperated or a file with all paths
		-nr1     | --normal-repfq1 [path,..]  : normal replicate fastq input - optional. single or first pair, comma seperated or a file with all paths
		-nr2     | --normal-repfq2 [path,..]  : normal replicate fastq input - optional. second pair, comma seperated or a file with all paths
		-t1      | --treat-fq1 [path,..]      : *IP-Seq fastq input - single or first pair, comma seperated or a file with all paths
		-t2      | --treat-fq2 [path,..]      : *IP-Seq fastq input - optional. second pair, comma seperated or a file with all paths
		-tr1     | --treat-repfq1 [path,..]   : *IP-Seq replicate fastq input - optional. single or first pair, comma seperated or a file with all paths
		-tr2     | --treat-repfq2 [path,..]   : *IP-Seq replicate fastq input - optional. second pair, comma seperated or a file with all paths
		-nm      | --normal-map [path,..]     : normal SAM/BAM input - comma seperated or a file with all paths (replaces fastq input)
		-nrm     | --normal-repmap [path,..]  : normal replicate SAM/BAM input - optional. comma seperated or a file with all paths (replaces fastq input)
		-tm      | --treat-map [path,..]      : *IP-Seq SAM/BAM input - comma seperated or a file with all paths (replaces fastq input)
		-trm     | --treat-repmap [path,..]   : *IP-Seq replicate SAM/BAM input - optional. comma seperated or a file with all paths (replaces fastq input)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input - optional. default: custom

		-f       | --fragmentsize [value]     : fragment size of sequenced mate pairs - default: 150
		-rx      | --regex [string]           : regex of read name identifier with grouped tile information - default: ^\S+:(\d+):(\d+):(\d+)\s*.*
		                                        NOTE: necessary for sucessful deduplication, if unavailable set to 'null'
		-no-rmd  | --no-removeduplicates      : disables removing duplicates - not recommended
		-no-macs | --no-macs                  : disables peak calling by macs
		-no-gem  | --no-gem                   : disables peak calling by gem

		DIFFERENTIAL EXPRESSION ANALYSIS OPTIONS
		-1       | --fq1 [path,..]            : fastq input - single or first pair, comma seperated or a file with all paths
		-2       | --fq2 [path,..]            : fastq input - optional. second pair, comma seperated or a file with all paths
		-m       | --mapped [path,..]         : SAM/BAM input - comma seperated or a file with all paths (replaces fastq input)
		-mn      | --mapper-name [string]     : name to use for output subdirectories in case of SAM/BAM input - optional. default: custom
		-rmd     | --removeduplicates         : enable removing duplicates - not recommended
		-rx      | --regex [string]           : regex of read name identifier with grouped tile information - default: ^\S+:(\d+):(\d+):(\d+)\s*.*
		                                        NOTE: necessary for sucessful deduplication, if unavailable set to 'null'
		-c       | --comparisons [path,..]    : experiment info file(s) for pairwise analyses according to column condition (primary factor)
		                                        format: 4 or more columns, seperated by tab or space(s)
		                                          sample   condition   [single-end|paired-end]   replicate   [factor1   factor2   ..]
		                                        NOTE: samples must match unique prefixes of input fastq or SAM/BAM basenames
		                                        example1: for input wt1.R1.fq wt1.R2.fq wt2.fq trA_1.fq trA_2.fq trB.n1.fq trB.n2_1.fq trB.n2_2.fq
		                                          wt1      wt   paired-end   N1   PE   female
		                                          wt2      wt   single-end   N2   SE   male
		                                          trA_1    A    single-end   N1   PE   female
		                                          trA_2    A    single-end   N2   PE   male
		                                          trB.n1   B    single-end   N1   SE   female
		                                          trB.n2   B    paired-end   N2   PE   male
		                                        output: wt_vs_A wt_vs_b A_vs_B (N=2 vs N=2 each)
		                                        example2:
		                                          wt1      wt   paired-end   N1   PE   wt   female
		                                          wt2      wt   single-end   N2   SE   wt   male
		                                          trA_1    tr   single-end   N1   PE   A    female
		                                          trA_2    tr   single-end   N2   PE   A    male
		                                          trB.n1   tr   single-end   N3   SE   B    female
		                                          trB.n2   tr   paired-end   N4   PE   B    male
		                                        output: wt_vs_tr (N=2 vs N=4)
		-no-dsj  | --no-diffsplicejunctions   : disables differential splice junction analysis
		-no-dea  | --no-diffexanalysis        : disables differential feature expression analysis plus downstream analyses
		-no-go   | --no-geneontology          : disables gene ontology enrichment analyses for differentially expressed features and co-expression clusters
		-no-clust| --no-clustering            : disables feature co-expression clustering
		-cf      | --clusterfilter [value]    : decide for a set of features by to be clustered for co-expression - default: 0
		                                        0 - padj <= 0.05 in at least one comparison defined in experiment info file (see -c)
		                                        	to take effect, this filter requires upstream performed differential expression analysis
		                                        1 - log2foldchange difference >= 0.5 in at least one comparison defined in experiment info file (see -c)
		                                        	to take effect, this filter requires upstream performed differential expression analysis
		                                        2 - basemean/TPM >=5 in at least one comparison/sample
		                                        	to take effect, this filter requires either
		                                        	- upstream performed differential expression analysis as defined in experiment info file (see -c)
		                                        	- upstream performed quantification and TPM calculation
		                                        3 - discard features within the lower 30% percentile of expression values
		                                        NOTE: filter values can be combined. e.g 01 (equals 10) or 023 or ..
		-cb      | --clusterbiotype [string]  : decide for features annotated by a gene_biotype tag (see -gtf) to be clustered for co-expression

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 0
}

options::developer() {
	cat <<- EOF
		DESCRIPTION
		In case of restarting or to resume an analysis use the identifiers below, listed in processing order

		DEVELOPER OPTIONS
		md5   : check for md5sums and if necessary trigger genome indexing
		qual  : quality analysis for input and trim, clip, cor, rrm
		trim  : trimming
		clip  : adapter clipping
		cor   : raw read correction
		rrm   : rRNA filtering
		sege  : Segemehl mapping
		star  : STAR mapping
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
	local arg=false
	case $1 in
		-h   | --help) options::usage;;
		-dev | --devel) options::developer;;

		-r   | --remove) CLEANUP=true;;
		-v   | --verbosity) arg=true; VERBOSITY=$2;;
		-t   | --threads) arg=true; THREADS=$2;;
		-mem | --memory) arg=true; MEMORY=$2;;
		-g   | --genome) arg=true; GENOME=$2;;
		-gtf | --gtf) arg=true; GTF=$2;;
		-o   | --out) arg=true; OUTDIR=$2;;
		-l   | --log) arg=true; LOG=$2;;
		-tmp | --tmp) arg=true; TMPDIR=$2;;
		-prevtmp | --previoustmp) arg=true; PREVIOUSTMPDIR=$2;;

		-1   | --fq1 | -n1 | --normal-fq1) arg=true; nfq1=$2;;
		-2   | --fq2 | -n2 | --normal-fq2) arg=true; nfq2=$2;;
		-nr1 | --normal-repfq1) arg=true; nrfq1=$2;;
		-nr2 | --normal-repfq2) arg=true; nrfq2=$2;;
		-t1  | --treat-fq1) arg=true; tfq1=$2;;
		-t2  | --treat-fq2) arg=true; tfq2=$2;;
		-tr1 | --treat-repfq1) arg=true; rfq1=$2;;
		-tr2 | --treat-repfq2) arg=true; rfq2=$2;;

		-m   | --mapped | -nm | --normal-mapped) arg=true; nmap=$2;;
		-nrm | --normal-repmapped) arg=true; nrmap=$2;;
		-tm  | --treat-mapped) arg=true; tmap=$2;;
		-trm | --treat-repmapped) arg=true; rmap=$2;;
		-mn  | --mapper-name) arg=true; MAPNAME=$2;;

		-rx  | --regex) arg=true; REGEX=$2;;
		-rip | --rna-ip) RIPSEQ=true;;
		-c   | --comparisons) arg=true; mapfile -t -d ',' COMPARISONS < <(printf '%s' "$2");;
		-a1  | --adapter1) arg=true; mapfile -t -d ',' ADAPTER1 < <(printf '%s' "$2");;
		-a2  | --adapter2) arg=true; mapfile -t -d ',' ADAPTER2 < <(printf '%s' "$2");;
		-d   | --distance) arg=true; DISTANCE=$2;;
		-f   | --fragmentsize) arg=true; FRAGMENTSIZE=$2;;
		-i   | --insertsize) arg=true; INSERTSIZE=$2;;
		-ql  | --quantifylevel) arg=true; QUANTIFYFLEVEL=$2;;
		-qt  | --quantifytag) arg=true; QUANTIFYTAG=$2;;
		-cf  | --clusterfilter) arg=true; CLUSTERFILTER=$2;;
		-cb  | --clusterbiotype) arg=true; CLUSTERBIOTYPE=$2;;

		-resume | --resume-from) arg=true; options::resume "$2";;
		-skip | --skip) arg=true; options::skip "$2";;
		-redo | --redo) arg=true; options::redo "$2";;

		-x  | --index) INDEX=true;;
		-no-qual  | --no-qualityanalysis) noqual=true;;
		-no-clip  | --no-clipping) noclip=true;;
		-no-trim  | --no-trimming) notrim=true;;
		-no-cor   | --no-correction) nocor=true;;
		-no-rrm   | --no-rrnafilter) norrm=true;;
		-no-split | --no-split) nosplitreads=true;;
		-no-sege  | --no-segemehl) nosege=true;;
		-no-star  | --no-star) nostar=true;;
		-no-uniq  | --no-uniqify) nouniq=true;;
		-no-sort  | --no-sort) nosort=true;;
		-no-idx   | --no-index) noidx=true;;
		-no-stats | --no-statistics) nostats=true;;
		-no-rmd   | --no-removeduplicates) normd=true;;
		-rmd      | --removeduplicates) normd=false;;
		-cmo      | --clipmateoverlaps) nocmo=false;;
		-no-macs  | --no-macs) nomacs=true;;
		-no-gem   | --no-gem) nogem=true;;
		-no-quant | --no-quantification) noquant=true; nodea=true; noclust=true; nogo=true;;
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
	for s in qual trim clip cor rrm sege star uniq sort rep rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
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
		for s in md5 qual trim clip cor rrm sege star uniq sort rep slice rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
			[[ "$x" == "$s" ]] && eval "S$s=true"
		done
	done
}

options::redo(){
	local x s
	declare -a mapdata
	mapfile -t -d ',' mapdata < <(printf '%s' "$1")
	for s in qual trim clip cor rrm sege star uniq sort rep rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
		eval "\${S$s:=true}" # unless (no|S)$s alredy set to false by -resume, do skip
	done
	for x in "${mapdata[@]}"; do
		for s in qual trim clip cor rrm sege star uniq sort rep rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
			[[ "$x" == "$s" ]] && eval "S$s=false"
		done
	done
}
