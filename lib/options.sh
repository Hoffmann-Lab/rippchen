#! /usr/bin/env bash
# (c) Konstantin Riege

		# -dexnew  | --dexnewprepare [value]  : switch DEXSeq annotation file preparation method - default: 0
		#                                       0 - use preparation script shipped with DEXSeq
		#                                       1 - do not merge transcripts and if present in GTF info and use protein_coding biotype only
		#                                       2 - see [1] and additionally if feature name matches UTR/utr, exclude them
		# -no-iso  | --no-isoformanalysis     : disables differential isoform analysis by DEXSeq and annotation preparation
		# -cf      | --clusterfilter [value]  : decide for a set of differntially expressed features to be clustered for co-expression - default: 0
		#                                       0 - padj <= 0.05 in at least one comparison defined via -c
		#                                       1 - padj <= 0.05 and a log2foldchange difference >= 0.5 in at least one comparison defined via -c
		#                                       2 - if present in GTF info, use features of protein_coding biotype only
		#                                       3 - discard features within the 30% percentile of lowest expression values
		#                                       20 - 2 + 0
		#                                       21 - 2 + 1
		#                                       23 - 2 + 3
		#                                       230 - 2 + 3 + 0
		#                                       231 - 2 + 3 + 1

options::usage() {
	cat <<- EOF
		DESCRIPTION
		RIPPCHEN are tasty!
		acquire a taste for peak calling from *IP-Seq experiments or for differential expression- and ontology analysis from RNA-Seq data

		VERSION
		$version

		SYNOPSIS PREPROCESSING
		$(basename $0) -1 ctr1.fq,ctr2.fq,treat1.fq,treat2.fq -g genome.fa -no-quant

		SYNOPSIS DIFFERENTIAL EXPRESSION ANALYSIS
		$(basename $0) -1 ctr1.fq,ctr2.fq,treat1.fq,treat2.fq -g genome.fa -c cmp.txt

		SYNOPSIS PEAK CALLING
		$(basename $0) -1 ctrA.fq,ctrB.fq -t1 treatA1.fq,treatB1.fq -r1 treatA2.fq,treatB2.fq -g genome.fa

		SYNOPSIS PEAK CALLING AND DIFFERENTIAL EXPRESSION ANALYSIS
		$(basename $0) -1 ctrA.fq,ctrB.fq -t1 treatA1.fq,treatB1.fq -r1 treatA2.fq,treatB2.fq -g genome.fa -c cmp.txt

		BASIC OPTIONS
		-h       | --help                   : prints this message
		-dev     | --devel                  : prints extended pipeline options
		-r       | --remove                 : clean up after successful termination
		-v       | --verbosity [value]      : set level of verbosity - default: 0
		                                      0 - get simple status updates
		                                      1 - get status updates and commands
		                                      2 - get full output
		-g       | --genome [path]          : genome fasta input, without only preprocessing is performed
		-gtf     | --gtf [path]             : annotation gtf input - optional, default: genome.fasta.gtf
		-a1      | --adapter1 [string,..]   : adapter sequence(s) - optional. single or first pair, comma seperated
		-a2      | --adapter2 [string,..]   : adapter sequence(s) - optional. second pair, comma seperated
		-o       | --out [path]             : output directory - default: $OUTDIR
		-l       | --log [path]             : output directory - default: $OUTDIR/run.log
		-tmp     | --tmp                    : temporary directory - default: $TMPDIR/tmp.XXXXXXXXXX.rippchen
		-t       | --threads [value]        : threads - predicted default: $THREADS
		-mem     | --memory [value]         : amout of memory for creating bam slices and processing them in parallel instances
		                                      available: $MAXMEMORY
		                                      default: 30000 (allows for $MTHREADS instances)
		                                      NOTE: needs to be raised in case of GCThreads, HeapSize or OutOfMemory errors
		-resume  | --resume-from [value]    : resume from a specific pipeline step - see -dev|--devel
		-skip    | --skip [value,..]        : skip specific pipeline step(s) - see -dev|--devel, comma seperated
		-redo    | --redo [value,..]        : just rerun specific pipeline step(s) - see -dev|--devel, comma seperated
		-no-qual | --no-qualityanalysis     : disables read quality analysis
		-no-clip | --no-clipping            : disables removal of adapter sequences if -a|--adapter is used
		-no-trim | --no-trimming            : disables quality trimming
		-no-cor  | --no-correction          : disables majority based raw read error correction
		-no-rrm  | --no-rrnafilter          : disables rRNA filter
		-no-sege | --no-segemehl            : disables mapping by Segemehl
		-no-star | --no-star                : disables mapping by STAR
		-no-stats| --no-statistics          : disables fastq preprocessing and mapping statistics

		ALIGNMENT OPTIONS
		-d       | --distance               : maximum read alignment edit distance in % - default: 5
		-i       | --insertsize             : maximum allowed insert for aligning mate pairs - default: 200000
		-no-split| --no-split               : disable split read mapping
		-no-uniq | --no-uniqify             : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                : disables sorting alignments
		-no-idx  | --no-index               : disables indexing alignments
		-cmo     | --clipmateoverlaps       : enable read clipping in case of an overlap with its mate

		QUANTIFICATION OPTIONS
		-ql      | --quantifylevel          : switch to other feature type for quantification - default: exon
		                                      NOTE: quantifying using a different feature will break differential expression analysis
		-qt      | --quantifytag            : switch to other feature tag for quantification - default: gene_id
		-no-quant| --no-quantification      : disables per feature read quantification plus downstream analyses

		PEAK CALLING OPTIONS
		-ip      | --iptype [chip|rip]      : type of *IP-Seq experiment - default: chip
		-n1      | --normalfq1 [path,..]    : normal fastq input - single or first pair, comma seperated
		-n2      | --normalfq2 [path,..]    : normal fastq input - optional. second pair, comma seperated
		-t1      | --treatmentfq1 [path,..] : *IP-Seq fastq input - single or first pair, comma seperated
		-t2      | --treatmentfq2 [path,..] : *IP-Seq fastq input - optional. second pair, comma seperated
		-nr1     | --normalrepfq1 [path,..] : normal replicate fastq input - optional. single or first pair, comma seperated
		-nr2     | --normalrepfq2 [path,..] : normal replicate fastq input - optional. second pair, comma seperated
		-tr1     | --treatrepfq1 [path,..]  : *IP-Seq replicate fastq input - optional. single or first pair, comma seperated
		-tr2     | --treatrepfq2 [path,..]  : *IP-Seq replicate fastq input - optional. second pair, comma seperated
		-f       | --fragmentsize           : fragment size of sequenced mate pairs - default: 150
		-rx      | --regex                  : regex of read name identifier with grouped tile information - default: ^\S+:(\d+):(\d+):(\d+)\s*.*
		                                      NOTE: necessary for sucessful deduplication, if unavailable set to 'null'
		-no-rmd  | --no-removeduplicates    : disables removing duplicates - not recommended
		-no-macs | --no-macs                : disables peak calling by macs
		-no-gem  | --no-gem                 : disables peak calling by gem

		DIFFERENTIAL EXPRESSION ANALYSIS OPTIONS
		-1       | --fq1 [path,..]          : fastq input - single or first pair, comma seperated
		-2       | --fq2 [path,..]          : fastq input - optional. second pair, comma seperated
		-m       | --mapped [path,..]       : SAM/BAM input - comma seperated (replaces -1 and -2)
		                                      NOTE: alignment postprocessing steps can be disabled (see ALIGNMENT OPTIONS)
		-rmd     | --removeduplicates       : enable removing duplicates - not recommended
		-c       | --comparisons [path,..]  : experiment info file(s) for pairwise analyses according to column condition (primary factor)
		                                      format: 4 or more columns, seperated by tab or space(s)
		                                        sample   condition   [single-end|paired-end]   replicate   [factor1   factor2   ..]
		                                      NOTE: samples must match unique prefixes of input fastq basenames
		                                      
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
		-no-dsj  | --no-diffsplicejunctions : disables differential splice junction analysis
		-no-dea  | --no-diffexanalysis      : disables differential feature expression analysis plus downstream analyses
		-no-clust| --no-clustering          : disables downstream feature co-expression clustering
		-no-go   | --no-geneontology        : disables downstream gene ontology enrichment analysis

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
		rep   : pooling/generating replicates
		sort  : sorting and indexing of sam/bam files
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

		-1   | --fq1 | -n1 | --normalfq1) arg=true; nfq1=$2;;
		-2   | --fq2 | -n2 | --normalfq2) arg=true; nfq2=$2;;
		-m   | --mapped) arg=true; nmap=$2; noqual=true; notrim=true; nocor=true; norrm=true; nosege=true; nostar=true; notop=1; nobwa=true;;
		-nr1 | --normalrepfq1) arg=true; nrfq1=$2;;
		-nr2 | --normalrepfq2) arg=true; nrfq2=$2;;
		-t1  | --treatmentfq1) arg=true; tfq1=$2;;
		-t2  | --treatmentfq2) arg=true; tfq2=$2;;
		-tr1 | --treatrepfq1) arg=true; rfq1=$2;;
		-tr2 | --treatrepfq2) arg=true; rfq2=$2;;
		-rx  | --regex) arg=true; REGEX=$2;;
		-ip  | --iptype) arg=true; IPTYPE=$2;;
		-c   | --comparisons) arg=true; mapfile -t -d ',' COMPARISONS <<< $2; COMPARISONS[-1]="$(sed -r 's/\s*\n*$//' <<< "${COMPARISONS[-1]}")";;
		-a1  | --adapter1) arg=true; mapfile -t -d ',' ADAPTER1 <<< $2; ADAPTER1[-1]="$(sed -r 's/\s*\n*$//' <<< "${ADAPTER1[-1]}")";;
		-a2  | --adapter2) arg=true; mapfile -t -d ',' ADAPTER2 <<< $2; ADAPTER2[-1]="$(sed -r 's/\s*\n*$//' <<< "${ADAPTER2[-1]}")";;
		-d   | --distance) arg=true; DISTANCE=$2;;
		-f   | --fragmentsize) arg=true; FRAGMENTSIZE=$2;;
		-i   | --insertsize) arg=true; INSERTSIZE=$2;;
		-ql  | --quantifylevel) arg=true; QUANTIFYFLEVEL=$2;;
		-qt  | --quantifytag) arg=true; QUANTIFYTAG=$2;;
		#-cf  | --clusterfilter) arg=true; CLUSTERFILTER=$2;;
		#-dexnew | --dexnewprepare) arg=true; DEXSEQNEW=$2;;

	   	-resume | --resume-from)
			arg=true
			local enable=false
			# don't Smd5, Sslice !
			for s in qual trim clip cor rrm sege star uniq rep sort rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
				eval "\${S$s:=true}" # unless S$s already set to false by -redo, do skip
				$enable || [[ "$2" == "$s" ]] && {
					enable=true
					eval "S$s=false"
				}
			done
		;;
	    -skip | --skip)
			arg=true
			mapfile -d ',' -t <<< $2
			for x in ${MAPFILE[@]}; do # do not quote!! "MAPFILE[@]" appends newline to last element
				for s in md5 qual trim clip cor rrm sege star uniq rep sort slice rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
					[[ "$x" == "$s" ]] && eval "S$s=true"
				done
			done
		;;
		-redo | --redo)
			arg=true
			# don't Smd5, Sslice !
			for s in qual trim clip cor rrm sege star uniq rep sort rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
				eval "\${S$s:=true}" # unless S$s alredy set to false by -resume, do skip
			done
			mapfile -d ',' -t <<< $2
			for x in ${MAPFILE[@]}; do # do not quote!! "MAPFILE[@]" appends newline to last element
				for s in qual trim clip cor rrm sege star uniq rep sort rmd cmo idx stats macs gem quant tpm dsj dea join clust go; do
					[[ "$x" == "$s" ]] && eval "S$s=false"
				done
			done
		;;

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
		-no-dea   | --no-diffsplicejunctions) nodsj=true;;
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
