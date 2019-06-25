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
		acquire a taste for peak calling from *IP-Seq experiments or for differential gene expression- and ontology analysis from RNA-Seq data

		VERSION
		$version

		SYNOPSIS DGE ANALYSIS
		$(basename $0) -1 ctr1.fq,ctr2.fq,treat1.fq,treat2.fq -g genome.fa -c info

		SYNOPSIS PEAK CALLING
		$(basename $0) -1 ctrA.fq,ctrB.fq -t1 treatA1.fq,treatB1.fq -r1 treatA2.fq,treatB2.fq -g genome.fa

		SYNOPSIS PEAK CALLING AND DGE ANALYSIS
		$(basename $0) -1 ctrA.fq,ctrB.fq -t1 treatA1.fq,treatB1.fq -r1 treatA2.fq,treatB2.fq -g genome.fa -c info

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
		-a       | --adapter [string,..]    : adapter sequence(s), comma seperated - optional
		-o       | --out [path]             : output directory - default: $OUTDIR
		-l       | --log [path]             : output directory - default: $OUTDIR/run.log
		-tmp     | --tmp                    : temporary directory - default: $TMPDIR/rippchen_tmp
		-t       | --threads [value]        : threads - predicted default: $THREADS
		-mem     | --memory [value]         : amout of memory for creating bam slices and processing them in parallel instances
		                                      available: $MAXMEMORY
		                                      default: 30000 (allows for $MTHREADS instances)
		                                      NOTE: needs to be raised in case of GCThreads, HeapSize or OutOfMemory errors

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
		                                      NOTE: necessary for sucessfully duplicates removal
		-no-rmd  | --no-removeduplicates    : disables removing duplicates - not recommended
		-no-macs | --no-macs                : disables peak calling by macs
		-no-gem  | --no-gem                 : disables peak calling by gem

		DIFF GENE EXPR ANALYSIS OPTIONS
		-1       | --fq1 [path,..]          : fastq input - single or first pair, comma seperated
		-2       | --fq2 [path,..]          : fastq input - optional. second pair, comma seperated
		-m       | --mapped [path,..]       : SAM/BAM input - comma seperated (replaces -1 and -2)
		                                      NOTE: alignment postprocessing steps can be disabled (see ALIGNMENT OPTIONS)
		-c       | --comparisons [path,..]  : experiments info file(s) for pairwise comparisons (according to column two)
		                                      - triggers differential gene expression analysis
		                                      - requires -gtf|--gtf (see BASIC OPTIONS)
		                                      - format: 4 or 5 tab-seperated columns (5 in case of paired analysis)
		                                        common_basename   experiment   [single-end|paired-end]   Nreplicate   pairs
		                                      - fastq files needs to have unique basenames up to the first '.' in column 1
		                                      - example for input A1.R1.fq, A1.R2.fq, A2.fq, B.R1.fq, B.R2.fq, C.fq:
		                                        A1   experimentA   paired-end   N1   patient1
		                                        A2   experimentA   single-end   N2   patient2
		                                        B    experimentB   paired-end   N1   patient1
		                                        C    experimentC   single-end   N1   patient2
		-ql      | --quantifylevel          : switch to other feature type for quantification - default: exon
		                                      NOTE: quantifying using a different feature will break differential gene expression analysis
		-qt      | --quantifytag            : switch to other feature tag for quantification - default: gene_id
		-no-quant| --no-quantification      : disables per feature read quantification plus downstream analyses
		-no-dea  | --no-diffexanalysis      : disables differential feature expression analysis plus downstream analyses
		-no-clust| --no-clustering          : disables feature co-expression clustering for samples defined via -c
		-no-go   | --no-geneontology        : disables gene ontology enrichment analysis

		GENERAL OPTIONS
		-resume  | --resume-from [value]    : resume from a specific pipeline step - see -dev|--devel
		-skip    | --skip [value,..]        : skip specific pipeline step(s) - see -dev|--devel, comma seperated
		-redo    | --redo [value]           : just rerun a specific pipeline step - see -dev|--devel, comma seperated
		-no-qual | --no-qualityanalysis     : disables quality analysis
		-no-clip | --no-clipping            : disables removal of adapter sequences if -a|--adapter is used
		-no-trim | --no-trimming            : disables quality trimming
		-no-cor  | --no-correction          : disables majority based raw read error correction
		-no-rrm  | --no-rrnafilter          : disables rRNA filter
		-no-sege | --no-segemehl            : disables mapping by Segemehl
		-no-star | --no-star                : disables mapping by STAR
		-no-stats| --no-statistics          : disables fastq preprocessing statistics

		ALIGNMENT OPTIONS
		-d       | --distance               : maximum read alignment edit distance in % - default: 5
		-i       | --insertsize             : maximum allowed insert for aligning mate pairs - default: 200000
		-no-split| --no-split               : disable split read mapping
		-no-uniq | --no-uniqify             : disables extraction of properly paired and uniquely mapped reads
		-no-sort | --no-sort                : disables sorting alignments

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
		qual  : quality analysis
		clip  : adapter clipping
		trim  : trimming
		cor   : raw read correction
		rrm   : rRNA filtering
		stats : proprocessing statistics
		sege  : Segemehl mapping
		star  : STAR mapping
		uniq  : extraction of properly paired and uniquely mapped reads
		rep   : pooling/generating replicates
		sort  : sorting and indexing of sam/bam files
		slice : slicing bams during rmd for parallelization
		rmd   : removing duplicates
		idx   : intermediate and final bam indexing
		stats : fastq preprocessing and mapping statistics
		macs  : peak calling by macs
		gem   : peak calling by gem
		quant : read quantification
		tpm   : TPM calculation
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
		-a   | --adapter) arg=true; mapfile -t -d ',' ADAPTER <<< $2; ADAPTER[-1]="$(sed -r 's/\s*\n*$//' <<< "${ADAPTER[-1]}")";;
		-d   | --distance) arg=true; DISTANCE=$2;;
		-f   | --fragmentsize) arg=true; FRAGMENTSIZE=$2;;
		-i   | --insertsize) arg=true; INSERTSIZE=$2;;
		-ql  | --quantifylevel) arg=true; QUANTIFYFLEVEL=$2;;
		-qt  | --quantifytag) arg=true; QUANTIFYTAG=$2;;
		#-cf  | --clusterfilter) arg=true; CLUSTERFILTER=$2;;
		#-dexnew | --dexnewprepare) arg=true; DEXSEQNEW=$2;;

	   	-resume | --resume-from)
			arg=true
			# don't Smd5, Sslice !
			for s in qual clip trim cor rrm stats sege star uniq rep sort rmd idx macs gem quant tpm dea join clust go; do
				[[ "$2" == "$s" ]] && break
				eval "S$s=true"
			done
		;;
	    -skip   | --skip) 
			arg=true
			mapfile -d ',' -t <<< $2
			for x in ${MAPFILE[@]}; do # do not quote!! "MAPFILE[@]" appends newline to last element
				for s in md5 qual clip trim cor rrm stats sege star uniq rep sort slice rmd idx macs gem quant tpm dea join clust go; do
					[[ "$x" == "$s" ]] && eval "S$s=true"
				done
			done
		;;
		-redo | --redo)
			arg=true
			# don't Smd5, Sslice !
			for s in qual clip trim cor rrm stats sege star uniq rep sort rmd idx macs gem quant tpm dea join clust go; do
				[[ "$2" == "$s" ]] && continue
				eval "S$s=true"
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
		-no-stats | --no-statistics) nostats=true;;
		-no-rmd   | --no-removeduplicates) normd=true;;
		-no-macs  | --no-macs) nomacs=true;;
		-no-gem   | --no-gem) nogem=true;;
		-no-quant | --no-quantification) noquant=true; nodea=true; noclust=true; nogo=true;;
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
