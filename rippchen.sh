#! /usr/bin/env bash
# (c) Konstantin Riege

die() {
	echo ":ERROR: $*" >&2
	exit 1
}

cleanup() {
	[[ -e $TMPDIR ]] && {
		find $TMPDIR -type f -name "cleanup.*" -exec rm -f {} \;
		find $TMPDIR -depth -type d -name "cleanup.*" -exec rm -rf {} \;
	}
	[[ $1 -eq 0 ]] && ${CLEANUP:=false} && {
		[[ -e $TMPDIR ]] && {
			find $TMPDIR -type f -exec rm -f {} \;
			find $TMPDIR -depth -type d -exec rm -rf {} \;
			rm -rf $TMPDIR
		}
		[[ -e $OUTDIR ]] && {
			for f in "${FASTQ1[@]}"; do
				readlink -e "$f" | file -f - | grep -qE '(gzip|bzip)' && b=$(basename $f | rev | cut -d '.' -f 3- | rev) || b=$(basename $f | rev | cut -d '.' -f 2- | rev)
				find $OUTDIR -depth -type d -name "$b*._STAR*" -exec rm -rf {} \;
				find $OUTDIR -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .sorted.bam).bam' \;
				find $OUTDIR -type f -name "$b*.*.gz" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .gz)' \;
			done
			local b
			for f in "${MAPPED[@]}"; do
				b=$(basename $f | rev | cut -d '.' -f 2- | rev)
				find $OUTDIR -depth -type d -name "$b*._STAR*" -exec rm -rf {} \;
				find $OUTDIR -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .sorted.bam).bam' \;
			done
		}
	}
}

# defines INSDIR and by sourcing bashbone it defines BASHBONEVERSION variable as well
source $(dirname $(readlink -e $0))/activate.sh -c true || die

trap 'configure::exit -p $$ -f cleanup $?' EXIT
trap 'die "killed"' INT TERM

VERSION=$version
CMD="$(basename $0) $*"
THREADS=$(grep -cF processor /proc/cpuinfo)
MAXMEMORY=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
MEMORY=30000
[[ MTHREADS=$((MAXMEMORY/MEMORY)) -gt $THREADS ]] && MTHREADS=$THREADS
VERBOSITY=0
OUTDIR=$PWD/results
TMPDIR=$OUTDIR
REGEX='\S+:(\d+):(\d+):(\d+)\s*.*'
DISTANCE=5
FRAGMENTSIZE=150
IPTYPE='chip'
# all idx of FASTQ1[.] are equal to MAPPER[.]
nidx=() #normal idx 
nridx=() #normal replicate idx 
tidx=() #treatment idx
ridx=() #treatment replicate idx
pidx=() #pool (2x0.5 pseudoppol and 2x1 fullpool) idx 
FASTQ1=()
FASTQ2=()
MAPPED=()

options::parse "$@" || die "parameterization issue"

mkdir -p $OUTDIR || die "cannot access $OUTDIR"
OUTDIR=$(readlink -e $OUTDIR)
if [[ $PREVIOUSTMPDIR ]]; then
	TMPDIR=$PREVIOUSTMPDIR
	mkdir -p $TMPDIR || die "cannot access $TMPDIR"
	TMPDIR=$(readlink -e $TMPDIR)
else
	mkdir -p $TMPDIR || die "cannot access $TMPDIR"
	TMPDIR=$(readlink -e $TMPDIR)
	TMPDIR=$(mktemp -d -p $TMPDIR rippchen.XXXXXXXXXX) || die "cannot access $TMPDIR"
fi


[[ ! $LOG ]] && LOG=$OUTDIR/run.log
[[ MTHREADS=$((MAXMEMORY/MEMORY)) -gt $THREADS ]] && MTHREADS=$THREADS
[[ $MTHREADS -eq 0 ]] && die "too less memory available ($MAXMEMORY)"
${INDEX:=false} || {
	[[ ! $nfq1 ]] && [[ ! $tfq1 ]] && [[ ! $nmap ]] && die "fastq or sam/bam file input missing"
}
[[ ! $nfq2 ]] && {
	[[ "$nocmo" == "false" ]] && {
		commander::warn "no second mate fastq file given - proceeding without mate overlap clipping"
		nocmo=true
	}
}
if [[ $GENOME ]]; then
	readlink -e $GENOME | file -f - | grep -qF ASCII || die "genome file does not exists or is compressed $GENOME"
else
	${INDEX:=false} && die "genome file missing"
	commander::warn "proceeding without genome file"
	nosege=true
	nostar=true
fi
if [[ $GTF ]]; then
	readlink -e $GTF | file -f - | grep -qF ASCII || die "annotation file does not exists or is compressed $GTF"
else
	readlink -e $GENOME.gtf | file -f - | grep -qF ASCII && {
		GTF=$GENOME.gtf
	} || {
		${INDEX:=false} && die "annotation file missing"
		commander::warn "proceeding without gtf file"
		noquant=true
	}
fi


IFS=','
i=-1
if [[ ! $nmap ]]; then
	for f in $nfq1; do
		readlink -e $f &> /dev/null || die "single or first mate fastq file does not exists $f"
		FASTQ1[((++i))]=$f
		nidx+=($i)
	done
	for f in $nrfq1; do
		readlink -e $f &> /dev/null || die "single or first mate normal replicate fastq file does not exists $f"
		FASTQ1[((++i))]=$f
		nridx+=($i)
	done
	for f in $tfq1; do
		readlink -e $f &> /dev/null || die "single or first mate treatment fastq file does not exists $f"
		FASTQ1[((++i))]=$f
		[[ $rfq1 ]] && tidx+=($i) || pidx+=($i) #necessary for pooling, make pseudo-replicates respectively
	done
	for f in $rfq1; do
		readlink -e $f &> /dev/null  || die "single or first mate treatment replicate fastq file does not exists $f"
		FASTQ1[((++i))]=$f
		ridx+=($i)
	done
	i=-1
	for f in {$nfq2,$nrfq2,$tfq2,$rfq2}; do
		readlink -e $f &> /dev/null || die "second mate fastq file does not exists $f"
		FASTQ2[((++i))]=$f
	done
	[[ $FASTQ2 ]] && [[ ${#FASTQ1[@]} -ne ${#FASTQ2[@]} ]] && die "unequal number of mate pairs"
else
	for f in $nmap; do
		readlink -e $f &> /dev/null || die "alignment file does not exists $f"
		MAPPED[((++i))]=$f
		nidx+=($i)
	done
	for f in $nrmap; do
		readlink -e $f &> /dev/null || die "normal replicate alignment file does not exists $f"
		MAPPED[((++i))]=$f
		nridx+=($i)
	done
	[[ $nridx ]] && [[ ${#nidx[@]} -ne ${#nridx[@]} ]] && die "unequal number of normal replicates"
	for f in $tmap; do
		readlink -e $f &> /dev/null || die "treatment alignment file does not exists $f"
		MAPPED[((++i))]=$f
		pidx+=($i)
	done
	for f in $rmap; do
		readlink -e $f &> /dev/null || die "treatment replicate alignment file does not exists $f"
		MAPPED[((++i))]=$f
		ridx+=($i)
	done
	[[ $ridx ]] && {
		tidx+=("${pidx[@]}")
		pidx=()
		[[ ${#ridx[@]} -ne ${#tidx[@]} ]] && die "unequal number of treatment replicates"
	}
fi
unset IFS

printf '' > $LOG || die "cannot access $LOG"
progress::log -v $VERBOSITY -o $LOG
commander::printinfo "rippchen $VERSION utilizing bashbone $BASHBONEVERSION started with command: $CMD" >> $LOG
commander::printinfo "temporary files go to: $HOSTNAME:$TMPDIR" >> $LOG

${Smd5:=false} || {
	[[ ! -s $GENOME.md5.sh ]] && cp $(dirname $(readlink -e $0))/bashbone/lib/md5.sh $GENOME.md5.sh
	source $GENOME.md5.sh
}
${INDEX:=false} && {
	pipeline::index 2> >(tee -ai $LOG >&2) >> $LOG || die
} || {
	if [[ $tfq1 || $tmap ]]; then
		${RIPSEQ:=false} || nosplit=true
		pipeline::callpeak 2> >(tee -ai $LOG >&2) >> $LOG || die
	else
		pipeline::dea 2> >(tee -ai $LOG >&2) >> $LOG || die
	fi
}
${Smd5:=false} || {
	commander::printinfo "finally updating genome and annotation md5 sums" >> $LOG
	thismd5genome=$(md5sum $GENOME | cut -d ' ' -f 1)
	[[ "$md5genome" != "$thismd5genome" ]] && sed -i "s/md5genome=.*/md5genome=$thismd5genome/" $GENOME.md5.sh
	thismd5gtf=$(md5sum $GTF | cut -d ' ' -f 1)
	[[ "$md5gtf" != "$thismd5gtf" ]] && sed -i "s/md5gtf=.*/md5gtf=$thismd5gtf/" $GENOME.md5.sh
}

commander::printinfo "success" >> $LOG
exit 0
