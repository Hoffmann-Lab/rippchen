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
[[ $MTHREADS -eq 0 ]] && die "too less memory available ($MAXMEMORY)"
VERBOSITY=0
OUTDIR=$PWD/results
TMPDIR=$OUTDIR
REGEX='\S+:(\d+):(\d+):(\d+)\s*.*'
DISTANCE=5
FRAGMENTSIZE=200
FASTQ1=() # all idx of FASTQ1[.] are equal to MAPPED[.]
FASTQ2=()
MAPPED=()
nidx=() #normal idx
nridx=() #normal replicate idx
tidx=() #treatment idx
ridx=() #treatment replicate idx
pidx=() #pool (2x0.5 pseudoppol and 2x1 fullpool) idx


options::parse "$@" || die "parameterization issue"


mkdir -p $OUTDIR || die "cannot access $OUTDIR"
OUTDIR=$(readlink -e $OUTDIR)
[[ ! $LOG ]] && LOG=$OUTDIR/run.log
mkdir -p $(dirname $LOG) || die "cannot access $LOG"
printf '' > $LOG || die "cannot access $LOG"

[[ $PREVIOUSTMPDIR ]] && {
	TMPDIR=$PREVIOUSTMPDIR
	mkdir -p $TMPDIR || die "cannot access $TMPDIR"
	TMPDIR=$(readlink -e $TMPDIR)
} || {
	mkdir -p $TMPDIR || die "cannot access $TMPDIR"
	TMPDIR=$(readlink -e $TMPDIR)
	TMPDIR=$(mktemp -d -p $TMPDIR rippchen.XXXXXXXXXX) || die "cannot access $TMPDIR"
}

${INDEX:=false} || {
	[[ ! $nfq1 ]] && [[ ! $tfq1 ]] && [[ ! $nmap ]] && die "fastq or sam/bam file input missing"
}

[[ ! $nfq2 ]] && [[ "$nocmo" == "false" ]] && {
	commander::warn "second mate fastq file missing. proceeding without mate overlap clipping"
	nocmo=true
}

[[ $GENOME ]] && {
	readlink -e $GENOME | file -f - | grep -qF ASCII || die "genome file does not exists or is compressed $GENOME"
	[[ ! -s $GENOME.md5.sh ]] && cp $(dirname $(readlink -e $0))/bashbone/lib/md5.sh $GENOME.md5.sh
	source $GENOME.md5.sh
} || {
	${INDEX:=false} && die "genome file missing"
	commander::warn "genome file missing. proceeding without mapping"
	Smd5=true
	nosege=true
	nostar=true
}

[[ $GTF ]] && {
	readlink -e $GTF | file -f - | grep -qF ASCII || die "annotation file does not exists or is compressed $GTF"
} || {
	readlink -e $GENOME.gtf | file -f - | grep -qF ASCII && {
		GTF=$GENOME.gtf
	} || {
		${INDEX:=false} && die "annotation file missing"
		commander::warn "gtf file missing. proceeding without quantification"
		noquant=true
	}
}


checkfile(){
	declare -n _idx=$2 _arr=$3
	local f ifs
	[[ $(readlink -e $1 | file -f - | grep ASCII) && -e $(readlink -e $(head -1 $1)) ]] && {
		ifs=$IFS
		unset IFS
		while read -r f; do
			readlink -e $f &> /dev/null || die "$4 $f"
			_arr[((++i))]=$f
			_idx+=($i)
		done < $1
		IFS=$ifs
	} || {
		f=$1
		readlink -e $f &> /dev/null || die "$4 $f"
		_arr[((++i))]=$f
		_idx+=($i)
	}
}

IFS=','
i=-1
if [[ ! $nmap ]]; then
	for f in $nfq1; do
		checkfile $f nidx FASTQ1 "single or first mate fastq file does not exist"
	done
	for f in $nrfq1; do
		checkfile $f nridx FASTQ1 "single or first mate normal replicate fastq file does not exists"
	done
	for f in $tfq1; do
		#idx depends on available replicates - i.e. trigger pooling or make pseudo-replicates
		checkfile $f $([[ $rfq1 ]] && echo tidx || echo pidx) FASTQ1 "single or first mate normal replicate fastq file does not exists"
	done
	for f in $rfq1; do
		checkfile $f ridx FASTQ1 "single or first mate treatment replicate fastq file does not exists"
	done
	i=-1
	for f in {$nfq2,$nrfq2,$tfq2,$rfq2}; do
		checkfile $f foo FASTQ2 "second mate fastq file does not exists"
	done
	[[ $FASTQ2 ]] && [[ ${#FASTQ1[@]} -ne ${#FASTQ2[@]} ]] && die "unequal number of mate pairs"
else
	for f in $nmap; do
		checkfile $f nidx MAPPED "alignment file does not exists"
	done
	for f in $nrmap; do
		checkfile $f nridx MAPPED "normal replicate alignment file does not exists"
	done
	[[ $nridx ]] && [[ ${#nidx[@]} -ne ${#nridx[@]} ]] && die "unequal number of normal replicates"
	for f in $tmap; do
		checkfile $f pidx MAPPED "treatment alignment file does not exists"
	done
	for f in $rmap; do
		checkfile $f ridx MAPPED "treatment replicate alignment file does not exists"
	done
	[[ $ridx ]] && {
		tidx+=("${pidx[@]}")
		pidx=()
		[[ ${#ridx[@]} -ne ${#tidx[@]} ]] && die "unequal number of treatment replicates"
	}
fi
unset IFS

[[ ${#nidx[@]} -lt 2 ]] && [[ "$noclust" == "false" ]] && {
	commander::warn "too few samples. proceeding without clustering"
	noclust=true
}

progress::log -v $VERBOSITY -o $LOG
commander::printinfo "rippchen $VERSION utilizing bashbone $BASHBONEVERSION started with command: $CMD" >> $LOG
commander::printinfo "temporary files go to: $HOSTNAME:$TMPDIR" >> $LOG

${INDEX:=false} && {
	pipeline::index 2> >(tee -ai $LOG >&2) >> $LOG || die
} || {
	if [[ $tfq1 || $tmap ]]; then
		${RIPSEQ:=false} || nosplit=true
		pipeline::callpeak 2> >(tee -ai $LOG >&2) >> $LOG
		[[ $? -gt 0 ]] && die
	else
		pipeline::dea 2> >(tee -ai $LOG >&2) >> $LOG
		[[ $? -gt 0 ]] && die
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
