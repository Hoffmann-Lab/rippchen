#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'die' INT TERM
trap 'sleep 1; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT
shopt -s extglob

die() {
	unset CLEANUP
	echo -ne "\e[0;31m"
	echo ":ERROR: $*" >&2
	echo -ne "\e[m"
	exit 1
}

cleanup() {
	if [[ $CLEANUP ]]; then
		local b e
		for f in "${FASTQ1[@]}"; do
			helper::basename -f "$f" -o b -e e
			f=$b
			[[ -e $TMPDIR ]] && find $TMPDIR -type f -name "$f*" -exec rm -f {} \;
			if [[ -e $OUTDIR ]]; then
				find $OUTDIR -type f -name "$f*.all" -exec rm -f {} \;
				find $OUTDIR -type f -name "$f*.sorted.bam" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .sorted.bam).bam' \;
				find $OUTDIR -type f -name "$f*.*.gz" -exec bash -c '[[ -s {} ]] && rm -f $(dirname {})/$(basename {} .gz)' \;
			fi
		done
	fi
}

[[ ! $OSTYPE =~ linux ]] && die "unsupported operating system"
[[ ${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4 ]] || [[ ${BASH_VERSINFO[0]} -gt 4 ]] || die "requieres bash version 4.4 or above"
[[ ! $RIPPCHEN ]] && die "cannot find installation. please run setup and/or do: export RIPPCHEN=/path/to/install/dir"
INSDIR=$RIPPCHEN
for f in {$INSDIR/latest/bashbone/lib/*.sh,$INSDIR/latest/rippchen/lib/*.sh}; do
	source $f
done
configure::environment -i $INSDIR

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
pidx=() #pool (2x0.5) idx 
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
	Sslice=false
	mkdir -p $TMPDIR || die "cannot access $TMPDIR"
	TMPDIR=$(readlink -e $TMPDIR)
	TMPDIR=$(mktemp -p $TMPDIR -d --suffix=.rippchen) || die "cannot access $TMPDIR"
fi


[[ ! $LOG ]] && LOG=$OUTDIR/run.log
[[ MTHREADS=$[MAXMEMORY/MEMORY] -gt $THREADS ]] && MTHREADS=$THREADS
[[ $MTHREADS -eq 0 ]] && die "too less memory available ($MAXMEMORY)"
[[ ! $nfq1 ]] && [[ ! $tfq1 ]] && [[ ! $nmap ]] && die "fastq file input missing - call "$(basename $0)" -h for help"
if [[ $GENOME ]]; then
	readlink -e $GENOME | file -f - | grep -qF ASCII || die "genome file does not exists or is compressed $GENOME"
else
	commander::warn "proceeding without genome file"
	Smd5=true
	nosege=true
	nostar=true
fi
if [[ $GTF ]]; then
	readlink -e $GTF | file -f - | grep -qF ASCII || die "annotation file does not exists or is compressed $GTF"
else
	readlink -e $GENOME.gtf | file -f - | grep -qF ASCII && {
		GTF=$GENOME.gtf
	} || {
		if [[ ! $tfq1 ]]; then
			commander::warn "proceeding without gtf file"
			noquant=true
		else
			[[ $comp ]] && die "annotation file needed"
		fi
	}
fi


i=-1
IFS=','
for f in $nfq1; do
	readlink -e $f &> /dev/null || die "single or first mate fastq file does not exists $f"
	FASTQ1[((++i))]=$f
	nidx+=($i)
done
for f in $nrfq1; do
	readlink -e $f &> /dev/null || die "single or first mate replicate fastq file does not exists $f"
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
	[[ $comp ]] && [[ $(basename ${FASTQ2[$i]} | cut -d '.' -f 1) != ${fqbaseprefix[$i]} ]] && die "'.'-separated basenames do not match"
done
for f in $nmap; do
	readlink -e $f &> /dev/null || die "alignment file does not exists $f"
	MAPPED+=($f)
done
unset IFS


commander::print "rippchen v$version started with command: $CMD" > $LOG || die "cannot access $LOG"
commander::print "temporary files go to: $HOST:$TMPDIR" >> $LOG
progress::log -v $VERBOSITY -o $LOG

if [[ $tfq1 ]]; then
	[[ $IPTYPE == 'chip' ]] && nosplit=true || IPTYPE='rip'
	pipeline::callpeak >> $LOG 2> >(tee -a $LOG >&2) || die
else
	normd=true
	pipeline::dea >> $LOG 2> >(tee -a $LOG >&2) || die
fi

commander::print "success" >> $LOG
exit 0
