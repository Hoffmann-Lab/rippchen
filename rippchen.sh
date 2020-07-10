#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'die' INT TERM
trap 'sleep 1; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT
shopt -s extglob
shopt -s expand_aliases

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

[[ ! $RIPPCHEN ]] && die "cannot find installation. please run setup and/or do: export RIPPCHEN=/path/to/install/dir"
INSDIR=$RIPPCHEN
source $(dirname $(readlink -e $0))/bashbone/activate.sh -i $RIPPCHEN -c true || die
BASHBONEVERSION=$version
for f in $(dirname $(readlink -e $0))/lib/*.sh; do
	source $f || die "unexpected error in source code - please contact developer"
done
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
${INDEX:=false} || {
	[[ ! $nfq1 ]] && [[ ! $tfq1 ]] && [[ ! $nmap ]] && die "fastq file input missing"
}
[[ ! $nfq2 ]] && {
	[[ "$nocomo" == "false" ]] && {
		commander::warn "no second mate fastq file given - proceeding without mate overlap clipping"
		nocmo=true
	}
}
if [[ $GENOME ]]; then
	readlink -e $GENOME | file -f - | grep -qF ASCII || die "genome file does not exists or is compressed $GENOME"
else
	${INDEX:=false} && die "genome file missing"
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
		${INDEX:=false} && die "annotation file missing"
		commander::warn "proceeding without gtf file"
		noquant=true
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
done
for f in $nmap; do
	readlink -e $f &> /dev/null || die "alignment file does not exists $f"
	MAPPED+=($f)
done
unset IFS

echo > $LOG || die "cannot access $LOG"
progress::log -v $VERBOSITY -o $LOG
commander::print "rippchen $VERSION utilizing bashbone $BASHBONEVERSION started with command: $CMD" > $LOG 
commander::print "temporary files go to: $HOSTNAME:$TMPDIR" >> $LOG

${Smd5:=false} || {
	[[ ! -s $GENOME.md5.sh ]] && cp $(dirname $(readlink -e $0))/bashbone/lib/md5.sh $GENOME.md5.sh
	source $GENOME.md5.sh
}
${INDEX:=false} && {
	pipeline::index >> $LOG 2> >(tee -a $LOG >&2) || die
} || {
	if [[ $tfq1 ]]; then
		[[ $IPTYPE == 'chip' ]] && nosplit=true || IPTYPE='rip'
		pipeline::callpeak >> $LOG 2> >(tee -a $LOG >&2) || die
	else
		pipeline::dea >> $LOG 2> >(tee -a $LOG >&2) || die
	fi
}
${Smd5:=false} || {
	commander::print "finally updating genome and annotation md5 sums" >> $LOG
	thismd5genome=$(md5sum $GENOME | cut -d ' ' -f 1)
	[[ "$md5genome" != "$thismd5genome" ]] && sed -i "s/md5genome=.*/md5genome=$thismd5genome/" $GENOME.md5.sh
	thismd5gtf=$(md5sum $GTF | cut -d ' ' -f 1)
	[[ "$md5gtf" != "$thismd5gtf" ]] && sed -i "s/md5gtf=.*/md5gtf=$thismd5gtf/" $GENOME.md5.sh
}

commander::print "success" >> $LOG
exit 0


