#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$(readlink -e "$0")")/activate.sh" -c true -x cleanup -a "$@" || exit 1

cleanup() {
	[[ -e "$LOG" ]] && {
		commander::printinfo "date: $(date)" | tee -ia "$LOG"
		[[ $1 -eq 0 ]] && commander::printinfo "success" | tee -ia "$LOG" || commander::printinfo "failed" | tee -ia "$LOG"
	}
	[[ -e "$CLEANUP_TMPDIR" ]] && {
		find -L "$CLEANUP_TMPDIR" -type f -name "cleanup.*" -exec rm -f "{}" \; &> /dev/null || true
		find -L "$CLEANUP_TMPDIR" -depth -type d -name "cleanup.*" -exec rm -rf "{}" \; &> /dev/null || true
	}
	[[ $1 -eq 0 ]] && ${CLEANUP:=false} || ${FORCECLEANUP:=false} && {
		[[ -e "$CLEANUP_TMPDIR" ]] && {
			find -L "$CLEANUP_TMPDIR" -type f -exec rm -f "{}" \; &> /dev/null || true
			find -L "$CLEANUP_TMPDIR" -depth -type d -exec rm -rf "{}" \; &> /dev/null || true
			rm -rf "$CLEANUP_TMPDIR"
		}
		[[ -e "$OUTDIR" ]] && {
			local b f
			for f in "${FASTQ1[@]}"; do
				readlink -e "$f" | file -f - | grep -qE '(gzip|bzip)' && b=$(basename "$f" | rev | cut -d '.' -f 3- | rev) || b=$(basename "$f" | rev | cut -d '.' -f 2- | rev)
				find -L "$OUTDIR" -depth -type d -name "$b*._STAR*" -exec rm -rf {} \; &> /dev/null || true
				find -L "$OUTDIR" -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s "{}" ]] && rm -f "$(dirname "{}")/$(basename "{}" .sorted.bam).bam"' \; &> /dev/null || true
				find -L "$OUTDIR" -type f -name "$b*.*.gz" -exec bash -c '[[ -s "{}" ]] && rm -f "$(dirname "{}")/$(basename "{}" .gz)"' \; &> /dev/null || true
			done
			for f in "${MAPPED[@]}"; do
				b=$(basename "$f" | rev | cut -d '.' -f 2- | rev)
				find -L "$OUTDIR" -depth -type d -name "$b*._STAR*" -exec rm -rf "{}" \; &> /dev/null || true
				find -L "$OUTDIR" -type f -name "$b*.sorted.bam" -exec bash -c '[[ -s {} ]] && rm -f "$(dirname "{}")/$(basename "{}" .sorted.bam).bam"' \; &> /dev/null || true
				find -L "$OUTDIR" -type f -name "$b*.*.gz" -exec bash -c '[[ -s {} ]] && rm -f "$(dirname "{}")/$(basename "{}" .gz)"' \; &> /dev/null || true
			done
		}
	}
	return 0
}

VERSION=$version
CMD="$(basename "$0") $*"
THREADS=$(grep -cF processor /proc/cpuinfo)
MAXMEMORY=$(grep -F MemTotal /proc/meminfo | awk '{printf("%d",$2*0.95/1024)}')
MEMORY=10000
[[ MTHREADS=$((MAXMEMORY/MEMORY)) -gt $THREADS ]] && MTHREADS=$THREADS
BASHBONE_ERROR="too less memory available ($MAXMEMORY)"
[[ $MTHREADS -eq 0 ]] && false
VERBOSITY=0
OUTDIR="$PWD/results"
TMPDIR="${TMPDIR:-$OUTDIR}"
DISTANCE=5
FRAGMENTSIZE=200
FASTQ1=() # all idx of FASTQ1[.] are equal to MAPPED[.]
FASTQ2=()
FASTQ3=()
MAPPED=()
nidx=() #normal idx
nridx=() #normal replicate idx
tidx=() #treatment idx
ridx=() #treatment replicate idx
pidx=() #pool (2x0.5 pseudoppol and 2x1 fullpool) idx

BASHBONE_ERROR="parameterization issue"
options::parse "$@"

BASHBONE_ERROR="cannot access $OUTDIR"
mkdir -p "$OUTDIR"
OUTDIR="$(readlink -e "$OUTDIR")"
[[ ! $LOG ]] && LOG="$OUTDIR/run.log"
BASHBONE_ERROR="cannot access $LOG"
mkdir -p "$(dirname "$LOG")"

BASHBONE_ERROR="cannot access $TMPDIR"
if [[ $PREVIOUSTMPDIR ]]; then
	TMPDIR="$(realpath -se "$PREVIOUSTMPDIR")"
else
	mkdir -p "$TMPDIR"
	TMPDIR="$(realpath -se "$TMPDIR")"
	TMPDIR="$(mktemp -d -p "$TMPDIR" rippchen.XXXXXXXXXX)"
fi
CLEANUP_TMPDIR="$TMPDIR"

${INDEX:=false} || {
	BASHBONE_ERROR="fastq or sam/bam file input missing"
	[[ ! $nfq1 && ! $tfq1 && ! $nmap && ! $tmap ]] && false
}

[[ $nfq2 ]] || {
	[[ "$nocmo" == "false" ]] && commander::warn "second mate fastq file missing. proceeding without mate overlap clipping"
	nocmo=true
}

if [[ $GENOME ]]; then
	BASHBONE_ERROR="genome file does not exists or is compressed $GENOME"
	readlink -e "$GENOME" | file -f - | grep -qF ASCII
	[[ ! -s "$GENOME.md5.sh" ]] && cp "$(dirname "$(readlink -e "$0")")/bashbone/lib/md5.sh" "$GENOME.md5.sh"
	source "$GENOME.md5.sh"
else
	BASHBONE_ERROR="genome file missing"
	! ${INDEX:=false}
	commander::warn "genome file missing. proceeding without mapping"
	nosege=true
	nostar=true
	nobwa=true
	Smd5=true
fi

if [[ $GTF ]]; then
	BASHBONE_ERROR="annotation file does not exists or is compressed $GTF"
	readlink -e "$GTF" | file -f - | grep -qF ASCII
else
	# readlink -e "$GENOME.gtf" | file -f - | grep -qF ASCII && {
	# 	GTF="$GENOME.gtf"
	# } || {
		if ! ${BISULFITE:=false}; then # does not require gtf, even for indexing
			if ${INDEX:=false}; then
				commander::warn "gtf file missing. star index generation without prior knowledge"
				#nostar=true
				commander::warn "gtf file missing. proceeding without diego"
				nodsj=true
			elif [[ $FUSIONS ]]; then
				commander::warn "gtf file missing. proceeding without arriba"
				noarr=true
			elif [[ $tfq1 || $tmap ]]; then
				commander::warn "gtf file missing. proceeding without gem"
				nogem=true
			else
				commander::warn "gtf file missing. proceeding without quantification"
				noquant=true
				nodsj=true
				noclust=true
				nogo=true
			fi
		fi
	# }
fi
[[ -s "$GTF.go" ]] || nogo=true

if [[ $COMPARISONS ]]; then
	for f in "${COMPARISONS[@]}"; do
		BASHBONE_ERROR="experiment summary file for pairwise comparisons does not exists or is compressed $f"
		readlink -e "$f" | file -f - | grep -qF ASCII
		BASHBONE_ERROR="experiment summary file for pairwise comparisons is not tab separated $f"
		perl -F'\t' -lanE 'if($#F>=3){$m{$#F}=1}; END{@m=(keys %m); if($#m!=0){exit 1}}' "$f"
	done
fi

checkfile(){
	declare -n _idx=$2 _arr=$3
	local f ifs
	if [[ $(readlink -e "$1" | file -f - | grep ASCII) && -e "$(readlink -e $(head -1 "$1"))" ]]; then
		ifs=$IFS
		unset IFS
		while read -r f; do
			readlink -e "$f" &> /dev/null || return 1
			_arr[((++i))]="$(cd -P "$(dirname "$f")"; echo "$PWD/$(basename "$f")")"
			_idx+=($i)
		done < $1
		IFS=$ifs
	else
		f=$1
		readlink -e "$f" &> /dev/null || return 1
		_arr[((++i))]="$(cd -P "$(dirname "$f")"; echo "$PWD/$(basename "$f")")"
		_idx+=($i)
	fi
	return 0
}

IFS=','
i=-1
if [[ $tmap || $nmap ]]; then
	for f in $nmap; do
		BASHBONE_ERROR="alignment file does not exists $f"
		checkfile "$f" nidx MAPPED
	done
	for f in $nrmap; do
		BASHBONE_ERROR="normal replicate alignment file does not exists $f"
		checkfile "$f" nridx MAPPED
	done
	BASHBONE_ERROR="unequal number of normal replicates"
	[[ $nridx ]] && { [[ ${#nidx[@]} -eq ${#nridx[@]} ]] || false; }
	for f in $tmap; do
		BASHBONE_ERROR="treatment alignment file does not exists $f"
		checkfile "$f" pidx MAPPED
	done
	for f in $rmap; do
		BASHBONE_ERROR="treatment replicate alignment file does not exists $f"
		checkfile "$f" ridx MAPPED
	done
else
	for f in $nfq1; do
		BASHBONE_ERROR="single or first mate fastq file does not exist $f"
		checkfile "$f" nidx FASTQ1
	done
	for f in $nrfq1; do
		BASHBONE_ERROR="single or first mate normal replicate fastq file does not exists $f"
		checkfile "$f" nridx FASTQ1
	done
	for f in $tfq1; do
		BASHBONE_ERROR="single or first mate treatment fastq file does not exists $f"
		#idx depends on available replicates - i.e. trigger pooling or make pseudo-replicates
		checkfile "$f" $([[ $rfq1 ]] && echo tidx || echo pidx) FASTQ1
	done
	for f in $rfq1; do
		BASHBONE_ERROR="single or first mate treatment replicate fastq file does not exists $f"
		checkfile "$f" ridx FASTQ1
	done
	i=-1
	for f in {$nfq2,$nrfq2,$tfq2,$rfq2}; do
		BASHBONE_ERROR="second mate fastq file does not exists $f"
		checkfile "$f" foo FASTQ2
	done
	BASHBONE_ERROR="unequal number of mate pairs"
	[[ $FASTQ2 ]] && { [[ ${#FASTQ1[@]} -eq ${#FASTQ2[@]} ]] || false; }
	i=-1
	for f in {$nfq3,$nrfq3,$tfq3,$rfq3}; do
		BASHBONE_ERROR="umi fastq file does not exists $f"
		checkfile "$f" foo FASTQ3
	done
	BASHBONE_ERROR="unequal number of read and umi fastq files"
	[[ $FASTQ3 ]] && { [[ ${#FASTQ1[@]} -eq ${#FASTQ3[@]} ]] || false; }
fi
unset IFS

if $noidr; then
	tidx=("${pidx[@]}")
fi

if [[ $ridx ]]; then
	tidx=("${pidx[@]}")
	pidx=()
	BASHBONE_ERROR="unequal number of treatment replicates"
	[[ ${#ridx[@]} -eq ${#tidx[@]} ]] || false
fi

if [[ $FASTQ1 ]] && ${RRBS:=false}; then
	BASHBONE_ERROR="rrbs data analysis requires adapter sequence input"
	[[ ! $ADAPTER1 ]] && false
fi

if [[ $FASTQ3 ]]; then
	normd=false
fi

commander::printinfo "rippchen $VERSION utilizing bashbone $BASHBONE_VERSION started with command: $CMD" | tee -i "$LOG"
commander::printinfo "temporary files go to: $HOSTNAME:$TMPDIR" | tee -ia "$LOG"
commander::printinfo "date: $(date)" | tee -ia "$LOG"
x=$(ulimit -Hn)
[[ $((x/100)) -lt $THREADS ]] && {
	commander::warn "detected a low user limit of open file descriptors (ulimit -Hn : $x) for too many threads ($THREADS)"
	commander::warn "in case of memory allocation errors, you may decrease the number of threads to $((x/100))." | tee -ia "$LOG"
	commander::warn "possible memory allocation errors are 'bash: fork: Cannot allocate memory', 'Failed to read from standard input', 'Failed to open -', 'Too many open files'" | tee -ia "$LOG"
}

if ${INDEX:=false}; then
	BASHBONE_ERROR="indexing failed"
	progress::log -v $VERBOSITY -o "$LOG" -f pipeline::index
else
	if [[ $tfq1 || $tmap ]]; then
		${RIPSEQ:=false} && nosplitreads=${nosplitreads:-false} || nosplitreads=${nosplitreads:-true}
		BASHBONE_ERROR="peak calling pipeline failed"
		progress::log -v $VERBOSITY -o "$LOG" -f pipeline::callpeak
	elif [[ $FUSIONS ]]; then
		BASHBONE_ERROR="fusion detection pipeline failed"
		progress::log -v $VERBOSITY -o "$LOG" -f pipeline::fusions
	elif ${BISULFITE:=false}; then
		BASHBONE_ERROR="methylation analysis pipeline failed"
		progress::log -v $VERBOSITY -o "$LOG" -f pipeline::bs
	else
		[[ ${#nidx[@]} -lt 2 && "$noclust" != "true" ]] && {
			commander::warn "too few samples. proceeding without clustering"
			noclust=true
		}
		BASHBONE_ERROR="expression analysis pipeline failed"
		progress::log -v $VERBOSITY -o "$LOG" -f pipeline::dea
	fi
fi
unset BASHBONE_ERROR

${Smd5:=false} || {
	commander::printinfo "finally updating genome and annotation md5 sums" | tee -ia "$LOG"
	thismd5genome=$(md5sum "$GENOME" | cut -d ' ' -f 1)
	[[ "$md5genome" != "$thismd5genome" ]] && sed -i "s/md5genome=.*/md5genome=$thismd5genome/" "$GENOME.md5.sh"
	[[ $GTF ]] && {
		thismd5gtf=$(md5sum "$GTF" | cut -d ' ' -f 1)
		[[ "$md5gtf" != "$thismd5gtf" ]] && sed -i "s/md5gtf=.*/md5gtf=$thismd5gtf/" "$GENOME.md5.sh"
	}
}

exit 0
