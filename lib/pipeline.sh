#! /usr/bin/env bash
# (c) Konstantin Riege

pipeline::_preprocess(){
	source $INSDIR/conda/bin/activate py2

	${Smd5:=false} || {
		[[ ! -s $GENOME.md5.sh ]] && cp $INSDIR/latest/bashbone/lib/md5.sh $GENOME.md5.sh
		source $GENOME.md5.sh
		thismd5genome=$(md5sum $GENOME | cut -d ' ' -f 1)
		[[ "$md5genome" != "$thismd5genome" ]] && sed -i "s/md5genome=.*/md5genome=$thismd5genome/" $GENOME.md5.sh
	}
	
	if [[ ! $MAPPED ]]; then
		declare -a qualdirs

		{	qualdirs+=("$OUTDIR/qualities/raw") && \
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Squal:=false} \
				-t $THREADS \
				-o $OUTDIR/qualities/raw \
				-1 FASTQ1 \
				-2 FASTQ2
		} || return 1
		${notrim:=false} || { 
			{	qualdirs+=("$OUTDIR/qualities/trimmed") && \
				preprocess::trimmomatic \
					-S ${notrim:=false} \
					-s ${Strim:=false} \
					-t $THREADS \
					-m $MEMORY \
					-o $OUTDIR/trimmed \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${noqual:=false} \
					-s ${Squal:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/trimmed \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		if [[ $ADAPTER1 ]]; then
			${noclip:=false} || {
				{	qualdirs+=("$OUTDIR/qualities/clipped") && \
					preprocess::cutadapt \
						-S ${noclip:=false} \
						-s ${Sclip:=false} \
						-a ADAPTER1 \
						-A ADAPTER2 \
						-t $THREADS \
						-o $OUTDIR/clipped \
						-1 FASTQ1 \
						-2 FASTQ2 && \
					preprocess::fastqc \
						-S ${noqual:=false} \
						-s ${Squal:=false} \
						-t $THREADS \
						-o $OUTDIR/qualities/clipped \
						-1 FASTQ1 \
						-2 FASTQ2
				} || return 1
			}
		fi
		${nocor:=false} || {
			{	# qualdirs+=("$OUTDIR/qualities/corrected") && \ no difference to trimmed in qcstats 
				preprocess::rcorrector \
					-S ${nocor:=false} \
					-s ${Scor:=false} \
					-t $THREADS \
					-o $OUTDIR/corrected \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${noqual:=false} \
					-s ${Squal:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/corrected \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		${norrm:=false} || {
			{	qualdirs+=("$OUTDIR/qualities/rrnafiltered") && \
				preprocess::sortmerna \
					-S ${norrm:=false} \
					-s ${Srrm:=false} \
					-t $THREADS \
					-m $MEMORY \
					-i $INSDIR \
					-o $OUTDIR/rrnafiltered \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${noqual:=false} \
					-s ${Squal:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/rrnafiltered \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		${nostats:=false} || {
			{	preprocess::qcstats \
					-S ${nostats:=false} \
					-s ${Sstats:=false} \
					-i qualdirs \
					-o $OUTDIR/stats \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}
		${nosege:=false} || {
			{	alignment::segemehl \
					-S ${nosege:=false} \
					-s ${Ssege:=false} \
					-5 ${Smd5:=false} \
					-1 FASTQ1 \
					-2 FASTQ2 \
					-o $OUTDIR/mapped \
					-t $THREADS \
					-a $((100-DISTANCE)) \
					-i ${INSERTSIZE:=200000} \
					-p ${nosplitreads:=false} \
					-g $GENOME \
					-x $GENOME.segemehl.idx \
					-r mapper && \
				alignment::add4stats -r mapper
			} || return 1
		}
	else
		custom=("${MAPPED[@]}")
		mapper+=(custom)
	fi

	return 0
}

pipeline::dea() {
	declare -a mapper coexpressions

	pipeline::_preprocess || return 1
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	{	alignment::postprocess \
			-S ${nouniq:=false} \
			-s ${Suniq:=false} \
			-j uniqify \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		alignment::add4stats -r mapper && \
		alignment::postprocess \
			-S ${nosort:=false} \
			-s ${Ssort:=false} \
			-j sort \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		alignment::postprocess \
			-S ${noidx:=false} \
			-s ${Sidx:=false} \
			-j index \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		alignment::bamstats \
			-S ${nostats:=false} \
			-s ${Sstats:=false} \
			-r mapper \
			-t $THREADS \
			-o $OUTDIR/stats && \
		quantify::featurecounts \
			-S ${noquant:=false} \
			-s ${Squant:=false} \
			-t $THREADS \
			-p $TMPDIR \
			-g $GTF \
			-l ${QUANTIFYFLEVEL:=exon} \
			-f ${QUANTIFYTAG:=gene_id} \
			-o $OUTDIR/counted \
			-r mapper && \
		quantify::tpm \
			-S ${noquant:=false} \
			-s ${Stpm:=false} \
			-t $THREADS \
			-g $GTF \
			-i $OUTDIR/counted \
			-r mapper
	} || return 1

	! $noquant && [[ $COMPARISONS ]] && {
		 {	expression::deseq \
				-S ${nodea:=false} \
				-s ${Sdea:=false} \
				-t $THREADS \
				-r mapper \
				-g $GTF \
				-c COMPARISONS \
				-i $OUTDIR/counted \
				-o $OUTDIR/deseq && \
			expression::joincounts \
				-S ${noquant:=false} \
				-s ${Sjoin:=false} \
				-t $THREADS \
				-r mapper \
				-c COMPARISONS \
				-i $OUTDIR/counted \
				-j $OUTDIR/deseq \
				-o $OUTDIR/counted \
				-p $TMPDIR && \
			cluster::coexpression \
				-S ${noclust:=false} \
				-s ${Sclust:=false} \
				-t $THREADS \
				-m $MEMORY \
				-r mapper \
				-c COMPARISONS \
				-z coexpressions \
				-i $OUTDIR/counted \
				-j $OUTDIR/deseq \
				-o $OUTDIR/coexpressed \
				-p $TMPDIR && \
			enrichment::go \
				-S ${nogo:=false} \
				-s ${Sgo:=false} \
				-t $THREADS \
				-r mapper \
				-c COMPARISONS \
				-l coexpressions \
				-g $GTF.go \
				-i $OUTDIR/deseq
		} || return 1
	}

	return 0
}

pipeline::_slice(){
	alignment::slice \
		-S $sliced \
		-s $(${Sslice:=false} && echo true || echo $1) \
		-t $THREADS \
		-m $MEMORY \
		-r mapper \
		-c slicesinfo \
		-p $TMPDIR || return 1

	$1 || sliced=true # i.e. if not skiptool: sliced=true and -S NOslice=true, else just by SKIPslices slicesinfo will be further updated

	return 0
}

pipeline::callpeak() {
	declare -a mapper
	declare -A slicesinfo

	pipeline::_preprocess || return 1
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	genome::mkdict \
		-S ${nodict:=false} \
		-s ${Sdict:=false} \
		-5 ${Smd5:=false} \
		-i $GENOME \
		-p $TMPDIR \
		-t $THREADS || return 1

	local sliced=false	
	
    {	alignment::postprocess \
			-S ${nouniq:=false} \
			-s ${Suniq:=false} \
			-j uniqify \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		alignment::add4stats -r mapper && \
		callpeak::mkreplicates && \
		alignment::postprocess \
			-S ${nosort:=false} \
			-s ${Ssort:=false} \
			-j sort \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		# alignment::postprocess \ <- applied by alignment::slice anyways
		# 	-S ${noidx:=false} \
		# 	-s ${Sidx:=false} \
		# 	-j index \
		# 	-t $THREADS \
		# 	-p $TMPDIR \
		# 	-o $OUTDIR/mapped \
		# 	-r mapper
		pipeline::_slice $($sliced || ${Srmd:=false} || ${normd:=false} && echo true || echo false) && \
		alignment::rmduplicates \
			-S ${normd:=false} \
			-s ${Srmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-r mapper \
			-c slicesinfo \
			-x "$REGEX" \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		alignment::add4stats -r mapper && \
		alignment::postprocess \
			-S ${noidx:=false} \
			-s ${Sidx:=false} \
			-j index \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		alignment::bamstats \
			-S ${nostats:=false} \
			-s ${Sstats:=false} \
			-r mapper \
			-t $THREADS \
			-o $OUTDIR/stats && \
		callpeak::macs_$IPTYPE && \
		callpeak::gem_$IPTYPE
	} || return 1

	return 0
}
