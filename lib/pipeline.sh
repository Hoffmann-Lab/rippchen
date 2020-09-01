#! /usr/bin/env bash
# (c) Konstantin Riege

pipeline::index(){
	{	unset NA1 NA2 && \
		alignment::segemehl \
			-S ${nosege:=false} \
			-s true \
			-t $THREADS \
			-g $GENOME \
			-x $GENOME.segemehl.idx \
			-o $TMPDIR \
			-r NA1 \
			-1 NA2 && \
		unset NA1 NA2 && \
		alignment::star \
			-S ${nostar:=false} \
			-s true \
			-t $THREADS \
			-g $GENOME \
			-x $GENOME-staridx \
			-o $TMPDIR \
			-r NA1 \
			-1 NA2 && \
		genome::mkdict \
			-t $THREADS \
			-i $GENOME \
			-p $TMPDIR && \
		unset NA1 NA2 && \
		expression::diego \
			-S ${nodsj:=false} \
			-s true \
			-t $THREADS \
			-r NA1 \
			-g $GTF \
			-c NA2 \
			-i $TMPDIR \
			-j $TMPDIR \
			-p $TMPDIR \
			-o $TMPDIR
	} || return 1

	return 0
}

pipeline::_slice(){
	alignment::slice \
		-S ${SLICED:-$1} \
		-s ${Sslice:-$2} \
		-t $THREADS \
		-m $MEMORY \
		-r mapper \
		-c slicesinfo \
		-p $TMPDIR || return 1
	! $1 && ! $2 && SLICED=true
	! $1 && ${Sslice:-false} && SLICED=true

	return 0
}

pipeline::_preprocess(){
	if [[ ! $MAPPED ]]; then
		declare -a qualdirs

		{	qualdirs+=("$OUTDIR/qualities/raw") && \
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Squal:=false} \
				-t $THREADS \
				-o $OUTDIR/qualities/raw \
				-p $TMPDIR \
				-1 FASTQ1 \
				-2 FASTQ2
		} || return 1

		${notrim:=false} || { 
			{	qualdirs+=("$OUTDIR/qualities/trimmed") && \
				preprocess::trimmomatic \
					-S ${notrim:=false} \
					-s ${Strim:=false} \
					-t $THREADS \
					-o $OUTDIR/trimmed \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${noqual:=false} \
					-s ${Squal:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/trimmed \
					-p $TMPDIR \
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
						-p $TMPDIR \
						-1 FASTQ1 \
						-2 FASTQ2
				} || return 1
			}
		fi

		{	preprocess::rcorrector \
				-S ${nocor:=false} \
				-s ${Scor:=false} \
				-t $THREADS \
				-o $OUTDIR/corrected \
				-p $TMPDIR \
				-1 FASTQ1 \
				-2 FASTQ2
		} || return 1

		${norrm:=false} || {
			{	qualdirs+=("$OUTDIR/qualities/rrnafiltered") && \
				preprocess::sortmerna \
					-S ${norrm:=false} \
					-s ${Srrm:=false} \
					-t $THREADS \
					-m $MEMORY \
					-o $OUTDIR/rrnafiltered \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \
				preprocess::fastqc \
					-S ${noqual:=false} \
					-s ${Squal:=false} \
					-t $THREADS \
					-o $OUTDIR/qualities/rrnafiltered \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2
			} || return 1
		}


		{	preprocess::qcstats \
					-S ${nostats:=false} \
					-s ${Sstats:=false} \
					-i qualdirs \
					-o $OUTDIR/stats \
					-p $TMPDIR \
					-1 FASTQ1 \
					-2 FASTQ2 && \

			alignment::segemehl \
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

			alignment::star \
				-S ${nostar:=false} \
				-s ${Sstar:=false} \
				-5 ${Smd5:=false} \
				-1 FASTQ1 \
				-2 FASTQ2 \
				-o $OUTDIR/mapped \
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-i ${INSERTSIZE:=200000} \
				-p ${nosplitreads:=false} \
				-g $GENOME \
				-f "$GTF" \
				-x $GENOME-staridx \
				-r mapper
		} || return 1
	else
		declare -g -a ${MAPNAME:=custom}
		declare -n _MAPNAME_rippchen=$MAPNAME
		_MAPNAME_rippchen=("${MAPPED[@]}")
		mapper+=($MAPNAME)
	fi

	alignment::add4stats -r mapper
	
	return 0
}

pipeline::dea() {
	declare -a mapper coexpressions
	declare -A slicesinfo

	pipeline::_preprocess || return 1
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	{	genome::mkdict \
			-S ${nodict:=false} \
			-s ${Sdict:=false} \
			-5 ${Smd5:=false} \
			-i $GENOME \
			-p $TMPDIR \
			-t $THREADS && \

		alignment::postprocess \
			-S ${nouniq:=false} \
			-s ${Suniq:=false} \
			-j uniqify \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		${nouniq:=false} || alignment::add4stats -r mapper && \

		alignment::postprocess \
			-S ${nosort:=false} \
			-s ${Ssort:=false} \
			-j sort \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \

		pipeline::_slice ${normd:=true} ${Srmd:=false}
		alignment::rmduplicates \
			-S ${normd:=true} \
			-s ${Srmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-r mapper \
			-c slicesinfo \
			-x "$REGEX" \
			-p $TMPDIR \
			-o $OUTDIR/mapped && \
		${normd:=true} || alignment::add4stats -r mapper && \

		pipeline::_slice ${nocmo:=true} ${Scmo:=false} && \
		alignment::clipmateoverlaps \
			-S ${nocmo:=true} \
			-s ${Scmo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-r mapper \
			-c slicesinfo \
			-o $OUTDIR/mapped && \
		${nocmo:=true} || alignment::add4stats -r mapper && \
		
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

	if [[ ! $noquant && $COMPARISONS ]]; then
		{	expression::diego \
				-S ${nodsj:=false} \
				-s ${Sdsj:=false} \
				-5 ${Smd5:=false} \
				-t $THREADS \
				-r mapper \
				-g $GTF \
				-c COMPARISONS \
				-i $OUTDIR/counted \
				-j $OUTDIR/mapped \
				-p $TMPDIR \
				-o $OUTDIR/diego && \

		 	expression::deseq \
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

			cluster::coexpression_deseq \
				-S ${noclust:=false} \
				-s ${Sclust:=false} \
				-f ${CLUSTERFILTER:=0} \
				-b ${CLUSTERBIOTYPE:=""} \
				-g $GTF \
				-t $THREADS \
				-m $MEMORY \
				-r mapper \
				-c COMPARISONS \
				-l coexpressions \
				-i $OUTDIR/counted \
				-j $OUTDIR/deseq \
				-o $OUTDIR/coexpressed \
				-p $TMPDIR
		} || return 1
	else
		{	cluster::coexpression \
				-S ${noclust:=false} \
				-s ${Sclust:=false} \
				-f ${CLUSTERFILTER:=0} \
				-b ${CLUSTERBIOTYPE:=""} \
				-g $GTF \
				-t $THREADS \
				-m $MEMORY \
				-r mapper \
				-l coexpressions \
				-i $OUTDIR/counted \
				-o $OUTDIR/coexpressed \
				-p $TMPDIR
		} || return 1
	fi

	{	enrichment::go \
		-S ${nogo:=false} \
		-s ${Sgo:=false} \
		-t $THREADS \
		-r mapper \
		-c COMPARISONS \
		-l coexpressions \
		-g $GTF.go \
		-i $OUTDIR/deseq
	} || return 1

	return 0
}

pipeline::callpeak() {
	declare -a mapper caller
	declare -A slicesinfo

	pipeline::_preprocess || return 1
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	{	genome::mkdict \
			-S ${nodict:=false} \
			-s ${Sdict:=false} \
			-5 ${Smd5:=false} \
			-i $GENOME \
			-p $TMPDIR \
			-t $THREADS && \

	
    	alignment::postprocess \
			-S ${nouniq:=false} \
			-s ${Suniq:=false} \
			-j uniqify \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \
		${nouniq:=false} || alignment::add4stats -r mapper && \

		alignment::postprocess \
			-S ${nosort:=false} \
			-s ${Ssort:=false} \
			-j sort \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/mapped \
			-r mapper && \

		alignment::mkreplicates \
			-S ${norep:=false} \
			-s ${Srep:=false} \
			-t $THREADS \
			-o $OUTDIR/mapped \
			-p $TMPDIR \
			-r mapper \
			-n nidx \
			-m nridx \
			-i tidx \
			-j ridx \
			-k pidx && \

		pipeline::_slice ${normd:=false} ${Srmd:=false} && \
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
		${normd:=false} || alignment::add4stats -r mapper && \

		pipeline::_slice ${nocmo:=true} ${Scmo:=false} && \
		alignment::clipmateoverlaps \
			-S ${nocmo:=true} \
			-s ${Scmo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-r mapper \
			-c slicesinfo \
			-o $OUTDIR/mapped && \
		${nocmo:=true} || alignment::add4stats -r mapper && \
		

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

		callpeak::macs \
			-S ${nomacs:=false} \
			-s ${Smacs:=false} \
			-q ${RIPSEQ:=false} \
			-f $FRAGMENTSIZE \
			-g $GENOME \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-c caller \
			-t $THREADS \
			-m $MEMORY \
			-p $TMPDIR \
			-o $OUTDIR/peaks && \

		callpeak::gem \
			-S ${nogem:=false} \
			-s ${Sgem:=false} \
			-q ${RIPSEQ:=false} \
			-f $FRAGMENTSIZE \
			-g $GENOME \
			-f $GTF \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-c caller \
			-t $THREADS \
			-p $TMPDIR \
			-o $OUTDIR/peaks
	} || return 1

	return 0
}
