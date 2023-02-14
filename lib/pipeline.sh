#! /usr/bin/env bash
# (c) Konstantin Riege

pipeline::index(){
	genome::mkdict \
		-F \
		-t $THREADS \
		-i "$GENOME" \
		-p "$TMPDIR"
	if [[ -s "$GTF.go" ]]; then
		genome::mkgodb \
			-S ${nogo:=false} \
			-t $THREADS \
			-p "$TMPDIR" \
			-g "$GTF.go"
	fi
	if [[ -s "$GTF" ]]; then
		genome::indexgtf \
			-t $THREADS \
			-i "$GTF" \
			-F
	fi

	unset NA1 NA2
	if ${BISULFITE:=false}; then
		bisulfite::segemehl \
			-S ${nosege:=false} \
			-s true \
			-t $THREADS \
			-g "$GENOME" \
			-x "$GENOME.segemehl.ctidx" \
			-y "$GENOME.segemehl.gaidx" \
			-o "$TMPDIR" \
			-p "$TMPDIR" \
			-F \
			-r NA \
			-1 NA
		unset NA1 NA2
		bisulfite::bwa \
			-S ${nobwa:=false} \
			-s true \
			-t $THREADS \
			-g "$GENOME" \
			-o "$TMPDIR" \
			-F \
			-r NA \
			-1 NA
	else
		alignment::segemehl \
			-S ${nosege:=false} \
			-s true \
			-t $THREADS \
			-g "$GENOME" \
			-x "$GENOME.segemehl.idx" \
			-o "$TMPDIR" \
			-F \
			-r NA1 \
			-1 NA2
		unset NA1 NA2
		alignment::star \
			-S ${nostar:=false} \
			-s true \
			-t $THREADS \
			-g "$GENOME" \
			-x "$GENOME.star.idx" \
			-f "$GTF" \
			-o "$TMPDIR" \
			-p "$TMPDIR" \
			-F \
			-r NA1 \
			-1 NA2
		unset NA1 NA2
		alignment::bwa \
			-S ${nobwa:=false} \
			-s true \
			-t $THREADS \
			-g "$GENOME" \
			-x "$GENOME.bwa.idx/bwa" \
			-o "$TMPDIR" \
			-F \
			-r NA1 \
			-1 NA2
		unset NA1 NA2 NA3
		expression::diego \
			-S ${nodsj:=false} \
			-s true \
			-t $THREADS \
			-r NA1 \
			-x NA2 \
			-g "$GTF" \
			-c NA3 \
			-e false \
			-i "$TMPDIR" \
			-p "$TMPDIR" \
			-o "$TMPDIR" \
			-F
	fi

	return 0
}

pipeline::_slice(){
	alignment::slice \
		-S ${SLICED:-$1} \
		-s ${Sslice:-$2} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-r mapper \
		-c slicesinfo \
		-p "$TMPDIR"
	! $1 && ! $2 && SLICED=true
	! $1 && ${Sslice:-false} && SLICED=true

	return 0
}

pipeline::_preprocess(){
	if [[ ! $MAPPED ]]; then
		declare -a qualdirs

		local params=""
		[[ $ADAPTER1 ]] || params="-a ADAPTER1 -A ADAPTER2"
		preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/raw" -1 FASTQ1 -2 FASTQ2
		preprocess::fastqc \
			-S ${noqual:=false} \
			-s ${Sfqual:=false} \
			-t "$THREADS" \
			-M $MAXMEMORY \
			-o "$OUTDIR/qualities/raw" \
			-p "$TMPDIR" \
			-1 FASTQ1 \
			-2 FASTQ2 \
			$params

		${RRBS:=false} && ! ${nomspi:=false} && {
			bisulfite::mspicut \
				-S ${nomspi:=false} \
				-s ${Smspi:=false} \
				-t $THREADS \
				-d $DIVERSITY \
				-c 0 \
				-o "$OUTDIR/mspicut" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/mspicut" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Sfqual:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/mspicut" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		${notrim:=false} || {
			preprocess::trimmomatic \
				-S ${notrim:=false} \
				-s ${Strim:=false} \
				-b ${RRBS:=false} \
				-t "$THREADS" \
				-M $MAXMEMORY \
				-o "$OUTDIR/trimmed" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/trimmed" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Sfqual:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/trimmed" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		if [[ $ADAPTER1 ]]; then
			${noclip:=false} || {
				preprocess::cutadapt \
					-S ${noclip:=false} \
					-s ${Sclip:=false} \
					-a ADAPTER1 \
					-A ADAPTER2 \
					-b ${RRBS:=false} \
					-t $THREADS \
					-o "$OUTDIR/adapterclipped" \
					-1 FASTQ1 \
					-2 FASTQ2
				preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/adapterclipped" -1 FASTQ1 -2 FASTQ2
				preprocess::fastqc \
					-S ${noqual:=false} \
					-s ${Sfqual:=false} \
					-t $THREADS \
					-M $MAXMEMORY \
					-o "$OUTDIR/qualities/adapterclipped" \
					-p "$TMPDIR" \
					-1 FASTQ1 \
					-2 FASTQ2
			}
		fi

		# do after adapter clipping to not trim Ns before -U 2 in rrbs mode
		${noclip:=false} || {
			preprocess::rmpolynt \
				-S ${noclip:=false} \
				-s ${Sclip:=false} \
				-t $THREADS \
				-o "$OUTDIR/polyntclipped" \
				-d $(${BISULFITE:=false} && echo false || echo true) \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/polyntclipped" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Sfqual:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/polyntclipped" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		preprocess::rcorrector \
			-S ${nocor:=false} \
			-s ${Scor:=false} \
			-t $THREADS \
			-o "$OUTDIR/corrected" \
			-p "$TMPDIR" \
			-1 FASTQ1 \
			-2 FASTQ2

		${norrm:=false} || {
			preprocess::sortmerna \
				-S ${norrm:=false} \
				-s ${Srrm:=false} \
				-t $THREADS \
				-o "$OUTDIR/rrnafiltered" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/rrnafiltered" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Sfqual:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/rrnafiltered" \
				-p "$TMPDIR" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		preprocess::qcstats \
			-S ${nostats:=false} \
			-s ${Sstats:=false} \
			-t $THREADS \
			-i qualdirs \
			-o "$OUTDIR/stats" \
			-p "$TMPDIR" \
			-1 FASTQ1 \
			-2 FASTQ2
	fi

	return 0
}

pipeline::_mapping(){
	if [[ ! $MAPPED ]]; then
		if ${BISULFITE:=false}; then
			bisulfite::segemehl \
				-S ${nosege:=false} \
				-s ${Ssege:=false} \
				-5 ${Smd5:=false} \
				-1 FASTQ1 \
				-2 FASTQ2 \
				-o "$OUTDIR/mapped" \
				-p "$TMPDIR" \
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-i ${INSERTSIZE:=200000} \
				-g "$GENOME" \
				-x "$GENOME.segemehl.ctidx" \
				-y "$GENOME.segemehl.gaidx" \
				-r mapper

			bisulfite::bwa \
				-S ${nobwa:=false} \
				-s ${Sbwa:=false} \
				-5 ${Smd5:=false} \
				-1 FASTQ1 \
				-2 FASTQ2 \
				-o "$OUTDIR/mapped" \
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-g "$GENOME" \
				-r mapper
		else
			alignment::segemehl \
				-S ${nosege:=false} \
				-s ${Ssege:=false} \
				-5 ${Smd5:=false} \
				-1 FASTQ1 \
				-2 FASTQ2 \
				-o "$OUTDIR/mapped" \
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-i ${INSERTSIZE:=200000} \
				-n ${nosplitreads:=false} \
				-g "$GENOME" \
				-x "$GENOME.segemehl.idx" \
				-r mapper

			alignment::star \
				-S ${nostar:=false} \
				-s ${Sstar:=false} \
				-5 ${Smd5:=false} \
				-1 FASTQ1 \
				-2 FASTQ2 \
				-o "$OUTDIR/mapped" \
				-p "$TMPDIR" \
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-i ${INSERTSIZE:=200000} \
				-n ${nosplitreads:=false} \
				-g "$GENOME" \
				-f "$GTF" \
				-x "$GENOME.star.idx" \
				-r mapper

			! ${nosplitreads:=false} || alignment::bwa \
				-S ${nobwa:=false} \
				-s ${Sbwa:=false} \
				-5 ${Smd5:=false} \
				-1 FASTQ1 \
				-2 FASTQ2 \
				-o "$OUTDIR/mapped" \
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-f true \
				-g "$GENOME" \
				-x "$GENOME.bwa.idx/bwa" \
				-r mapper
		fi
	else
		declare -g -a ${MAPNAME:=custom}
		declare -n _MAPNAME_rippchen=$MAPNAME
		_MAPNAME_rippchen=("${MAPPED[@]}")
		mapper+=($MAPNAME)
	fi

	[[ ${#mapper[@]} -eq 0 ]] && return 0

	alignment::add4stats -r mapper
	alignment::bamqc \
		-S ${noqual:=false} \
		-s ${Smqual:=false} \
		-t $THREADS \
		-r mapper

	alignment::postprocess \
		-S ${nouniq:=false} \
		-s ${Suniq:=false} \
		-j uniqify \
		-t $THREADS \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped" \
		-r mapper
	! ${nouniq:=false} && ${nosort:=false} && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::postprocess \
		-S ${nosort:=false} \
		-s ${Ssort:=false} \
		-j sort \
		-t $THREADS \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped" \
		-r mapper
	alignment::postprocess \
		-S ${nosort:=false} \
		-s ${Ssort:=false} \
		-j index \
		-t $THREADS \
		-p "$TMPDIR" \
		-o "$OUTDIR/mapped" \
		-r mapper
	(${nouniq:=false} && ! ${nosort:=false}) || (! ${nouniq:=false} && ! ${nosort:=false}) && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	return 0
}

pipeline::fusions(){
	declare -a mapper

	pipeline::_preprocess

	fusions::arriba \
			-S ${noarr:=false} \
			-s ${Sarr:=false} \
			-5 ${Smd5:=false} \
			-t $THREADS \
			-g "$GENOME" \
			-v "$FUSIONS" \
			-a "$GTF" \
			-o "$OUTDIR/fusions" \
			-p "$TMPDIR" \
			-f $FRAGMENTSIZE \
			-1 FASTQ1 \
			-2 FASTQ2

	fusions::starfusion \
			-S ${nosfus:=false} \
			-s ${Ssfus:=false} \
			-5 ${Smd5:=false} \
			-t $THREADS \
			-g "$GENOME" \
			-o "$OUTDIR/fusions" \
			-p "$TMPDIR" \
			-1 FASTQ1 \
			-2 FASTQ2

	return 0
}

pipeline::bs(){
	declare -a mapper
	declare -A slicesinfo strandness

	pipeline::_preprocess
	pipeline::_mapping
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	genome::mkdict \
		-S ${normd:=true} \
		-s ${Srmd:=false} \
		-5 ${Smd5:=false} \
		-i "$GENOME" \
		-p "$TMPDIR" \
		-t $THREADS

	pipeline::_slice ${normd:=false} ${Srmd:=false}
	${normd:=false} || {
		alignment::rmduplicates \
			-S ${normd:=false} \
			-s ${Srmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-3 FASTQ3 \
			-c slicesinfo \
			-x "$REGEX" \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${normd:=false} \
			-s ${Srmd:=false} \
			-j index \
			-t $THREADS \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	pipeline::_slice ${nocmo:=false} ${Scmo:=false}
	${nocmo:=false} || {
		alignment::clipmateoverlaps \
			-S ${nocmo:=false} \
			-s ${Scmo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${nocmo:=false} \
			-s ${Scmo:=false} \
			-j index \
			-t $THREADS \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::qcstats \
		-S ${nostats:=false} \
		-s ${Sstats:=false} \
		-r mapper \
		-t $THREADS \
		-o "$OUTDIR/stats"

	bisulfite::methyldackel \
		-S ${nomedl:=false} \
		-s ${Smedl:=false} \
		-t $THREADS \
		-g "$GENOME" \
		-r mapper \
		-o "$OUTDIR/mecall"

	bisulfite::haarz \
		-S ${nohaarz:=false} \
		-s ${Shaarz:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-r mapper \
		-o "$OUTDIR/mecall" \
		-p "$TMPDIR"

	if ! ${nomedl:=false} || ! ${nohaarz:=false} && [[ $COMPARISONS ]]; then
		local params=""
		${nomedl:=false} || params+=" -d methyldackel"
		${nohaarz:=false} || params+=" -d haarz"
		bisulfite::metilene \
			-S ${nodma:=false} \
			-s ${Sdma:=false} \
			-t $THREADS \
			-c COMPARISONS \
			-m ${MINDATA:=0.8} \
			-u ${MINDATACAP:=999999} \
			-r mapper \
			-i "$OUTDIR/mecall" \
			-o "$OUTDIR/metilene" \
			-p "$TMPDIR" \
			$params
	fi

	return 0
}

pipeline::dea(){
	declare -a mapper coexpressions
	declare -A slicesinfo strandness

	pipeline::_preprocess
	pipeline::_mapping
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	genome::mkdict \
		-S ${normd:=true} \
		-s ${Srmd:=false} \
		-5 ${Smd5:=false} \
		-i "$GENOME" \
		-p "$TMPDIR" \
		-t $THREADS

	pipeline::_slice ${normd:=true} ${Srmd:=false}
	${normd:=true} || {
		alignment::rmduplicates \
			-S ${normd:=true} \
			-s ${Srmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-3 FASTQ3 \
			-c slicesinfo \
			-x "$REGEX" \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${normd:=true} \
			-s ${Srmd:=false} \
			-j index \
			-t $THREADS \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	pipeline::_slice ${nocmo:=true} ${Scmo:=false}
	${nocmo:=true} || {
		alignment::clipmateoverlaps \
			-S ${nocmo:=true} \
			-s ${Scmo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${nocmo:=true} \
			-s ${Scmo:=false} \
			-j index \
			-t $THREADS \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::qcstats \
		-S ${nostats:=false} \
		-s ${Sstats:=false} \
		-r mapper \
		-t $THREADS \
		-o "$OUTDIR/stats"

	alignment::inferstrandness \
		-S $( (${noquant:=false} && ${nodsj:=false}) && echo true || echo false) \
		-s $([[ $STRANDNESS ]] && echo true || { (${nodsj:=false} || ${Sdsj:=false}) && (${noquant:=false} || ${Squant:=false}) && echo true || echo false; }) \
		-d "$STRANDNESS" \
		-t $THREADS \
		-r mapper \
		-x strandness \
		-g "$GTF" \
		-l ${QUANTIFYFLEVEL:=exon} \
		-p "$TMPDIR"

	quantify::featurecounts \
		-S ${noquant:=false} \
		-s ${Squant:=false} \
		-t $THREADS \
		-p "$TMPDIR" \
		-g "$GTF" \
		-l ${QUANTIFYFLEVEL:=exon} \
		-f ${QUANTIFYFEATURE:=gene} \
		-o "$OUTDIR/counted" \
		-r mapper \
		-x strandness

	quantify::tpm \
		-S ${noquant:=false} \
		-s ${Stpm:=false} \
		-t $THREADS \
		-g "$GTF" \
		-l ${QUANTIFYFLEVEL:=exon} \
		-f ${QUANTIFYFEATURE:=gene} \
		-i "$OUTDIR/counted" \
		-r mapper

	if [[ $COMPARISONS ]]; then
		expression::join \
			-S $(${nodea:=false} && echo false || echo true) \
			-s ${Sjoin:=false} \
			-t $THREADS \
			-r mapper \
			-c COMPARISONS \
			-f ${QUANTIFYFEATURE:=gene} \
			-i "$OUTDIR/counted" \
			-o "$OUTDIR/counted"

		expression::diego \
			-S ${nodsj:=false} \
			-s ${Sdsj:=false} \
			-5 ${Smd5:=false} \
			-t $THREADS \
			-r mapper \
			-x strandness \
			-g "$GTF" \
			-c COMPARISONS \
			-i "$OUTDIR/counted" \
			-p "$TMPDIR" \
			-o "$OUTDIR/diego"

	 	expression::deseq \
			-S ${nodea:=false} \
			-s ${Sdea:=false} \
			-t $THREADS \
			-r mapper \
			-g "$GTF" \
			-f ${QUANTIFYFEATURE:=gene} \
			-c COMPARISONS \
			-i "$OUTDIR/counted" \
			-o "$OUTDIR/deseq"

		expression::join_deseq \
			-S ${nodea:=false} \
			-s ${Sjoin:=false} \
			-t $THREADS \
			-r mapper \
			-c COMPARISONS \
			-g "$GTF" \
			-f ${QUANTIFYFEATURE:=gene} \
			-i "$OUTDIR/counted" \
			-j "$OUTDIR/deseq" \
			-o "$OUTDIR/counted" \
			-p "$TMPDIR"

		cluster::coexpression_deseq \
			-S ${noclust:=false} \
			-s ${Sclust:=false} \
			-f ${CLUSTERFILTER:=04} \
			-b ${CLUSTERBIOTYPE:="."} \
			-g "$GTF" \
			-f ${QUANTIFYFEATURE:=gene} \
			-t $THREADS \
			-M $MAXMEMORY \
			-r mapper \
			-c COMPARISONS \
			-l coexpressions \
			-i "$OUTDIR/counted" \
			-j "$OUTDIR/deseq" \
			-o "$OUTDIR/coexpressed" \
			-p "$TMPDIR"

		enrichment::go \
			-S ${nogo:=false} \
			-s ${Sgo:=false} \
			-t $THREADS \
			-r mapper \
			-c COMPARISONS \
			-l coexpressions \
			-g "$GTF.go" \
			-i "$OUTDIR/deseq" \
			-j "$OUTDIR/counted"
	else
		cluster::coexpression \
			-S ${noclust:=false} \
			-s ${Sclust:=false} \
			-f ${CLUSTERFILTER:=4} \
			-b ${CLUSTERBIOTYPE:="."} \
			-g "$GTF" \
			-f ${QUANTIFYFEATURE:=gene} \
			-t $THREADS \
			-M $MAXMEMORY \
			-r mapper \
			-l coexpressions \
			-i "$OUTDIR/counted" \
			-o "$OUTDIR/coexpressed" \
			-p "$TMPDIR"

		enrichment::go \
			-S ${nogo:=false} \
			-s ${Sgo:=false} \
			-t $THREADS \
			-r mapper \
			-l coexpressions \
			-g "$GTF.go" \
			-j "$OUTDIR/counted"
	fi

	return 0
}

pipeline::callpeak() {
	declare -a mapper
	declare -A slicesinfo strandness

	pipeline::_preprocess
	pipeline::_mapping
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	${noidr:=false} || {
		alignment::mkreplicates \
			-S ${norep:=false} \
			-s ${Srep:=false} \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-p "$TMPDIR" \
			-r mapper \
			-n nidx \
			-m nridx \
			-i tidx \
			-j ridx \
			-k pidx
	}

	genome::mkdict \
		-S ${normd:=false} \
		-s ${Srmd:=false} \
		-5 ${Smd5:=false} \
		-i "$GENOME" \
		-p "$TMPDIR" \
		-t $THREADS

	pipeline::_slice ${normd:=false} ${Srmd:=false}
	${normd:=false} || {
		alignment::rmduplicates \
			-S ${normd:=false} \
			-s ${Srmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-3 FASTQ3 \
			-c slicesinfo \
			-x "$REGEX" \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${normd:=false} \
			-s ${Srmd:=false} \
			-j index \
			-t $THREADS \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	pipeline::_slice ${nocmo:=true} ${Scmo:=false}
	${nocmo:=true} || {
		alignment::clipmateoverlaps \
			-S ${nocmo:=true} \
			-s ${Scmo:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-o "$OUTDIR/mapped"
		alignment::postprocess \
			-S ${nocmo:=true} \
			-s ${Scmo:=false} \
			-j index \
			-t $THREADS \
			-p "$TMPDIR" \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Smqual:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::qcstats \
		-S ${nostats:=false} \
		-s ${Sstats:=false} \
		-r mapper \
		-t $THREADS \
		-o "$OUTDIR/stats"

	if ${noidr:=false}; then
		peaks::macs \
			-S ${nomacs:=false} \
			-s ${Smacs:=false} \
			-q ${RIPSEQ:=false} \
			-f $FRAGMENTSIZE \
			-g "$GENOME" \
			-a nidx \
			-i tidx \
			-r mapper \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-p "$TMPDIR" \
			-o "$OUTDIR/peaks" \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		${RIPSEQ:=false} && {
			alignment::inferstrandness \
				-S ${nogem:=false} \
				-s $([[ $STRANDNESS ]] && echo true || echo ${Sgem:=false}) \
				-d "$STRANDNESS" \
				-t $THREADS \
				-r mapper \
				-x strandness \
				-g "$GTF" \
				-p "$TMPDIR"
		}

		peaks::gem \
			-S ${nogem:=false} \
			-s ${Sgem:=false} \
			-q ${RIPSEQ:=false} \
			-g "$GENOME" \
			-a nidx \
			-i tidx \
			-r mapper \
			-x strandness \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-p "$TMPDIR" \
			-o "$OUTDIR/peaks" \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::peakachu \
			-S ${nopeaka:=false} \
			-s ${Speaka:=false} \
			-f $FRAGMENTSIZE \
			-a nidx \
			-i tidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-z ${STRICTPEAKS:=false}

		${RIPSEQ:=false} && [[ $tidx ]] && {
			peaks::matk \
				-S ${nomatk:=true} \
				-s ${Smatk:=false} \
				-a nidx \
				-i tidx \
				-r mapper \
				-t $THREADS \
				-M $MAXMEMORY \
				-p "$TMPDIR" \
				-o "$OUTDIR/peaks"

			peaks::m6aviewer \
				-S ${nom6a:=true} \
				-s ${Sm6a:=false} \
				-f $FRAGMENTSIZE \
				-a nidx \
				-i tidx \
				-r mapper \
				-t $THREADS \
				-m $MEMORY \
				-M $MAXMEMORY \
				-p "$TMPDIR" \
				-o "$OUTDIR/peaks"
		}
	else
		peaks::macs_idr \
			-S ${nomacs:=false} \
			-s ${Smacs:=false} \
			-q ${RIPSEQ:=false} \
			-f $FRAGMENTSIZE \
			-g "$GENOME" \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-p "$TMPDIR" \
			-o "$OUTDIR/peaks" \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		${RIPSEQ:=false} && {
			alignment::inferstrandness \
				-S ${nogem:=false} \
				-s $([[ $STRANDNESS ]] && echo true || echo ${Sgem:=false}) \
				-d "$STRANDNESS" \
				-t $THREADS \
				-r mapper \
				-x strandness \
				-g "$GTF" \
				-p "$TMPDIR"
		}

		peaks::gem_idr \
			-S ${nogem:=false} \
			-s ${Sgem:=false} \
			-q ${RIPSEQ:=false} \
			-g "$GENOME" \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-x strandness \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-p "$TMPDIR" \
			-o "$OUTDIR/peaks" \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::peakachu_idr \
			-S ${nopeaka:=false} \
			-s ${Speaka:=false} \
			-f $FRAGMENTSIZE \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-z ${STRICTPEAKS:=false}

		${RIPSEQ:=false} && {
			peaks::matk_idr \
				-S ${nomatk:=true} \
				-s ${Smatk:=false} \
				-a nidx \
				-b nridx \
				-i tidx \
				-j ridx \
				-k pidx \
				-r mapper \
				-t $THREADS \
				-M $MAXMEMORY \
				-p "$TMPDIR" \
				-o "$OUTDIR/peaks"

			peaks::m6aviewer_idr \
				-S ${nom6a:=true} \
				-s ${Sm6a:=false} \
				-f $FRAGMENTSIZE \
				-a nidx \
				-b nridx \
				-i tidx \
				-j ridx \
				-k pidx \
				-r mapper \
				-t $THREADS \
				-m $MEMORY \
				-M $MAXMEMORY \
				-p "$TMPDIR" \
				-o "$OUTDIR/peaks"
		}
	fi

	return 0
}
