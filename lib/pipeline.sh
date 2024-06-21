#! /usr/bin/env bash
# (c) Konstantin Riege

function pipeline::index(){
	genome::mkdict \
		-F \
		-t $THREADS \
		-i "$GENOME"
	if [[ -s "$GO" ]]; then
		genome::mkgodb \
			-S ${nogo:=false} \
			-t $THREADS \
			-g "$GO"
	fi
	if [[ -s "$GTF" ]]; then
		genome::indexgtf \
			-t $THREADS \
			-i "$GTF" \
			-F
	fi

	if ${BISULFITE:=false}; then
		unset NA1 NA2
		bisulfite::segemehl \
			-S ${nosege:=false} \
			-s true \
			-t $THREADS \
			-g "$GENOME" \
			-x "$GENOME.segemehl.ctidx" \
			-y "$GENOME.segemehl.gaidx" \
			-o "$TMPDIR" \
			-F \
			-r NA1 \
			-1 NA2
		unset NA1 NA2
		bisulfite::bwa \
			-S ${nobwa:=false} \
			-s true \
			-t $THREADS \
			-g "$GENOME" \
			-o "$TMPDIR" \
			-F \
			-r NA1 \
			-1 NA2
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
			-F \
			-r NA1 \
			-1 NA2
		unset NA1 NA2
		quantify::salmon \
			-S ${nosalm:=true} \
			-s true \
			-t $THREADS \
			-M $MAXMEMORY \
			-g "$GENOME" \
			-a "$GTF" \
			-i ${TRANSCRIPTOME:=false} \
			-x "$GENOME.salmon.idx" \
			-o "$TMPDIR" \
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
			-o "$TMPDIR" \
			-F
	fi

	return 0
}

function pipeline::_slice(){
	alignment::slice \
		-S ${SLICED:-$1} \
		-s ${Sslice:-$2} \
		-t $THREADS \
		-m $MEMORY \
		-M $MAXMEMORY \
		-r mapper \
		-c slicesinfo
	! $1 && ! $2 && SLICED=true
	! $1 && ${Sslice:-false} && SLICED=true

	return 0
}

function pipeline::_preprocess(){
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
				-s ${Smspi:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/mspicut" \
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
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/trimmed" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Strim:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/trimmed" \
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
					-s ${Sclip:=false} \
					-t $THREADS \
					-M $MAXMEMORY \
					-o "$OUTDIR/qualities/adapterclipped" \
					-1 FASTQ1 \
					-2 FASTQ2
			}
		fi

		# do after adapter clipping to not trim Ns before -U 2 in rrbs mode
		${nopclip:=false} || {
			preprocess::rmpolynt \
				-S ${nopclip:=false} \
				-s ${Spclip:=false} \
				-t $THREADS \
				-o "$OUTDIR/polyntclipped" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/polyntclipped" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Spclip:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/polyntclipped" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		preprocess::rcorrector \
			-S ${nocor:=false} \
			-s ${Scor:=false} \
			-t $THREADS \
			-o "$OUTDIR/corrected" \
			-1 FASTQ1 \
			-2 FASTQ2

		${norrm:=false} || {
			preprocess::sortmerna \
				-S ${norrm:=false} \
				-s ${Srrm:=false} \
				-t $THREADS \
				-o "$OUTDIR/rrnafiltered" \
				-1 FASTQ1 \
				-2 FASTQ2
			preprocess::add4stats -r qualdirs -a "$OUTDIR/qualities/rrnafiltered" -1 FASTQ1 -2 FASTQ2
			preprocess::fastqc \
				-S ${noqual:=false} \
				-s ${Srrm:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-o "$OUTDIR/qualities/rrnafiltered" \
				-1 FASTQ1 \
				-2 FASTQ2
		}

		preprocess::qcstats \
			-S ${nostats:=false} \
			-s ${Sstats:=false} \
			-t $THREADS \
			-i qualdirs \
			-o "$OUTDIR/stats" \
			-1 FASTQ1 \
			-2 FASTQ2
	fi

	return 0
}

function pipeline::_mapping(){
	if [[ ! $MAPPED ]]; then
		if ${BISULFITE:=false}; then
			bisulfite::segemehl \
				-S ${nosege:=false} \
				-s ${Ssege:=false} \
				-5 ${Smd5:=false} \
				-1 FASTQ1 \
				-2 FASTQ2 \
				-o "$OUTDIR/mapped" \
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-i ${INSERTSIZE:=1000} \
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
				-t $THREADS \
				-a $((100-DISTANCE)) \
				-i ${INSERTSIZE:=200000} \
				-n ${nosplitreads:=false} \
				-g "$GENOME" \
				-f "$GTF" \
				-x "$GENOME.star.idx" \
				-c $(${EMQUANT:-false} && ! ${TRANSCRIPTOME:-false} && echo true || echo false) \
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
		-M $MAXMEMORY \
		-o "$OUTDIR/mapped" \
		-r mapper
	! ${nouniq:=false} && ${nosort:=false} && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Suniq:=false} \
			-t $THREADS \
			-r mapper
	}

	alignment::postprocess \
		-S ${nosort:=false} \
		-s ${Ssort:=false} \
		-j sort \
		-t $THREADS \
		-M $MAXMEMORY \
		-o "$OUTDIR/mapped" \
		-r mapper
	(${nouniq:=false} && ! ${nosort:=false}) || (! ${nouniq:=false} && ! ${nosort:=false}) && {
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Ssort:=false} \
			-t $THREADS \
			-r mapper
	}

	${noblist:=true} || {
		alignment::postprocess \
			-S ${noblist:=false} \
			-s ${Sblist:=false} \
			-j blacklist \
			-f "$BLACKLIST" \
			-t $THREADS \
			-M $MAXMEMORY \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Sblist:=false} \
			-t $THREADS \
			-r mapper
	}

	${nofsel:=true} || {
		alignment::postprocess \
			-S ${nofsel:=false} \
			-s ${Sfsel:=false} \
			-j sizeselect \
			-f "$FRAGMENTSIZERANGE" \
			-t $THREADS \
			-M $MAXMEMORY \
			-o "$OUTDIR/mapped" \
			-r mapper
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Sfsel:=false} \
			-t $THREADS \
			-r mapper
	}

	return 0
}

function pipeline::fusions(){
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
			-f $FRAGMENTSIZE \
			-d "$STRANDNESS" \
			-1 FASTQ1 \
			-2 FASTQ2

	fusions::starfusion \
			-S ${nosfus:=false} \
			-s ${Ssfus:=false} \
			-5 ${Smd5:=false} \
			-t $THREADS \
			-g "$GENOME" \
			-o "$OUTDIR/fusions" \
			-1 FASTQ1 \
			-2 FASTQ2

	return 0
}

function pipeline::bs(){
	declare -a mapper
	declare -A slicesinfo strandness

	pipeline::_preprocess
	pipeline::_mapping
	[[ ${#mapper[@]} -eq 0 ]] && return 0

	pipeline::_slice ${normd:=false} ${Srmd:=false}
	${normd:=false} || {
		bisulfite::rmduplicates \
			-S ${normd:=false} \
			-s ${Srmd:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-g "$GENOME" \
			-r mapper \
			-3 FASTQ3 \
			-c slicesinfo \
			-o "$OUTDIR/mapped"
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Srmd:=false} \
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
		# alignment::add4stats -r mapper
		# alignment::bamqc \
		# 	-S ${noqual:=false} \
		# 	-s ${Scmo:=false} \
		# 	-t $THREADS \
		# 	-r mapper
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
		-x "$CONTEXT" \
		-r mapper \
		-o "$OUTDIR/mecall" \

	bisulfite::haarz \
		-S ${nohaarz:=false} \
		-s ${Shaarz:=false} \
		-t $THREADS \
		-M $MAXMEMORY \
		-g "$GENOME" \
		-x "$CONTEXT" \
		-r mapper \
		-o "$OUTDIR/mecall"

	if ! ${nomedl:=false} || ! ${nohaarz:=false} && [[ $COMPARISONS ]]; then
		local params=""
		${nomedl:=false} || params+=" -d methyldackel"
		${nohaarz:=false} || params+=" -d haarz"

		bisulfite::join \
			-S ${nojoin:=false} \
			-s ${Sjoin:=false} \
			-t $THREADS \
			-r mapper \
			-c COMPARISONS \
			-x "$CONTEXT" \
			-i "$OUTDIR/mecall" \
			-o "$OUTDIR/mecall" \
			$params

		bisulfite::metilene \
			-S ${nodma:=false} \
			-s ${Sdma:=false} \
			-t $THREADS \
			-c COMPARISONS \
			-x "$CONTEXT" \
			-m ${MINDATA:=0.8} \
			-u ${MINDATACAP:=999999} \
			-r mapper \
			-i "$OUTDIR/mecall" \
			-o "$OUTDIR/metilene" \
			$params
	fi

	return 0
}

function pipeline::dea(){
	declare -a mapper coexpressions
	declare -A slicesinfo strandness

	pipeline::_preprocess
	pipeline::_mapping
	if [[ ${#mapper[@]} -gt 0 ]]; then

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
				-o "$OUTDIR/mapped"
			alignment::add4stats -r mapper
			alignment::bamqc \
				-S ${noqual:=false} \
				-s ${Srmd:=false} \
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
			# alignment::add4stats -r mapper
			# alignment::bamqc \
			# 	-S ${noqual:=false} \
			# 	-s ${Scmo:=false} \
			# 	-t $THREADS \
			# 	-r mapper
		}

		alignment::qcstats \
			-S ${nostats:=false} \
			-s ${Sstats:=false} \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/stats"

		if ! ${EMQUANT:=false}; then
			#-S $( (${noquant:-false} && ${nodsj:-false}) && echo true || echo false) \
			#-s $([[ $STRANDNESS ]] && echo true || { (${nodsj:-false} || ${Sdsj:-false}) && (${noquant:-false} || ${Squant:-false}) && echo true || echo false; }) \
			alignment::inferstrandness \
				-S ${noquant:=false} \
				-s ${Squant:=false} \
				-d "$STRANDNESS" \
				-t $THREADS \
				-r mapper \
				-x strandness \
				-g "$GTF" \
				-l ${QUANTIFYFLEVEL:=exon}

			quantify::featurecounts \
				-S ${noquant:=false} \
				-s ${Squant:=false} \
				-t $THREADS \
				-M $MAXMEMORY \
				-g "$GTF" \
				-i ${TRANSCRIPTOME:=false} \
				-l ${QUANTIFYFLEVEL:=exon} \
				-f ${QUANTIFYFEATURE:=gene} \
				-o "$OUTDIR/counted" \
				-r mapper \
				-x strandness

			if ! ${nosplitreads:=false} && [[ $COMPARISONS ]]; then
				expression::diego \
					-S ${nodsj:=false} \
					-s ${Sdsj:=false} \
					-5 ${Smd5:=false} \
					-t $THREADS \
					-r mapper \
					-e $(${noquant:=false} && echo false || echo true) \
					-x strandness \
					-g "$GTF" \
					-c COMPARISONS \
					-i "$OUTDIR/counted" \
					-o "$OUTDIR/diego"
			fi
		fi
	fi

	# run after diego and misuse mapper array to store count file paths in
	if ${EMQUANT:=false}; then
		alignment::postprocess \
			-S ${nosaln:=true} \
			-s ${Ssaln:=false} \
			-j collate \
			-t $THREADS \
			-M $MAXMEMORY \
			-o "$OUTDIR/mapped" \
			-r mapper

		quantify::salmon \
			-S ${nosaln:=true} \
			-s ${Ssaln:=false} \
			-5 ${Smd5:=false} \
			-o "$OUTDIR/counted" \
			-t $THREADS \
			-M $MAXMEMORY \
			-f ${QUANTIFYFEATURE:=gene} \
			-g "$GENOME" \
			-a "$GTF" \
			-i ${TRANSCRIPTOME:=false} \
			-x "$GENOME.salmon.idx" \
			-d "$STRANDNESS" \
			-r mapper
	fi

	quantify::salmon \
		-S ${nosalm:=true} \
		-s ${Ssalm:=false} \
		-5 ${Smd5:=false} \
		-1 FASTQ1 \
		-2 FASTQ2 \
		-o "$OUTDIR/counted" \
		-t $THREADS \
		-M $MAXMEMORY \
		-f ${QUANTIFYFEATURE:=gene} \
		-g "$GENOME" \
		-a "$GTF" \
		-i ${TRANSCRIPTOME:=false} \
		-x "$GENOME.salmon.idx" \
		-r mapper

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
			-S $(${noquant:-false} && echo true || { ${nodea:-false} && echo false || echo true; }) \
			-s ${Sjoin:=false} \
			-t $THREADS \
			-r mapper \
			-c COMPARISONS \
			-f ${QUANTIFYFEATURE:=gene} \
			-i "$OUTDIR/counted" \
			-o "$OUTDIR/counted"

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
			-S $( ${noquant:-false} && echo true || { ${nodea:-false} && echo true || echo false; }) \
			-s ${Sjoin:=false} \
			-t $THREADS \
			-r mapper \
			-c COMPARISONS \
			-g "$GTF" \
			-f ${QUANTIFYFEATURE:=gene} \
			-i "$OUTDIR/counted" \
			-j "$OUTDIR/deseq" \
			-o "$OUTDIR/counted"

		cluster::coexpression_deseq \
			-S ${noclust:=false} \
			-s ${Sclust:=false} \
			-n ${CLUSTERFILTER:=04} \
			-b ${FEATUREBIOTYPE:="protein_coding"} \
			-g "$GTF" \
			-f ${QUANTIFYFEATURE:=gene} \
			-t $THREADS \
			-M $MAXMEMORY \
			-r mapper \
			-c COMPARISONS \
			-l coexpressions \
			-i "$OUTDIR/counted" \
			-j "$OUTDIR/deseq" \
			-o "$OUTDIR/coexpressed"

		enrichment::go \
			-S ${nogo:=false} \
			-s ${Sgo:=false} \
			-t $THREADS \
			-r mapper \
			-c COMPARISONS \
			-l coexpressions \
			-d "$GO" \
			-g "$GTF" \
			-b ${FEATUREBIOTYPE:="protein_coding"} \
			-f ${QUANTIFYFEATURE:=gene} \
			-i "$OUTDIR/deseq" \
			-j "$OUTDIR/counted"
	else
		cluster::coexpression \
			-S ${noclust:=false} \
			-s ${Sclust:=false} \
			-n ${CLUSTERFILTER:=4} \
			-b ${FEATUREBIOTYPE:="protein_coding"} \
			-g "$GTF" \
			-f ${QUANTIFYFEATURE:=gene} \
			-t $THREADS \
			-M $MAXMEMORY \
			-r mapper \
			-l coexpressions \
			-i "$OUTDIR/counted" \
			-o "$OUTDIR/coexpressed"

		enrichment::go \
			-S ${nogo:=false} \
			-s ${Sgo:=false} \
			-t $THREADS \
			-r mapper \
			-l coexpressions \
			-d "$GO" \
			-g "$GTF" \
			-b ${FEATUREBIOTYPE:="protein_coding"} \
			-f ${QUANTIFYFEATURE:=gene} \
			-j "$OUTDIR/counted"
	fi

	return 0
}

function pipeline::callpeak(){
	declare -a mapper
	declare -A slicesinfo strandness

	pipeline::_preprocess
	pipeline::_mapping
	[[ ${#mapper[@]} -eq 0 ]] && return 0

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
			-o "$OUTDIR/mapped"
		alignment::add4stats -r mapper
		alignment::bamqc \
			-S ${noqual:=false} \
			-s ${Srmd:=false} \
			-t $THREADS \
			-r mapper
	}

	pipeline::_slice ${noctn5:=true} ${Sctn5:=false}
	${noctn5:=true} || {
		alignment::clip \
			-S ${noctn5:=true} \
			-s ${Sctn5:=false} \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-r mapper \
			-c slicesinfo \
			-g $GENOME \
			-5 4 \
			-3 0 \
			-5 5 \
			-3 0 \
			-o "$OUTDIR/mapped"
		# alignment::add4stats -r mapper
		# alignment::bamqc \
		# 	-S ${noqual:=false} \
		# 	-s ${Sctn5:=false} \
		# 	-t $THREADS \
		# 	-r mapper
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
		alignment::add4stats -r mapper
		# alignment::bamqc \
		# 	-S ${noqual:=false} \
		# 	-s ${Scmo:=false} \
		# 	-t $THREADS \
		# 	-r mapper
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
			-o "$OUTDIR/peaks" \
			-w ${BROAD:=false} \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		${RIPSEQ:=false} && {
			alignment::inferstrandness \
				-S ${nogem:=false} \
				-s $([[ $STRANDNESS ]] && echo true || echo ${Sgem:-false}) \
				-d "$STRANDNESS" \
				-t $THREADS \
				-r mapper \
				-x strandness \
				-g "$GTF"
		}

		peaks::gem \
			-S ${nogem:=false} \
			-s ${Sgem:=false} \
			-q ${RIPSEQ:=false} \
			-g "$GENOME" \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-f $FRAGMENTSIZE \
			-a nidx \
			-i tidx \
			-r mapper \
			-x strandness \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-o "$OUTDIR/peaks" \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::homer \
			-S ${nohomer:=false} \
			-s ${Shomer:=false} \
			-q ${RIPSEQ:=false} \
			-g "$GENOME" \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-f $FRAGMENTSIZE \
			-a nidx \
			-i tidx \
			-r mapper \
			-x strandness \
			-t $THREADS \
			-m $MEMORY \
			-M $MAXMEMORY \
			-o "$OUTDIR/peaks" \
			-w ${BROAD:=false} \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::genrich \
			-S ${norich:=false} \
			-s ${Srich:=false} \
			-f $FRAGMENTSIZE \
			-a nidx \
			-i tidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-q ${RIPSEQ:=false} \
			-w ${BROAD:=false} \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::seacr \
			-S ${noseacr:=false} \
			-s ${Sseacr:=false} \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-m $MEMORY \
			-M $MAXMEMORY \
			-f $FRAGMENTSIZE \
			-g "$GENOME" \
			-a nidx \
			-i tidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-z ${STRICTPEAKS:=false}

		peaks::gopeaks \
			-S ${nogopeaks:=false} \
			-s ${Sgopeaks:=false} \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-m $MEMORY \
			-M $MAXMEMORY \
			-f $FRAGMENTSIZE \
			-g "$GENOME" \
			-a nidx \
			-i tidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-w ${BROAD:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::peakachu \
			-S ${nopeaka:=false} \
			-s ${Speaka:=false} \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-g "$GENOME" \
			-m $MEMORY \
			-M $MAXMEMORY \
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
				-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
				-m $MEMORY \
				-f $FRAGMENTSIZE \
				-g "$GENOME" \
				-a nidx \
				-i tidx \
				-r mapper \
				-t $THREADS \
				-M $MAXMEMORY \
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
				-o "$OUTDIR/peaks"
		}
	else
		alignment::mkreplicates \
			-S ${norep:=false} \
			-s ${Srep:=false} \
			-t $THREADS \
			-o "$OUTDIR/mapped" \
			-r mapper \
			-n nidx \
			-m nridx \
			-i tidx \
			-j ridx \
			-k pidx

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
			-o "$OUTDIR/peaks" \
			-w ${BROAD:=false} \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		${RIPSEQ:=false} && {
			alignment::inferstrandness \
				-S ${nogem:=false} \
				-s $([[ $STRANDNESS ]] && echo true || echo ${Sgem:-false}) \
				-d "$STRANDNESS" \
				-t $THREADS \
				-r mapper \
				-x strandness \
				-g "$GTF"
		}

		peaks::gem_idr \
			-S ${nogem:=false} \
			-s ${Sgem:=false} \
			-q ${RIPSEQ:=false} \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-f $FRAGMENTSIZE \
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
			-o "$OUTDIR/peaks" \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::homer_idr \
			-S ${nohomer:=false} \
			-s ${Shomer:=false} \
			-q ${RIPSEQ:=false} \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-f $FRAGMENTSIZE \
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
			-o "$OUTDIR/peaks" \
			-w ${BROAD:=false} \
			-y ${POINTYPEAKS:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::genrich_idr \
			-S ${norich:=false} \
			-s ${Srich:=false} \
			-f $FRAGMENTSIZE \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-w ${BROAD:=false} \
			-q ${RIPSEQ:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::seacr_idr \
			-S ${noseacr:=false} \
			-s ${Sseacr:=false} \
			-m $MEMORY \
			-M $MAXMEMORY \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-f $FRAGMENTSIZE \
			-g "$GENOME" \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-z ${STRICTPEAKS:=false}

		peaks::gopeaks_idr \
			-S ${nogopeaks:=false} \
			-s ${Sgopeaks:=false} \
			-m $MEMORY \
			-M $MAXMEMORY \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-f $FRAGMENTSIZE \
			-g "$GENOME" \
			-a nidx \
			-b nridx \
			-i tidx \
			-j ridx \
			-k pidx \
			-r mapper \
			-t $THREADS \
			-o "$OUTDIR/peaks" \
			-w ${BROAD:=false} \
			-z ${STRICTPEAKS:=false}

		peaks::peakachu_idr \
			-S ${nopeaka:=false} \
			-s ${Speaka:=false} \
			-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
			-g "$GENOME" \
			-m $MEMORY \
			-M $MAXMEMORY \
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
				-m $MEMORY \
				-M $MAXMEMORY \
				-c "$(${nomacs:-false} || echo "$OUTDIR/peaks")" \
				-f $FRAGMENTSIZE \
				-g "$GENOME" \
				-a nidx \
				-b nridx \
				-i tidx \
				-j ridx \
				-k pidx \
				-r mapper \
				-t $THREADS \
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
				-o "$OUTDIR/peaks"
		}
	fi

	return 0
}
