#! /usr/bin/env bash
# (c) Konstantin Riege

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	{	compile::bashbone -i "$insdir" -t $threads && \
		compile::rippchen -i "$insdir" -t $threads && \
		compile::conda -i "$insdir" -t $threads && \
		compile::java -i "$insdir" -t $threads && \
		compile::trimmomatic -i "$insdir" -t $threads && \
		compile::sortmerna -i "$insdir" -t $threads && \
		compile::segemehl -i "$insdir" -t $threads && \
		compile::knapsack -i "$insdir" -t $threads && \
		compile::dexseq -i "$insdir" -t $threads && \
		compile::wgcna -i "$insdir" -t $threads && \
		compile::dgca -i "$insdir" -t $threads && \
		compile::revigo -i "$insdir" -t $threads && \
		compile::gem -i "$insdir" -t $threads && \
		compile::idr -i "$insdir" -t $threads
	} || return 1

	return 0
}

compile::rippchen() {
	local insdir threads src=$(dirname $(dirname $(readlink -e ${BASH_SOURCE[0]})))
	compile::_parse -r insdir -s threads "$@" || return 1

	local version bashboneversion
	source $src/bashbone/lib/version.sh
	bashboneversion=$version
	source $src/lib/version.sh
	shopt -s extglob

	commander::printinfo "installing rippchen"
	{	rm -rf $insdir/rippchen-$version && \
		mkdir -p $insdir/rippchen-$version && \
		cp -r $src/!(bashbone|setup*) $insdir/rippchen-$version && \
		mkdir -p $insdir/latest && \
		ln -sfn $insdir/rippchen-$version $insdir/latest/rippchen && \
		ln -sfn $insdir/bashbone-$bashboneversion $insdir/rippchen-$version/bashbone
	} || return 1

	return 0
}

compile::upgrade(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	{	compile::bashbone -i "$insdir" -t $threads && \
		compile::rippchen -i "$insdir" -t $threads
	} || return 1
	
	return 0
}

compile::conda() {
	local insdir threads url version
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing conda and tools"
	{	url='https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && \
		wget -q -O $insdir/miniconda.sh $url && \
		version=$(bash $insdir/miniconda.sh -h | grep -F Installs | cut -d ' ' -f 3) && \
		rm -rf $insdir/conda && \
		mkdir -p $insdir/conda && \
		bash $insdir/miniconda.sh -b -f -p $insdir/conda && \
		rm $insdir/miniconda.sh && \
		source $insdir/conda/bin/activate && \

		conda env remove -y -n py2 && \
		conda env remove -y -n py3 && \
		conda create -y -n py2 python=2 && \
		conda create -y -n py2r python=2 && \
		conda create -y -n py3 python=3 && \

		conda install -n py2 -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			perl perl-threaded perl-db-file perl-dbi perl-app-cpanminus perl-list-moreutils perl-try-tiny perl-set-intervaltree perl-uri \
			numpy scipy pysam cython matplotlib \
			datamash \
			fastqc rcorrector \
			star star-fusion bwa hisat2 macs2 \
			samtools picard bamutil bedtools ucsc-facount khmer \
		chmod 755 $insdir/conda/envs/py2/bin/run_rcorrector.pl && \
		conda list -n py2 -f "fastqc|rcorrector|star|star-fusion|bwa|hisat2|macs2|samtools|picard|bamutil|bedtools|ucsc-facount|khmer" | grep -v '^#' > $insdir/condatools.txt && \

		conda install -n py3 -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			numpy scipy pysam cython matplotlib \
			cutadapt rseqc htseq diego bedtools && \
		conda list -n py3 -f "cutadapt|rseqc|diego|bedtools|htseq" | grep -v '^#' >> $insdir/condatools.txt && \

		conda install -n py2r -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			r-devtools bioconductor-biocinstaller bioconductor-biocparallel \
			bioconductor-genomicfeatures bioconductor-genefilter \
			subread r-wgcna bioconductor-deseq2 bioconductor-dexseq bioconductor-gseabase bioconductor-clusterprofiler \
			r-dplyr r-ggplot2 r-gplots r-rcolorbrewer r-svglite r-pheatmap r-ggpubr r-treemap r-rngtools && \
		conda list -n py2r -f "subread|r-wgcna|bioconductor-deseq2|bioconductor-dexseq" | grep -v '^#' >> $insdir/condatools.txt && \

		conda clean -y -a
	} || return 1

	return 0
}
