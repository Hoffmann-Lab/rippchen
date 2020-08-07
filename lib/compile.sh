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
	local insdir threads
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
		
		# macs2, tophat2/hisat2 and some R stuff needs python2 whereas cutadapt,idr,rseqc need python3 env
		# star-fusion needs perl-set-intervaltree perl-db-file perl-set-intervaltree perl-uri perl-io-gzip
		#   installation might be fixed manually via perl-app-cpanminus and execution of cpanm Set::IntervalTree URI ...
		conda install -n py2 -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			perl perl-threaded perl-db-file perl-dbi perl-app-cpanminus perl-list-moreutils perl-try-tiny perl-set-intervaltree perl-uri \
			numpy scipy pysam cython matplotlib \
			datamash \
			fastqc rcorrector \
			star star-fusion bwa hisat2 macs2 \
			samtools picard bamutil ucsc-facount khmer \
		chmod 755 $insdir/conda/envs/py2/bin/run_rcorrector.pl && \
		conda list -n py2 -f "fastqc|rcorrector|star|star-fusion|bwa|hisat2|macs2|samtols|picard|bamutil|ucsc-facount|khmer" | grep -v '^#' > $insdir/condatools.txt && \

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

compile::gem() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing gem"
	{	source $insdir/conda/bin/activate py2 && \
		url='https://groups.csail.mit.edu/cgs/gem/download/gem.v3.4.tar.gz' && \
		version=$(basename $url | sed -E 's/.+v([0-9]+.+).tar.gz/\1/') && \
		wget -q $url -O $insdir/gem.tar.gz && \
		tar -xzf $insdir/gem.tar.gz -C $insdir && \
		mv $insdir/gem $insdir/gem-$version
		rm $insdir/gem.tar.gz && \
		cd $insdir/gem-$version && \
		mkdir -p bin && \
		wget -q -O bin/Read_Distribution_default.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_default.txt && \
		wget -q -O bin/Read_Distribution_CLIP.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_CLIP.txt && \
		compile::_javawrapper bin/gem $(readlink -e gem.jar) $insdir/latest/java/java && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/gem
	} || return 1

	return 0
}

compile::idr() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing idr"
	{	source $insdir/conda/bin/activate py3 && \
		url='https://github.com/kundajelab/idr/archive/2.0.4.2.tar.gz' && \
		wget -q $url -O $insdir/idr.tar.gz && \
		tar -xzf $insdir/idr.tar.gz -C $insdir && \
		rm $insdir/idr.tar.gz && \
		cd $(ls -vd $insdir/idr*/ | tail -1) && \
		pip install numpy matplotlib && \
		python setup.py install && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/idr
	} || return 1

	return 0
}

### OLD STUFF

compile::m6aviewer() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing m6aviewer"
	{	source $insdir/conda/bin/activate py2 && \
		url='http://dna2.leeds.ac.uk/m6a/m6aViewer_1_6_1.jar' && \
		mkdir -p $insdir/m6aViewer/bin && \
		wget -q $url -O $insdir/m6aViewer/m6aViewer_1_6_1.jar && \
		cd $insdir/m6aViewer && \
		compile::_javawrapper $PWD/bin/m6aViewer $(readlink -e m6aViewer_1_6_1.jar) $insdir/latest/java/java && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/m6aViewer
	} || return 1

	return 0
}

compile::metpeak() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing metpeak"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('compgenomics/MeTPeak', build_opts = c('--no-resave-data', '--no-manual'), threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::zerone() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing zerone"
	{	source $insdir/conda/bin/activate py2 && \
		cd $insdir && \
		rm -rf zerone && \
		git clone https://github.com/nanakiksc/zerone.git && \
		cd zerone && \
		make clean; true && \
		make -j $threads && \
		mkdir bin && \
		mv zerone bin && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/zerone
	} || return 1

	return 0
}

compile::dpgpc() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing dp_gp_cluster"
	{	source $insdir/conda/bin/activate py2 && \
		cd $insdir && \
		rm -rf DP_GP_cluster && \
		git clone https://github.com/PrincetonUniversity/DP_GP_cluster.git && \
		cd DP_GP_cluster && \
		sed -i -r '18,19{s/^#\s*//}' bin/DP_GP_cluster.py && \
		pip install GPy pandas numpy scipy matplotlib cython sklearn && \
		python setup.py install && \
		touch bin/DP_GP_cluster && \
		chmod 755 bin/* && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/DP_GP_cluster
	} || return 1

	cat <<- 'EOF' > $insdir/latest/DP_GP_cluster/DP_GP_cluster || return 1
		#!/usr/bin/env bash
		export PYTHONPATH=$CONDA_PREFIX/lib/python2.7/site-packages/:$PYTHONPATH
		$(cd $(dirname \$0) && echo $PWD)/DP_GP_cluster.py $*
	EOF
	return 0
}

compile::webgestalt() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing webgestalt"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WebGestaltR', threads=$threads, force=T)"
	} || return 1

	return 0
}
