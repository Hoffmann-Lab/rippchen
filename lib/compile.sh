#! /usr/bin/env bash
# (c) Konstantin Riege

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::tools -i "$insdir" -t $threads
	compile::rippchen -i "$insdir" -t $threads
	compile::conda -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads
	compile::java -i "$insdir" -t $threads
	compile::trimmomatic -i "$insdir" -t $threads
	compile::sortmerna -i "$insdir" -t $threads
	compile::segemehl -i "$insdir" -t $threads
	compile::starfusion -i "$insdir" -t $threads
	compile::preparedexseq -i "$insdir" -t $threads
	compile::revigo -i "$insdir" -t $threads
	compile::gem -i "$insdir" -t $threads
	compile::m6aviewer -i "$insdir" -t $threads
	compile::idr -i "$insdir" -t $threads

	return 0
}

compile::rippchen() {
	local insdir threads version bashboneversion src=$(dirname $(dirname $(readlink -e ${BASH_SOURCE[0]})))
	commander::printinfo "installing rippchen"
	compile::_parse -r insdir -s threads "$@"
	#source $src/bashbone/lib/version.sh
	source $insdir/latest/bashbone/lib/version.sh
	bashboneversion=$version
	source $src/lib/version.sh
	rm -rf $insdir/rippchen-$version
	mkdir -p $insdir/rippchen-$version
	cp -r $src/!(bashbone|setup*) $insdir/rippchen-$version
	mkdir -p $insdir/latest
	ln -sfn $insdir/rippchen-$version $insdir/latest/rippchen
	ln -sfn $insdir/bashbone-$bashboneversion $insdir/rippchen-$version/bashbone

	return 0
}

compile::upgrade(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::tools -i "$insdir" -t $threads
	compile::rippchen -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads -u true

	return 0
}

compile::conda_tools() {
	local insdir threads upgrade=false url version tool n bin doclean=false
	declare -A envs

	compile::_parse -r insdir -s threads -c upgrade "$@"
	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	while read -r tool; do
		envs[$tool]=true
	done < <(conda info -e | awk -v prefix="^"$insdir '$NF ~ prefix {print $1}')

	for tool in fastqc cutadapt rcorrector star bwa rseqc subread htseq arriba picard bamutil macs2 peakachu diego metilene; do
		n=${tool/=*/}
		n=${n//[^[:alpha:]]/}
		$upgrade && ${envs[$n]:=false} || {
			doclean=true

			commander::printinfo "setup conda $n env"
			conda create -y -n $n #python=3
			conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		}
		# link commonly used base binaries into env
		for bin in perl bgzip samtools bedtools; do
			conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	done
	chmod 755 "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl" # necessary fix

	# manual setup of requirements from bioconda meta.yaml (see compile::starfusion) due to non-latest installation via conda
	# note: recent star is not compatible with CTAT plug-n-play genome index as of CTAT for star-fusion v1.9 (star indexer 2.7.1a)
	tool=star-fusion
	n=${tool/=*/}
	n=${n//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $n env"
		conda create -y -n $n #python=3
		# propably enought: perl perl-set-intervaltree perl-carp perl-carp-assert perl-db-file perl-io-gzip perl-json-xs perl-uri \
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			perl perl-file-path perl-getopt-long perl-set-intervaltree perl-carp perl-carp-assert perl-data-dumper perl-findbin perl-db-file perl-io-gzip perl-json-xs perl-uri perl-list-moreutils perl-list-util perl-storable \
			igv-reports star gmap bowtie bbmap samtools blast
	}
	for bin in perl bgzip samtools bedtools; do
		conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	done

	$doclean && {
		commander::printinfo "conda clean up"
		conda clean -y -a
	}

	conda deactivate
	return 0
}
