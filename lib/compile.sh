#! /usr/bin/env bash
# (c) Konstantin Riege

compile::all(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::rippchen -i "$insdir" -t $threads
	compile::conda -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads
	compile::java -i "$insdir" -t $threads
	compile::trimmomatic -i "$insdir" -t $threads
	compile::sortmerna -i "$insdir" -t $threads
	compile::segemehl -i "$insdir" -t $threads
	compile::preparedexseq -i "$insdir" -t $threads
	compile::revigo -i "$insdir" -t $threads
	compile::gem -i "$insdir" -t $threads
	compile::idr -i "$insdir" -t $threads

	return 0
}

compile::rippchen() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads version bashboneversion src=$(dirname $(dirname $(readlink -e ${BASH_SOURCE[0]})))
	commander::printinfo "installing rippchen"
	compile::_parse -r insdir -s threads "$@"
	source $src/bashbone/lib/version.sh
	bashboneversion=$version
	source $src/lib/version.sh
	shopt -s extglob
	rm -rf $insdir/rippchen-$version
	mkdir -p $insdir/rippchen-$version
	cp -r $src/!(bashbone|setup*) $insdir/rippchen-$version
	mkdir -p $insdir/latest
	ln -sfn $insdir/rippchen-$version $insdir/latest/rippchen
	ln -sfn $insdir/bashbone-$bashboneversion $insdir/rippchen-$version/bashbone

	return 0
}

compile::upgrade(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::rippchen -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads -u true

	return 0
}

compile::conda_tools() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads upgrade=false url version tool n bin doclean=false
	declare -A envs

	compile::_parse -r insdir -s threads -c upgrade "$@"
	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	while read -r tool; do
		envs[$tool]=true
	done < <(conda info -e | awk -v prefix="^"$insdir '$NF ~ prefix {print $1}')

	# python 3 envs
	for tool in fastqc cutadapt rcorrector star bwa rseqc subread arriba star-fusion picard bamutil macs2 diego; do
		n=${tool//[^[:alpha:]]/}
		$upgrade && ${envs[$n]:=false} && continue
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		# link commonly used base binaries into env
		for bin in perl samtools bedtools; do
			[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	done
	chmod 755 "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl" # necessary fix

	$doclean && {
		commander::printinfo "conda clean up"
		conda clean -y -a
	}

	conda deactivate
	return 0
}
