#! /usr/bin/env bash
# (c) Konstantin Riege

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	{	compile::bashbone -i "$insdir" -t $threads && \
		compile::rippchen -i "$insdir" -t $threads && \
		compile::conda -i "$insdir" -t $threads && \
		compile::java -i "$insdir" -t $threads && \
		compile::sortmerna -i "$insdir" -t $threads && \
		compile::segemehl -i "$insdir" -t $threads && \
		compile::dexseq -i "$insdir" -t $threads && \
		compile::wgcna -i "$insdir" -t $threads && \
		compile::dgca -i "$insdir" -t $threads && \
		compile::revigo -i "$insdir" -t $threads && \
		compile::gem -i "$insdir" -t $threads && \
		compile::idr -i "$insdir" -t $threads && \
		compile::knapsack -i "$insdir" -t $threads
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

	commander::print "installing rippchen"
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
