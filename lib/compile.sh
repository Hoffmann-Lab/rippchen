#! /usr/bin/env bash
# (c) Konstantin Riege

compile::bashbone() {
	local insdir threads version bashboneversion src="$(dirname "$(readlink -e "$0")")"
	commander::printinfo "installing bashbone"
	compile::_parse -r insdir -s threads "$@"
	source "$src/lib/version.sh"
	rm -rf "$insdir/bashbone-$version"
	mkdir -p "$insdir/bashbone-$version"
	cp -r "$src"/* "$insdir/bashbone-$version"
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"

	bashboneversion=$version
	src="$(dirname "$src")"

	commander::printinfo "installing rippchen"
	compile::_parse -r insdir -s threads "$@"
	source "$src/lib/version.sh"
	rm -rf "$insdir/rippchen-$version"
	mkdir -p "$insdir/rippchen-$version"
	cp -r "$src"/!(bashbone|setup*) "$insdir/rippchen-$version"
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/rippchen-$version" "$insdir/latest/rippchen"
	ln -sfn "$insdir/bashbone-$bashboneversion" "$insdir/rippchen-$version/bashbone"

	return 0
}

compile::rippchen() {
	compile::bashbone
}
