#! /usr/bin/env bash
# (c) Konstantin Riege

function compile::bashbone(){
	local insdir threads cfg version bashboneversion src="$(dirname "$(dirname "$(readlink -e "$0")")")"
	commander::printinfo "installing bashbone"
	compile::_parse -r insdir -s threads -f cfg "$@"
	source "$src/lib/version.sh"
	rm -rf "$insdir/bashbone-$version"
	mkdir -p "$insdir/bashbone-$version"
	cp -r "$src"/* "$insdir/bashbone-$version"
	rm -f "$insdir/bashbone-$version/scripts/"+(setup|test).sh
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"

	bashboneversion=$version
	src="$(dirname "$src")"

	commander::printinfo "installing rippchen"
	source "$src/lib/version.sh"
	rm -rf "$insdir/rippchen-$version"
	mkdir -p "$insdir/rippchen-$version"
	cp -r "$src"/!(bashbone) "$insdir/rippchen-$version"
	rm -f "$insdir/rippchen-$version/scripts/"+(setup|test).sh
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/rippchen-$version" "$insdir/latest/rippchen"
	ln -sfn "$insdir/bashbone-$bashboneversion" "$insdir/rippchen-$version/bashbone"

	return 0
}
