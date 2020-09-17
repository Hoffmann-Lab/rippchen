#! /usr/bin/env bash
# (c) Konstantin Riege
set -o pipefail
trap 'trap - ERR; trap - RETURN' RETURN

unset ERROR
trap 'e=$?; echo ":ERROR: ${ERROR:-"..an unexpected one"} (exit $e) @ $(basename "${BASH_SOURCE[0]}") (line: $LINENO) $BASH_COMMAND" >&2; exit $e' ERR
ERROR="$(basename "${BASH_SOURCE[0]}") script needs to be sourced"
[[ "${BASH_SOURCE[0]}" != "$0" ]]
trap 'e=$?; echo ":ERROR: ${ERROR:-"..an unexpected one"} (exit $e) @ $(basename "${BASH_SOURCE[0]}") (line: $LINENO) $BASH_COMMAND" >&2; return $e' ERR

ERROR="requires a bash shell"
[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]]
ERROR="unsupported operating system"
[[ $OSTYPE == "linux" ]]
ERROR="requieres bash >= v4.4"
[[ ${BASH_VERSINFO[0]} -gt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4) ]]

unset ERROR
INSDIR_PIPELINE="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
INSDIR_TOOLS="$(dirname "$INSDIR_PIPELINE")"
unset OPTIND activate
while getopts :i:c: arg; do
	case $arg in
		i) INSDIR_PIPELINE="$OPTARG";;
		c) activate="$OPTARG";;
		:) ERROR="argument missing"; false;;
	esac
done

f="$INSDIR_PIPELINE/bashbone/activate.sh"
ERROR="file not found $f"
[[ -s "$f" ]]
source "$f" -i "$INSDIR_PIPELINE" -c ${activate:-false}
BASHBONEVERSION=$version
INSDIR="$INSDIR_PIPELINE"

_IFS=$IFS
IFS=$'\n'
for f in "$INSDIR/lib/"*.sh; do
	ERROR="file not found $f"
	[[ -s "$f" ]]
	source "$f"
done
IFS=$_IFS

trap 'trap - ERR; trap - RETURN' RETURN
trap 'configure::err -x $? -s "${BASH_SOURCE[0]}" -l $LINENO -e "$ERROR" -c "$BASH_COMMAND"; return $?' ERR

ERROR="environment setup failed. use -c false to disable tools and conda activation"
configure::environment -i "$INSDIR_TOOLS" -b "$INSDIR" -c ${activate:-false}
[[ $activate ]] || {
	commander::printinfo {COMMANDER[0]}<<- EOF
		to activate conda environment do
		source $(basename "${BASH_SOURCE[0]}") -c true
	EOF
}
