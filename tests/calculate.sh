#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load samtools 2>/dev/null >/dev/null
for each in $(find . -name '*.bam' -printf "%p "); do
	samtools flagstat ${each} | tr '\n' '\t'
	echo
	printf "$(basename "${each}")="
	samtools view ${each} | md5sum
done | sort | uniq | tr '\t' '\n'

# get a listing of the output files
ls -1
