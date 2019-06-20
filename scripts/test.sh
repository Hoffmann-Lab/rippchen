#! /usr/bin/env bash

THREADS=$(grep -cF processor /proc/cpuinfo)
[[ $RIPPCHEN ]] || {
	echo -n "where is rippchen? [/path/to/install/dir]: "
	read -r p
	export RIPPCHEN=$(readlink -e ${p/#*~/$HOME})
}

### preprocessing
prep(){
	$RIPPCHEN/latest/rippchen/rippchen.sh \
		-v 2 \
		-t $THREADS \
		-1 $RIPPCHEN/latest/rippchen/data/test.SE.fq.gz \
		-o /ssd/tmp/rippchen_test_se

	$RIPPCHEN/latest/rippchen/rippchen.sh \
		-v 2 \
		-t $THREADS \
		-1 $RIPPCHEN/latest/rippchen/data/test.R1.fastq.gz \
		-2 $RIPPCHEN/latest/rippchen/data/test.R2.fastq.gz \
		-o /ssd/tmp/rippchen_test_pe
}
echo -n "test preprocessing? [y|n]: "
read -r yn
[[ $yn == "y" ]] && prep

### mapping
map(){
	$RIPPCHEN/latest/rippchen/rippchen.sh \
		-v 2 \
		-t $THREADS \
		-g $RIPPCHEN/latest/rippchen/data/test.fa \
		-1 $RIPPCHEN/latest/rippchen/data/test.SE.fq.gz \
		-o /ssd/tmp/rippchen_test_se \
		-resume sege

	$RIPPCHEN/latest/rippchen/rippchen.sh \
		-v 2 \
		-t $THREADS \
		-g $RIPPCHEN/latest/rippchen/data/test.fa \
		-1 $RIPPCHEN/latest/rippchen/data/test.R1.fastq.gz \
		-2 $RIPPCHEN/latest/rippchen/data/test.R2.fastq.gz \
		-o /ssd/tmp/rippchen_test_pe \
		-resume sege
}
echo -n "test mapping? [y|n]: "
read -r yn
[[ $yn == "y" ]] && map

quant(){
	$RIPPCHEN/latest/rippchen/rippchen.sh \
		-v 2 \
		-t $THREADS \
		-g $RIPPCHEN/latest/rippchen/data/test.fa \
		-gtf $RIPPCHEN/latest/rippchen/data/test.gtf \
		-1 $RIPPCHEN/latest/rippchen/data/test.SE.fq.gz \
		-o /ssd/tmp/rippchen_test_se \
		-resume quant

	$RIPPCHEN/latest/rippchen/rippchen.sh \
		-v 2 \
		-t $THREADS \
		-g $RIPPCHEN/latest/rippchen/data/test.fa \
		-gtf $RIPPCHEN/latest/rippchen/data/test.gtf \
		-1 $RIPPCHEN/latest/rippchen/data/test.R1.fastq.gz \
		-2 $RIPPCHEN/latest/rippchen/data/test.R2.fastq.gz \
		-o /ssd/tmp/rippchen_test_pe \
		-resume quant
}
echo -n "test quantification? [y|n]: "
read -r yn
[[ $yn == "y" ]] && quant

coex(){
	source $RIPPCHEN/conda/bin/activate py2
	INSDIR=$RIPPCHEN
	for f in {$INSDIR/latest/bashbone/lib/*.sh,$INSDIR/latest/rippchen/lib/*.sh}; do
		source $f
	done
	configure::environment -i $INSDIR

	mkdir -p /ssd/tmp/rippchen_test_counted/segemehl
	cp /misc/paras/data/kons/holger_tumorDriver/RNA-Seq/results/counted/segemehl/*counts* /ssd/tmp/rippchen_test_counted/segemehl
	chmod 644 /ssd/tmp/rippchen_test_counted/segemehl/*

	mapper=(segemehl)
	comparisons=(/misc/paras/data/kons/holger_tumorDriver/RNA-Seq/COMPARISONS)
	coexpressions=()
	expression::deseq \
		-S false \
		-s false \
		-t $THREADS \
		-r mapper \
		-c comparisons \
		-i /ssd/tmp/rippchen_test_counted \
		-o /ssd/tmp/rippchen_test_deseq

	expression::joincounts \
		-S false \
		-s false \
		-t $THREADS \
		-r mapper \
		-c comparisons \
		-i /ssd/tmp/rippchen_test_counted \
		-j /ssd/tmp/rippchen_test_deseq \
		-o /ssd/tmp/rippchen_test_counted \
		-p /ssd/tmp/rippchen_test_tmp

	cluster::coexpression \
		-S false \
		-s false \
		-t $THREADS \
		-m 30000 \
		-r mapper \
		-c comparisons \
		-z coexpressions \
		-i /ssd/tmp/rippchen_test_counted \
		-j /ssd/tmp/rippchen_test_deseq \
		-o /ssd/tmp/rippchen_test_coexpressed \
		-p /ssd/tmp/rippchen_test_tmp

	enrichment::go \
		-S false \
		-s false \
		-t $THREADS \
		-r mapper \
		-c comparisons \
		-l coexpressions \
		-i /ssd/tmp/rippchen_test_deseq \
		-g /misc/paras/data/genomes/GRCm38.p6/GRCm38.p6.fa.gtf.go

}
echo -n "test deseq and coexpression? [y|n]: "
read -r yn
[[ $yn == "y" ]] && coex

slice(){
	source $RIPPCHEN/conda/bin/activate py2
	INSDIR=$RIPPCHEN
	for f in {$INSDIR/latest/bashbone/lib/*.sh,$INSDIR/latest/rippchen/lib/*.sh}; do
		source $f
	done
	configure::environment -i $INSDIR

	declare -a segemehl=("/misc/paras/data/kons/holger_tumorDriver/RNA-Seq/results/mapped/segemehl/wt_1597.N1.unique.sorted.bam" "/misc/paras/data/kons/holger_tumorDriver/RNA-Seq/results/mapped/segemehl/tu_550.N3.unique.sorted.bam")
	declare -a mapper=(segemehl)

	declare -A slicesinfo
	alignment::slice \
		-S false \
		-s false \
		-t $THREADS \
		-m 30000 \
		-r mapper \
		-c slicesinfo \
		-p /ssd/tmp/rippchen_test_slices

	for f in "${segemehl[@]}"; do
		echo "$f"
		cat ${slicesinfo["$f"]}
	done

	alignment::rmduplicates \
		-S false \
		-s false \
		-t $THREADS \
		-m 30000 \
		-r mapper \
		-c slicesinfo \
		-p /ssd/tmp/rippchen_test_tmp \
		-o /ssd/tmp/rippchen_test_rmdup
}
echo -n "test slicing? [y|n]: "
read -r yn
[[ $yn == "y" ]] && slice
