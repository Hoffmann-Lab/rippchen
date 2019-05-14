#! /usr/bin/env bash
$RIPPCHEN/bin/rippchen/rippchen -v -t 8 -1 $RIPPCHEN/bin/rippchen/data/test.SE.fq.gz -o test_se
$RIPPCHEN/bin/rippchen/rippchen -v -t 8 -1 $RIPPCHEN/bin/rippchen/data/test.R1.fastq.gz -2 $RIPPCHEN/bin/rippchen/data/test.R2.fastq.gz -o test_pe
