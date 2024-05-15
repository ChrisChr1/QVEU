#!/bin/bash

# SampleOrg.sh
# Usage: sh SampleOrg.sh [ARGS]

dir = /vp4dms/analysis/scRNAseq_analysis/Cell_Ranger/AAC2YHVHV_12/outs/fastq_path/AAC*
for i in S*;
do;
echo $i;

for f in `ls /hpcdata/lvd_qve/Projects/vp4dms/analysis/scRNAseq_analysis/Cell_Ranger/AAC2YHVHV_12/outs/fastq_path/AAC*/.fastq.gz`;
do echo $f; 
done; #enddir
done; #endsamples
