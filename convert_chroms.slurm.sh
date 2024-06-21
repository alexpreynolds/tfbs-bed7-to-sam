#!/bin/bash

chrom=${1}
in_fn=${2}
threshold=${3}
scratch_dir=${4}
out_fn_prefix=${5}
cwd=${6}

${cwd}/bed7_to_sam.py ${in_fn} ${out_fn_prefix} ${threshold}

module add samtools/1.14
out_fn_sam=${out_fn_prefix}.sam
out_fn_bam=${out_fn_prefix}.bam
out_fn_bai=${out_fn_prefix}.bam.bai

# sort by coordinate
samtools sort ${out_fn_sam} -o ${out_fn_bam}

# index
samtools index ${out_fn_bam} ${out_fn_bai}
