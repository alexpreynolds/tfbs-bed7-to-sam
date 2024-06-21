#!/bin/bash

in_fn=${1}
threshold=${2}

SCRATCH_PREFIX=/net/seq/data/projects/fsv/scratch

slurm_partition=queue1
slurm_cpus_per_task=1
slurm_mem_per_cpu="10G"
slurm_ntasks=1

for chr in `seq 1 22` X Y
do
    chrom="chr${chr}"
    root_out_fn=${in_fn%.*}.${chrom}.bed
    out_fn=${PWD}/$(basename ${root_out_fn})
    echo "${chrom} | ${out_fn}"
    slurm_job_name="extract_tfbs_${chrom}_${threshold}"
    slurm_out=${PWD}/${slurm_job_name}.out
    slurm_err=${PWD}/${slurm_job_name}.err
    scratch_dir=${SCRATCH_PREFIX}/tfbs/${threshold}/${chrom}
    mkdir -p ${scratch_dir}
    sbatch \
        --partition=${slurm_partition} \
        --cpus-per-task=${slurm_cpus_per_task} \
        --ntasks=${slurm_ntasks} \
        --mem-per-cpu=${slurm_mem_per_cpu} \
        --job-name=${slurm_job_name} \
        -o ${slurm_out} \
        -e ${slurm_err} \
        ${PWD}/extract_chroms.slurm.sh ${chrom} ${in_fn} ${out_fn}
done
