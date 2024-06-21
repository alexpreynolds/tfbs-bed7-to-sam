#!/bin/bash

chrom=${1}
in_fn=${2}
out_fn=${3}

module add bedops

bedextract ${chrom} ${in_fn} > ${out_fn}

