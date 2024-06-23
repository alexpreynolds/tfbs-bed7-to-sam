#!/usr/bin/env python

import sys
import pyfaidx

in_fn = sys.argv[1]
out_fn = sys.argv[2]

seq_db = pyfaidx.Fasta('/net/seq/data/genomes/human/hg38/july/hg38.fa', strict_bounds=True)

with open(in_fn, 'r') as ifh, open(out_fn, 'w') as ofh:
  for line in ifh:
    in_elems = line.split('\t')
    chrom = in_elems[0]
    start = int(in_elems[1])
    end = int(in_elems[2])
    name = '{}%%{}'.format(in_elems[3], in_elems[6])
    score = float(in_elems[4])
    strand = in_elems[5]
    seq = seq_db[chrom][(start + 1):end].seq
    out_elems = [chrom, str(start), str(end), name, str(score), strand, seq]
    ofh.write('{}\n'.format('\t'.join(out_elems)))
