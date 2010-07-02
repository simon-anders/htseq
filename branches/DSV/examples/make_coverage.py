# This example demonstrates how to read in the output of an aligner and
# produce a wiggle file.

fasta_file = "Dmel_BDGP5.4.55.fa"
aligner_file = "some_data.bowtie_out"
track_file = "some_data.wig"

import sys

from HTSeq import *

chromlens = dict( (seq.name, len(seq)) 
   for seq in FastaReader( fasta_file ) )

coverage = GenomicArray( chromlens, False, 'i' )

bwt = BowtieReader( aligner_file )

for read in bwt:

   # Extend each read to the fragment length, say 200
   read.iv.length = 200
   # In case that we went over the chromosome boundaries with this
   # we better correct for it.
   if read.iv.start < 0:
      read.iv.start = 0
   if read.iv.end >= chromlens[ read.iv.chrom ]:
      read.iv.end = chromlens[ read.iv.chrom ] - 1

   coverage.add_value( 1, read.iv )
   
coverage.write_bedgraph_file( sys.argv[1] + ".wig" )   
