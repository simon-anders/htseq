import HTSeq
import numpy
from matplotlib import pyplot

bamfile = HTSeq.BAM_Reader( "SRR001432_head.bam" )
gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.56_chrom1.gtf" )
halfwinwidth = 3000
fragmentsize = 200

tsspos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
for feature in gtffile:
   if feature.type == "exon" and feature.attr["exon_number"] == "1":
      p = feature.iv.start_d_as_pos
      window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
      tsspos[ window ] += p

profile = numpy.zeros( 2*halfwinwidth, dtype="i" )
for almnt in bamfile:
   if almnt.aligned:
      almnt.iv.length = fragmentsize
      s = set()
      for step_iv, step_set in tsspos[ almnt.iv ].steps():
         s |= step_set
      for p in s:
         if p.strand == "+":
            start_in_window = almnt.iv.start - p.pos + halfwinwidth
            end_in_window   = almnt.iv.end   - p.pos + halfwinwidth
         else:
            start_in_window = p.pos + halfwinwidth - almnt.iv.end
            end_in_window   = p.pos + halfwinwidth - almnt.iv.start
         start_in_window = max( start_in_window, 0 )
         end_in_window = min( end_in_window, 2*halfwinwidth )
         profile[ start_in_window : end_in_window ] += 1

