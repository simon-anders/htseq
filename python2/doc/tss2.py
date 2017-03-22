import HTSeq
import numpy
from matplotlib import pyplot

sortedbamfile = HTSeq.BAM_Reader( "SRR001432_head_sorted.bam" )
gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.56_chrom1.gtf" )
halfwinwidth = 3000
fragmentsize = 200

tsspos = set()
for feature in gtffile:
   if feature.type == "exon" and feature.attr["exon_number"] == "1":
      tsspos.add( feature.iv.start_d_as_pos )

profile = numpy.zeros( 2*halfwinwidth, dtype='i' )   
for p in tsspos:
   window = HTSeq.GenomicInterval( p.chrom, 
       p.pos - halfwinwidth - fragmentsize, p.pos + halfwinwidth + fragmentsize, "." )
   for almnt in sortedbamfile[ window ]:
      almnt.iv.length = fragmentsize
      if p.strand == "+":
         start_in_window = almnt.iv.start - p.pos + halfwinwidth 
         end_in_window   = almnt.iv.end   - p.pos + halfwinwidth 
      else:
         start_in_window = p.pos + halfwinwidth - almnt.iv.end
         end_in_window   = p.pos + halfwinwidth - almnt.iv.start
      start_in_window = max( start_in_window, 0 )
      end_in_window = min( end_in_window, 2*halfwinwidth )
      if start_in_window >= 2*halfwinwidth or end_in_window < 0:
         continue
      profile[ start_in_window : end_in_window ] += 1

