import HTSeq
import itertools

def GenomicInterval_xrange( gi, step ):
   for pos in xrange( gi.start, gi.end, step ):
      yield HTSeq.GenomicPosition( gi.chrom, pos, gi.strand, gi.genome )
      
def GenomicInterval_xranged( gi, step ):
   if gi.strand == "-":
      step *= -1
   for pos in xrange( gi.start_d, gi.end_d, step ):
      yield HTSeq.GenomicPosition( gi.chrom, pos, gi.strand, gi.genome ) 
      
def GenomicArray_get_steps_convert_iv( step_iter, chrom, strand ):
   for start, end, value in step_iter:
      yield HTSeq.GenomicInterval( chrom, start, end, strand ), value
      
def Genomic_array_get_all_steps( self, values_only = False ):
   if values_only:
      for chrom in self.step_vectors:
	 if self.stranded:
	    for a in self.step_vectors[chrom]["+"].get_steps( values_only=True):
	       yield a
	    for a in self.step_vectors[chrom]["-"].get_steps( values_only=True):
	       yield a
	 else:
	    for a in self.step_vectors[chrom].get_steps( values_only=True):
	       yield a
   else:
      for chrom in self.step_vectors:
	 if self.stranded:
            for a in GenomicArray_get_steps_convert_iv( 
        	  self.step_vectors[chrom]["+"].get_steps(), chrom, "+" ):
               yield a         
            for a in GenomicArray_get_steps_convert_iv( 
        	  self.step_vectors[chrom]["-"].get_steps(), chrom, "-" ):
               yield a         
	 else:
            for a in GenomicArray_get_steps_convert_iv( 
        	  self.step_vectors[chrom].get_steps(), chrom, "." ):
               yield a
