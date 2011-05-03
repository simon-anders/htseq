import HTSeq
import itertools
import numpy

def GenomicInterval_xrange( gi, step ):
   for pos in xrange( gi.start, gi.end, step ):
      yield HTSeq.GenomicPosition( gi.chrom, pos, gi.strand, gi.genome )
      
def GenomicInterval_xranged( gi, step ):
   if gi.strand == "-":
      step *= -1
   for pos in xrange( gi.start_d, gi.end_d, step ):
      yield HTSeq.GenomicPosition( gi.chrom, pos, gi.strand, gi.genome ) 
                   
def ChromVector_steps( cv ):
   if isinstance( cv.array, numpy.ndarray ):
      start = cv.iv.start
      prev_val = None
      for i in xrange( cv.iv.start, cv.iv.end ):
         val = cv.array[ i - cv.offset ]
         if prev_val is None or val != prev_val:
            if prev_val is not None:
               yield ( HTSeq.GenomicInterval( cv.iv.chrom, start, i, cv.iv.strand ), prev_val )
            prev_val = val
            start = i
      yield ( HTSeq.GenomicInterval( cv.iv.chrom, start, cv.iv.end, cv.iv.strand ), prev_val )
   elif isinstance( cv.array, HTSeq.StepVector.StepVector ):
      for start, stop, value in cv.array[cv.iv.start:cv.iv.end].get_steps():
         yield ( HTSeq.GenomicInterval( cv.iv.chrom, start, stop, cv.iv.strand ), value )
   else:
      raise SystemError, "Unknown array type."
      
def GenomicArray_steps( ga ):
   for a in ga.chrom_vectors.values():
      for cv in a.values():
         for iv, val in cv.steps():
            yield iv, val

