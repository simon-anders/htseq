
cdef class GenomicInterval:

   cdef public str chrom
   cdef public long start
   cdef public long end
   cdef str _strand

   cpdef is_contained_in( self, GenomicInterval iv )
   cpdef contains( self, GenomicInterval iv )
   cpdef overlaps( self, GenomicInterval iv )
   cpdef extend_to_include( GenomicInterval self, GenomicInterval iv )

cdef class GenomicPosition( GenomicInterval ):
   pass


cdef class Sequence( object ):
   cdef public str seq
   cdef public str name
   cdef public str descr
   cpdef Sequence get_reverse_complement( self )      

cdef class SequenceWithQualities( Sequence ):
   cdef readonly str _qualstr
   cdef readonly str _qualscale
   cdef readonly object _qualarr

      
cdef class Alignment( object ):
   pass
   
cdef class AlignmentWithSequenceReversal( Alignment ):   
   cdef public SequenceWithQualities read_as_aligned
   cdef public SequenceWithQualities _read_as_sequenced
   cdef public GenomicInterval iv

cdef class SAM_Alignment( AlignmentWithSequenceReversal ):
   cdef public list cigar
   cdef public int aQual
   cdef str _tags
   
      
