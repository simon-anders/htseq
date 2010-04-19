cimport numpy

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
   cpdef object add_bases_to_count_array( Sequence self, numpy.ndarray count_array_ )
   cpdef Sequence trim_left_end( Sequence self, Sequence pattern, float mismatch_prop = ? )   
   cpdef Sequence trim_right_end( Sequence self, Sequence pattern, float mismatch_prop = ? )   
   

cdef class SequenceWithQualities( Sequence ):
   cdef readonly str _qualstr
   cdef readonly str _qualstr_phred
   cdef readonly str _qualscale
   cdef readonly object _qualarr
   cdef _fill_qual_arr( SequenceWithQualities self )
   cpdef object add_qual_to_count_array( SequenceWithQualities self, numpy.ndarray count_array_ )
   cpdef SequenceWithQualities trim_left_end_with_quals( SequenceWithQualities self, 
         Sequence pattern, int max_mm_qual = ? )
   cpdef SequenceWithQualities trim_right_end_with_quals( SequenceWithQualities self, 
         Sequence pattern, int max_mm_qual = ? )

      
cdef class Alignment( object ):
   pass
   
cdef class AlignmentWithSequenceReversal( Alignment ):   
   cdef public SequenceWithQualities read_as_aligned
   cdef public SequenceWithQualities _read_as_sequenced
   cdef public GenomicInterval iv

cdef class SAM_Alignment( AlignmentWithSequenceReversal ):
   cdef public int flags
   cdef public list cigar
   cdef public int aQual
   cdef str _tags
   cdef public GenomicPosition mate_start
   cdef public str pe_which
   cdef public int inferred_insert_size
 
