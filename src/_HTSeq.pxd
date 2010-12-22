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
   cdef public bytes seq
   cdef public str name
   cdef public str descr
   cpdef Sequence get_reverse_complement( self )      
   cpdef object add_bases_to_count_array( Sequence self, numpy.ndarray count_array_ )
   cpdef Sequence trim_left_end( Sequence self, Sequence pattern, float mismatch_prop = ? )   
   cpdef Sequence trim_right_end( Sequence self, Sequence pattern, float mismatch_prop = ? )   
   

cdef class SequenceWithQualities( Sequence ):
   cdef readonly bytes _qualstr
   cdef readonly bytes _qualstr_phred
   cdef readonly str _qualscale
   cdef readonly object _qualarr
   cdef _fill_qual_arr( SequenceWithQualities self )
   cpdef object add_qual_to_count_array( SequenceWithQualities self, numpy.ndarray count_array_ )
   cpdef SequenceWithQualities trim_left_end_with_quals( SequenceWithQualities self, 
         Sequence pattern, int max_mm_qual = ? )
   cpdef SequenceWithQualities trim_right_end_with_quals( SequenceWithQualities self, 
         Sequence pattern, int max_mm_qual = ? )

      
cdef class Alignment( object ):
   cdef public SequenceWithQualities _read
   cdef public GenomicInterval iv
   
cdef class AlignmentWithSequenceReversal( Alignment ):   
   cdef public SequenceWithQualities read_as_aligned
   cdef public SequenceWithQualities _read_as_sequenced

cdef class SAM_Alignment( AlignmentWithSequenceReversal ):
   cdef public list cigar
   cdef public int aQual
   cdef list _optional_fields
   cdef public GenomicPosition mate_start
   cdef public str pe_which
   cdef public int inferred_insert_size
   cdef public bint proper_pair
   cdef public bint not_primary_alignment
   cdef public bint failed_platform_qc
   cdef public bint pcr_or_optical_duplicate
   cdef readonly str original_sam_line
   
 
