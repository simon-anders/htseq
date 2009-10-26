cdef enum strand_enum:
   strand_nostrand = 0
   strand_plus = 1
   strand_minus = 2

cdef class Strand:
   cdef strand_enum se

cdef class GenomicInterval:

   cdef public str chrom
   cdef public long start
   cdef public long end
   cdef Strand _strand
   cdef public object genome

   cdef Strand _get_strand( self )

   cpdef is_contained_in( self, GenomicInterval iv )
   cpdef contains( self, GenomicInterval iv )
   cpdef overlaps( self, GenomicInterval iv )
   cpdef extend_to_include( GenomicInterval self, GenomicInterval iv )


cdef class GenomicPosition( GenomicInterval ):
   pass
