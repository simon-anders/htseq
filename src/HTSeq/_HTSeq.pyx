import sys
import os
import math
import re
import csv
import gzip
import itertools 
import collections
import cStringIO
import warnings

import numpy
cimport numpy

import StepVector
import _HTSeq_internal


###########################
##   GenomicInterval 
###########################

cdef str strand_plus = intern( "+" )
cdef str strand_minus = intern( "-" )
cdef str strand_nostrand = intern( "." )

cdef class GenomicInterval:

   """A GenomicInterval specifies an interval (i.e., a range of 
   consecutive positions) on a reference genome.
   
   A GenomicInterval object has the following slots, some of which 
   are calculated from the other:
     
      chrom: The name of a sequence (i.e., chromosome, contig, or 
         the like). If 'genome' is given, it should know 'chrom'.
      start: The start of the interval. Even on the reverse strand,
         this is always the smaller of the two values 'start' and 'end'.
         Note that all positions should be given as 0-based value!
      end: The end of the interval. Following Python convention for 
         ranges, this in one more than the coordinate of the last base
         that is considered part of the sequence.
      strand: The strand, as a single character, '+' or '-'. '.' indicates
         that the strand is irrelavant. (Alternatively, pass a Strand object.)
      length: The length of the interval, i.e., end - start
      start_d: The "directional start" position. This is the position of the
        first base of the interval, taking the strand into account. Hence, 
        this is the same as 'start' except when strand == '-', in which 
        case it is end-1.
      end_d: The "directional end": Usually, the same as 'end', but for 
        strand=='-1', it is start+1.
   """
   
   def __init__( GenomicInterval self, str chrom, long start, long end, 
         str strand = strand_nostrand ):
      """See the class docstring for the meaning of the slots. Note that 
      there is also a factory function, 'from_directional', to be used if
      you wish to specify start_d and length.
      """
      self.chrom = intern( chrom )
      self.start = start
      self.end = end
      self.strand = strand
      if self.start > self.end:
         raise ValueError, "start is larger than end"
            
   property strand:
      def __set__( self, strand ):
         strand = intern( strand )
         if not( strand is strand_plus or strand is strand_minus or 
               strand is strand_nostrand ):
            raise ValueError, "Strand must be'+', '-', or '.'."
         self._strand = strand
      def __get__( self ):
         return self._strand

   def __reduce__( GenomicInterval self ):
      return GenomicInterval, ( self.chrom, self.start, self.end, 
         self.strand )
         
   def __copy__( self ):
      constr, args = self.__reduce__()
      return constr( *args )   
   
   def __repr__( GenomicInterval self ):
      return "<%s object '%s', [%d,%s), strand '%s'>" % \
         ( self.__class__.__name__, self.chrom, self.start, 
           str(self.end) if self.end != sys.maxint else "Inf", self.strand )
         
   def __str__( GenomicInterval self ):
         return "%s:[%d,%s)/%s" % \
            ( self.chrom, self.start, str(self.end) if self.end != sys.maxint else "Inf", self.strand )

   property length:

      """The length is calculated as end - start. If you set the length, 
      'start_d' will be preserved, i.e., 'end' is changed, unless the strand
      is '-', in which case 'start' is changed."""   
            
      def __get__( GenomicInterval self ):
         return self.end - self.start

      def __set__( GenomicInterval self, long newLength ):
         if self._strand is not strand_minus:
            self.end = self.start + newLength
         else:
            self.start = self.end - newLength

   property start_d:
      """See the class docstring for the meaning of the 'directional start'.
      Note that if you set 'start_d', both the start and the end are changed, 
      such the interval gets the requested new directional start and its
      length stays unchanged."""

      def __get__( GenomicInterval self ):
         if self._strand is not strand_minus:
            return self.start
         else:
            return self.end - 1

      def __set__( GenomicInterval self, long newStartd ):
         if self._strand is not strand_minus:
            self.end = newStartd + self.length
            self.start = newStartd
         else:
            self.start = newStartd + 1 - self.length
            self.end = newStartd + 1
         
   property end_d:
   
      def __get__( GenomicInterval self ):
         if self._strand is not strand_minus:
            return self.end
         else:
            return self.start + 1
            
   property start_as_pos:
      def __get__( GenomicInterval self ):
         return GenomicPosition( self.chrom, self.start, self. strand )            

   property end_as_pos:
      def __get__( GenomicInterval self ):
         return GenomicPosition( self.chrom, self.end, self. strand )            

   property start_d_as_pos:
      def __get__( GenomicInterval self ):
         return GenomicPosition( self.chrom, self.start_d, self. strand )            

   property end_d_as_pos:
      def __get__( GenomicInterval self ):
         return GenomicPosition( self.chrom, self.end_d, self. strand )            

         
   def __richcmp__( GenomicInterval self, GenomicInterval other, int op ):
      if op == 2:  # ==
         if other == None:
            return False
         return self._strand is other._strand and \
            self.start == other.start and self.end == other.end
      elif op == 3:  # !=
         return not ( self == other )
      else:
         raise NotImplementedError
         
   def __hash__( GenomicInterval self ):
      return hash( ( self.chrom, self.start, self.end, self.strand ) )
           
   cpdef is_contained_in( GenomicInterval self, GenomicInterval iv ):
      """Returns a boolean value indicating whether the 'self' interval 
      is fully within the 'iv' interval.
      
      This is deemed the case if
        - both are on the same chromosome, and    
        - both are on the same strand, or at least one of them is
           not stranded (i.e., has strand == '.'), and
        - self.start >= iv.start, and
        - self.end <= iv.end
      """
      if iv == None:
         return False
      if self.chrom != iv.chrom:
         return False
      if self._strand is not strand_nostrand and iv.strand is not strand_nostrand and \
            self.strand is not iv._strand:
         return False 
      if self.start < iv.start or self.end > iv.end:
         return False
      return True

   cpdef contains( GenomicInterval self, GenomicInterval iv ):
      """Returns a boolean value indicating whether the 'self' interval 
      fully contains the 'iv' interval.

      See 'is_contained_in' for the exact criteria.
      """
      if iv == None:
            return False
      return iv.is_contained_in( self )

   cpdef overlaps( GenomicInterval self, GenomicInterval iv ):
      """Returns a boolean value indicating whether the 'self' interval 
      overlaps the 'iv' interval.
      
      This is deemed the case if
        - both are on the same chromosome, and    
        - both are on the same strand, or at least one of them is
           not stranded (i.e., has strand == '.'), and
        - the actual intervals overlap
      """
      if iv == None:
         return False
      if self.chrom != iv.chrom:
         return False
      if self.strand is not strand_nostrand and iv.strand is not strand_nostrand and \
            self.strand is not iv.strand:
         return False 
      if self.start <= iv.start:
         return self.end > iv.start
      else: 
         return iv.end > self.start
      
   def xrange( GenomicInterval self, long int step = 1 ):
      """Generate an iterator over the GenomicPositions covered by the interval,
      running from start to end.
      """
      return _HTSeq_internal.GenomicInterval_xrange( self, step )

   def xrange_d( GenomicInterval self, long int step = 1 ):
      """Generate an iterator over the GenomicPositions covered by the interval.
      running from start_d to end_d.
      """
      return _HTSeq_internal.GenomicInterval_xranged( self, step )
      
   cpdef extend_to_include( GenomicInterval self, GenomicInterval iv ):
      """Extend the interval such that it includes iv."""
      if iv is None:
         raise TypeError, "Cannot extend an interval to include None."
      if self.chrom != iv.chrom:
         raise ValueError, "Cannot extend an interval to include an interval on another chromosome."
      if self.strand.se is not strand_nostrand and iv.strand is not strand_nostrand and \
            self.strand is not iv.strand:
         raise ValueError, "Cannot extend an interval to include an interval on another strand."
      self.start = min( self.start, iv.start )
      self.end = max( self.end, iv.end )
      
   def copy( self ):
      return GenomicInterval( self.chrom, self.start, self.end, self.strand )


def GenomicInterval_from_directional( str chrom, long int start_d, long int length, str strand="." ):
   strand = intern( strand )
   if strand.se is not strand_minus:
      return GenomicInterval( chrom, start_d, start_d+length, strand )
   else:
      return GenomicInterval( chrom, start_d-length+1, start_d+1, strand )


cdef class GenomicPosition( GenomicInterval ):

   """A GenomicPosition specifies the position of a nucleotide or
   base pair on a reference genome.
   
   It has the following slots:
      chrom: The name of a sequence (i.e., chromosome, contig, or 
         the like). 
      pos: The position on the sequence specified by seqname.
         The position should always be given as 0-based value!
      strand: The strand, as a single character, '+' or '-'. '.' indicates
         that the strand is irrelavant.

   The GenomicPosition class is derived from GenomicInterval. Hence,
   a GenomicPosition is always a GenomicInterval of length 1. Do not tinker
   with the exposed GenomeInterval slots.
   """

   def __init__( self, str chrom, long int pos, str strand='.' ):
      GenomicInterval.__init__( self, chrom, pos, pos+1, strand )
      
   property pos:
   
      """As GenomicPosition is a subclass of GenomicInterval, 'pos' is actually
      just an alias for 'start_d'.
      """
   
      def __get__( self ):
         return self.start_d
         
      def __set__( self, long newValue ):
         self.start_d = newValue
      
   property end:
   
      def __get__( self ):
         return self.start + 1

   property length:
   
      def __get__( self ):
         return 1
      
   def __repr__( self ):
      return "<%s object '%s':%d, strand '%s'>" % \
         ( self.__class__.__name__, self.chrom, self.pos, self.strand )               
            
   def __str__( self ):
      return "%s:%d/%s" % ( self.chrom, self.pos, self.strand )            

   def __reduce__( GenomicPosition self ):
      return GenomicPosition, ( self.chrom, self.pos, self.strand )
   
   def copy( self ):
      return GenomicPosition( self.chrom, self.pos, self.strand )
      

cdef class ChromVector( object ):

   cdef public object array 
   cdef public GenomicInterval iv
   cdef public int offset
   cdef public bint is_vector_of_sets
   cdef public str _storage

   @classmethod 
   def create( cls, GenomicInterval iv, str typecode, str storage, str memmap_dir = "" ):
      ncv = cls()
      ncv.iv = iv
      if storage == "ndarray":
         if typecode != 'O':
            ncv.array = numpy.zeros( shape = ( iv.length, ), dtype = typecode )
         else:
            ncv.array = numpy.empty( shape = ( iv.length, ), dtype = typecode )
            ncv.array[:] = None
      elif storage == "memmap":
         ncv.array = numpy.memmap( shape = ( iv.length, ), dtype = typecode, 
            filename = os.path.join( memmap_dir, iv.chrom + iv.strand + ".nmm" ), mode='w+' )
      elif storage == "step":
         ncv.array = StepVector.StepVector.create( typecode = typecode )
      else:
         raise ValueError, "Illegal storage mode."
      ncv._storage = storage
      # TODO: Test whether offset works properly
      ncv.offset = iv.start
      ncv.is_vector_of_sets = False
      return ncv
   
   @classmethod 
   def _create_view( cls, ChromVector vec, GenomicInterval iv ):
      v = cls()
      v.iv = iv
      v.array = vec.array
      v.offset = vec.offset
      v.is_vector_of_sets = vec.is_vector_of_sets
      v._storage = vec._storage      
      return v

   def __getitem__( self, index ):
      cdef slice index_slice
      cdef long int index_int
      cdef long int start, stop
      cdef GenomicInterval iv
      if isinstance( index, int ):
         index_int = index
         if index_int < self.iv.start or index_int >= self.iv.end:
            raise IndexError
         return self.array[ index_int - self.offset ]
      elif isinstance( index, slice ):
         index_slice = index
         if index_slice.start is not None:
            start = index_slice.start
            if start < self.iv.start:
               raise IndexError, "start too small"
         else:
            start = self.iv.start
         if index_slice.stop is not None:
            stop = index_slice.stop
            if stop > self.iv.end:
               raise IndexError, "stop too large"
         else:
            stop = self.iv.end
         iv = GenomicInterval( self.iv.chrom, start, stop, self.iv.strand )
         if not self.iv.contains( iv ):
            raise IndexError
         return ChromVector._create_view( self, iv )
      elif isinstance( index, GenomicInterval ):
         if not self.iv.contains( index ):
            raise IndexError
         if self.iv.strand is strand_nostrand and \
               index.strand is not strand_nostrand:
            iv = iv.copy()
            iv.strand = strand_nostrand
         return ChromVector._create_view( self, iv )
      else:
         raise TypeError, "Illegal index type"
   
   def __setitem__( self, index, value ):
      cdef slice index_slice
      cdef long int start, stop
      if isinstance( value, ChromVector ):       
         if self.array is value.array and value.iv.start == index.start and \
               value.iv.end == index.stop and ( index.step is None or index.step == 1 ):
            return
         else:
            raise NotImplementedError, "Required assignment signature not yet implemented."
      if isinstance( index, int ):
         self.array[ index - self.iv.start ] = value
      elif isinstance( index, slice ):
         index_slice = index
         if index_slice.start is not None:
            start = index_slice.start
            if start < self.iv.start:
               raise IndexError, "start too small"
         else:
            start = self.iv.start
         if index_slice.stop is not None:
            stop = index_slice.stop
            if stop > self.iv.end:
               raise IndexError, "stop too large"
         else:
            stop = self.iv.end
         self.array[ start - self.offset : stop - self.iv.start : index.step ] = value
      elif isinstance( index, GenomicInterval ):
         if index.chrom != self.iv.chrom:
            raise KeyError, "Chromosome name mismatch."
         if self.iv.strand is not strand_nostrand and \
               self.iv.strand is not self.index.strand:
            raise KeyError, "Strand mismatch."
         self.array[ index.iv.start - self.iv.start, 
            index.iv.end - self.iv.start ] = value
      else:
         raise TypeError, "Illegal index type"

   def __iadd__( self, value ):
      if not self.is_vector_of_sets:
         self.array[ self.iv.start - self.offset : self.iv.end - self.offset ].__iadd__( value )
      else:
         def addval( x ):
            y = x.copy()
            y.add( value )
            return y
         self.apply( addval )
      return self
      
   def __iter__( self ):
      return self.values()
      
   def values( self ):
      return iter( self.array[ self.iv.start - self.offset : self.iv.end - self.offset ] )
   
   def steps( self ):
      return _HTSeq_internal.ChromVector_steps( self )
      
   def apply( self, fun ):
      for iv, value in self.steps():
         self.array[ iv.start - self.offset : iv.end - self.offset ] = fun( value )
         
   def __repr__( self ):
      return "<%s object, %s, %s>" % ( self.__class__.__name__, str(self.iv), self._storage )
         
   def __reduce__( self ):
      assert self.__class__ is ChromVector
      return( _ChromVector_unpickle, 
         ( self.array, self.iv, self.offset, self.is_vector_of_sets, self._storage ) )
         
def _ChromVector_unpickle( array, iv, offset, is_vector_of_sets, _storage ):
   cv = ChromVector()
   cv.array =  array 
   cv.iv = iv
   cv.offset = offset
   cv.is_vector_of_sets = is_vector_of_sets
   cv._storage = _storage
   return cv
   
cdef class GenomicArray( object ):
   
   cdef public dict chrom_vectors
   cdef readonly bint stranded
   cdef readonly str typecode
   cdef public bint auto_add_chroms
   cdef readonly str storage
   cdef readonly str memmap_dir
   
   def __init__( self, object chroms, bint stranded=True, str typecode='d',
         str storage='step', str memmap_dir = "" ):
      self.chrom_vectors = {}
      self.stranded = stranded
      self.typecode = typecode
      self.auto_add_chroms = chroms == "auto"
      if self.auto_add_chroms:      
         chroms = []
         if storage != 'step':
            raise TypeError, "Automatic adding of chromosomes can " + \
               " only be used with storage type 'StepVector'."
      elif isinstance( chroms, list ):
         if storage != 'step':
            raise TypeError, "Indefinite-length chromosomes can " + \
               " only be used with storage type 'StepVector'."
         chroms = dict( [ ( c, sys.maxint ) for c in chroms ] )
      elif not isinstance( chroms, dict ):
         raise TypeError, "'chroms' must be a list or a dict or 'auto'."
      self.storage = storage
      self.memmap_dir = memmap_dir
      
      for chrom in chroms:
         self.add_chrom( chrom, chroms[chrom] )
            
   def __getitem__( self, index ):
      if isinstance( index, GenomicInterval ):
         if self.stranded and index.strand not in ( strand_plus, strand_minus ):
            raise KeyError, "Non-stranded index used for stranded GenomicArray."
         if self.auto_add_chroms and index.chrom not in self.chrom_vectors:
            self.add_chrom( index.chrom )
         if isinstance( index, GenomicPosition ):
            if self.stranded:
               return self.chrom_vectors[ index.chrom ][ index.strand ][ index.pos ]
            else:
               return self.chrom_vectors[ index.chrom ][ strand_nostrand ][ index.pos ]
         else:
            if self.stranded:
               return self.chrom_vectors[ index.chrom ][ index.strand ][ index.start : index.end ]
            else:
               return self.chrom_vectors[ index.chrom ][ strand_nostrand ][ index.start : index.end ]
      else:
         return self.chrom_vectors[ index ]

   def __setitem__( self, index, value ):
      cdef GenomicInterval index2
      if isinstance( value, ChromVector ): 
         if not isinstance( index, GenomicInterval ):
            raise NotImplementedError, "Required assignment signature not yet implemented."
         index2 = index.copy()
         if not self.stranded:
            index2.strand = strand_nostrand
         if self.chrom_vectors[ index2.chrom ][ index2.strand ].array is value.array and index2 == value.iv:
            return
         raise NotImplementedError, "Required assignment signature not yet implemented."
      if isinstance( index, GenomicInterval ):
         if self.stranded and index.strand not in ( strand_plus, strand_minus ):
            raise KeyError, "Non-stranded index used for stranded GenomicArray."
         if self.auto_add_chroms and index.chrom not in self.chrom_vectors:
            self.add_chrom( index.chrom )
         if self.stranded:
            self.chrom_vectors[ index.chrom ][ index.strand ][ index.start : index.end ] = value
         else:
            self.chrom_vectors[ index.chrom ][ strand_nostrand ][ index.start : index.end ] = value
      else:
         raise TypeError, "Illegal index type."
            
   def add_chrom( self, chrom, length = sys.maxint, start_index = 0 ):
      cdef GenomicInterval iv 
      if length == sys.maxint:
         iv = GenomicInterval( chrom, start_index, sys.maxint, "." )
      else:
         iv = GenomicInterval( chrom, start_index, start_index + length, "." )
      if self.stranded:
         self.chrom_vectors[ chrom ] = {}
         iv.strand = "+"
         self.chrom_vectors[ chrom ][ strand_plus ] = \
            ChromVector.create( iv, self.typecode, self.storage, self.memmap_dir )
         iv = iv.copy()
         iv.strand = "-"
         self.chrom_vectors[ chrom ][ strand_minus ] = \
            ChromVector.create( iv, self.typecode, self.storage, self.memmap_dir )
      else:   
         self.chrom_vectors[ chrom ] = {
            strand_nostrand:  ChromVector.create( iv, self.typecode, self.storage ) }
   
   def __reduce__( self ):
      return ( _GenomicArray_unpickle, ( self.stranded, self.typecode, self.chrom_vectors ) )
      
   def write_bedgraph_file( self, file_or_filename, strand=".", track_options="" ):
      if ( not self.stranded ) and strand != ".":
         raise ValueError, "Strand specified in unstranded GenomicArray."
      if self.stranded and strand not in ( strand_plus, strand_minus ):
         raise ValueError, "Strand must be specified for stranded GenomicArray."
      if hasattr( file_or_filename, "write" ):
         f = file_or_filename
      else:
         f = open( file_or_filename, "w" )
      if track_options == "":
         f.write( "track type=bedGraph\n" )
      else:
         f.write( "track type=bedGraph %s\n" % track_options )
      for chrom in self.chrom_vectors:
         for iv, value in self.chrom_vectors[ chrom ][ strand ].steps():
            if iv.start == -sys.maxint-1 or iv.end == sys.maxint:
               continue
            f.write( "%s\t%d\t%d\t%f\n" % (iv.chrom, iv.start, iv.end, value) )    
      if not hasattr( file_or_filename, "write" ):
         f.close()
         
   def steps( self ):
      return _HTSeq_internal.GenomicArray_steps( self )
      
       
def _GenomicArray_unpickle( stranded, typecode, chrom_vectors ):
   ga = GenomicArray( {}, stranded, typecode )
   ga.chrom_vectors = chrom_vectors
   return ga
   
   
   
###########################
##   Sequences
###########################
   
   
def _make_translation_table_for_complementation( ):
   t = [ chr(i) for i in xrange(256) ]
   t[ ord('A') ] = 'T'
   t[ ord('T') ] = 'A'
   t[ ord('C') ] = 'G'
   t[ ord('G') ] = 'C'
   t[ ord('a') ] = 't'
   t[ ord('t') ] = 'a'
   t[ ord('c') ] = 'g'
   t[ ord('g') ] = 'c'
   return ''.join( t )
   
cdef bytes _translation_table_for_complementation = _make_translation_table_for_complementation( )

cpdef bytes reverse_complement( bytes seq ):
   """Returns the reverse complement of DNA sequence 'seq'. Does not yet
   work with extended IUPAC nucleotide letters or RNA."""

   return seq[ ::-1 ].translate( _translation_table_for_complementation )

base_to_column = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4 }

cdef class Sequence( object ):
   """A Sequence, typically of DNA, with a name.
   """
   
   def __init__( self, bytes seq, str name="unnamed" ):
      self.seq = seq
      self.name = name
      self.descr = None
      
   cpdef Sequence get_reverse_complement( self ):
      return Sequence(
         reverse_complement( self.seq ),
         "revcomp_of_" + self.name )
         
   def __str__( self ):
      return self.seq
      
   def __repr__( self ):
      return "<%s object '%s' (length %d)>" % ( 
         self.__class__.__name__, self.name, len( self.seq ) )

   def __len__( self ):
      return len( self.seq )
      
   def __getitem__( self, item ):
      if self.name.endswith( "[part]" ):
         new_name = self.name
      else:
         new_name = self.name + "[part]"
      return Sequence( self.seq[ item ], new_name )
   
   def write_to_fasta_file( self, fasta_file ):
      if self.descr is not None:
         fasta_file.write( ">%s %s\n" % ( self.name, self.descr ) )
      else:
         fasta_file.write( ">%s\n" % self.name )
      i = 0
      while i*70 < len(self.seq):
         fasta_file.write( self.seq[ i*70 : (i+1)*70 ] + "\n" )
         i += 1
         
   cpdef object add_bases_to_count_array( Sequence self, numpy.ndarray count_array_ ):
      
      cdef numpy.ndarray[ numpy.int_t, ndim=2 ] count_array = count_array_
      cdef int seq_length = len( self.seq )

      if numpy.PyArray_DIMS( count_array )[0] < seq_length:
         raise ValueError, "'count_array' too small for sequence."
      if numpy.PyArray_DIMS( count_array )[1] < 5:
         raise ValueError, "'count_array' has too few columns."
      
      cdef numpy.npy_intp i
      cdef char b
      cdef char* seq_cstr = self.seq
      for i in xrange( seq_length ):
         b = seq_cstr[i]
         if b == 'A' or b == 'a':
            count_array[ i, 0 ] += 1
         elif b == 'C' or b == 'c':
            count_array[ i, 1 ] += 1
         elif b == 'G' or b == 'g':
            count_array[ i, 2 ] += 1
         elif b == 'T' or b == 't':
            count_array[ i, 3 ] += 1
         elif b == 'N' or b == 'n' or b == ".":
            count_array[ i, 4 ] += 1
         else:
            raise ValueError, "Illegal base letter encountered."
         
      return None
      
   cpdef Sequence trim_left_end( Sequence self, Sequence pattern, float mismatch_prop = 0. ):
      cdef int seqlen = len( self.seq )
      cdef int patlen = len( pattern.seq )
      cdef int minlen
      if seqlen < patlen:
         minlen = seqlen
      else: 
         minlen = patlen
      cdef char * seq_cstr = self.seq
      cdef char * pat_cstr = pattern.seq
      cdef int match = 0
      cdef int i, j
      cdef int num_mismatches
      for i in xrange( 1, minlen+1 ):
         num_mismatches = 0
         for j in xrange( i ):
            if seq_cstr[ j ] != pat_cstr[ patlen - i + j ]:
               num_mismatches += 1
               if num_mismatches > mismatch_prop * i:
                  break
         else:
            match = i
      return self[ match : seqlen ]

   cpdef Sequence trim_right_end( Sequence self, Sequence pattern, float mismatch_prop = 0. ):
      cdef int seqlen = len( self.seq )
      cdef int patlen = len( pattern.seq )
      cdef int minlen
      if seqlen < patlen:
         minlen = seqlen
      else: 
         minlen = patlen
      cdef char * seq_cstr = self.seq
      cdef char * pat_cstr = pattern.seq
      cdef int match = 0
      cdef int i, j
      cdef int num_mismatches
      for i in xrange( 1, minlen+1 ):
         num_mismatches = 0
         for j in xrange( i ):
            if seq_cstr[ seqlen - i + j ] != pat_cstr[ j ]:
               num_mismatches += 1
               if num_mismatches > mismatch_prop * i:
                  break
         else:
            match = i
      return self[ 0 : seqlen-match ]


cdef class SequenceWithQualities( Sequence ):
   """A Sequence with base-call quality scores.
   It now has property  'qual', an integer NumPy array of Sanger/Phred 
   quality scores of the  base calls.
   """

   def __init__( self, bytes seq, str name, bytes qualstr, str qualscale="phred" ):
      """ Construct a SequenceWithQuality object.
      
        seq       - The actual sequence.
        name      - The sequence name or ID
        qualstr   - The quality string. Must have the same length as seq
        qualscale - The encoding scale of the quality string. Must be one of
                      "phred", "solexa", "solexa-old" )
      """
      if len( seq ) != len( qualstr ):
         raise ValueError, "'seq' and 'qualstr' do not have the same length."
      Sequence.__init__( self, seq, name )
      self._qualstr = qualstr
      self._qualscale = qualscale
      self._qualarr = None
      self._qualstr_phred = b''

   cdef _fill_qual_arr( SequenceWithQualities self ):
      cdef int seq_len = len( self.seq )
      if seq_len != len( self._qualstr ):
         raise ValueError, "Quality string has not the same length as sequence."
      cdef numpy.ndarray[ numpy.int_t, ndim=1 ] qualarr = numpy.empty( ( seq_len, ), numpy.int )
      cdef int i
      cdef char * qualstr = self._qualstr
      if self._qualscale == "phred":
         for i in xrange( seq_len ):
            qualarr[i] = qualstr[i] - 33
      elif self._qualscale == "solexa":
         for i in xrange( seq_len ):
            qualarr[i] = qualstr[i] - 64
      elif self._qualscale == "solexa-old":
         for i in xrange( seq_len ):
            qualarr[i] = 10 * math.log10( 1 + 10 ** ( qualstr[i] - 64 ) / 10.0 )
      else:
         raise ValueError, "Illegal quality scale '%s'." % self._qualscale
      self._qualarr = qualarr

   property qual:
      def __get__( self ):
         if self._qualarr is None:
            self._fill_qual_arr()
         return self._qualarr
      def __set__( self, newvalue ):
         if not ( isinstance( newvalue, numpy.ndarray ) and newvalue.dtype == numpy.int ) :
            raise TypeError, "qual can only be assigned a numpy array of type numpy.int"
         if not ( newvalue.shape == ( len(self.seq), ) ) :
            raise TypeError, "assignment to qual with illegal shape"
         self._qualarr = newvalue
         self._qualstr = b""
         self._qualscale = "none"
         self._qualstr_phred = b""
              
   def __repr__( self ):
      return "<%s object '%s'>" % ( self.__class__.__name__, self.name )

   def __getitem__( self, item ):
      if self.name.endswith( "[part]" ):
         new_name = self.name
      else:
         new_name = self.name + "[part]"
      return SequenceWithQualities( 
         self.seq[ item ], new_name, self.qualstr[ item ] )

   @property
   def qualstr( self ):
      cdef int seqlen
      cdef char * qualstr_phred_cstr = self._qualstr_phred
      cdef int i
      cdef numpy.ndarray[ numpy.int_t, ndim=1 ] qual_array
      if qualstr_phred_cstr[0] == 0:
         if self._qualscale == "phred":
            self._qualstr_phred = self._qualstr
         else:
            seqlen = len( self.seq )
            self._qualstr_phred = <bytes>(' ') * seqlen
            qualstr_phred_cstr = self._qualstr_phred
            if self._qualarr is None:
               self._fill_qual_arr()
            qual_array = self._qualarr               
            for i in xrange( seqlen ):
               qualstr_phred_cstr[i] = 33 + qual_array[i]            
      return self._qualstr_phred
      

   def write_to_fastq_file( self, fastq_file ):
      if hasattr( self, "descr" ) and self.descr is not None:
         fastq_file.write( "@%s %s\n" % ( self.name, self.descr ) )
      else:
         fastq_file.write( "@%s\n" % self.name )
      fastq_file.write( self.seq + "\n" )
      fastq_file.write( "+\n" )
      fastq_file.write( self.qualstr + "\n" )

   def get_fastq_str( self, bint convert_to_phred=False ):
      sio = cStringIO.StringIO()
      self.write_to_fastq_file( sio, convert_to_phred )
      return sio.getvalue()

   cpdef SequenceWithQualities get_reverse_complement( self ):
      cdef SequenceWithQualities res
      res = SequenceWithQualities(
         reverse_complement( self.seq ),
         "revcomp_of_" + self.name,
         self._qualstr[::-1],
         self._qualscale )
      if self._qualarr is not None:
         res._qualarr = self._qualarr[::-1]
      return res
      
   cpdef object add_qual_to_count_array( SequenceWithQualities self,
         numpy.ndarray count_array_ ):
      
      cdef numpy.ndarray[ numpy.int_t, ndim=2 ] count_array = count_array_
      if self._qualarr is None:
         self._fill_qual_arr()   
      cdef numpy.ndarray[ numpy.int_t, ndim=1 ] qual_array = self._qualarr

      cdef numpy.npy_intp seq_length = numpy.PyArray_DIMS( qual_array  )[0]
      cdef numpy.npy_intp qual_size  = numpy.PyArray_DIMS( count_array )[1]

      if seq_length > numpy.PyArray_DIMS( count_array )[0]:
         raise ValueError, "'count_array' too small for sequence."
      
      cdef numpy.npy_intp i
      cdef numpy.npy_int q
      for i in xrange( seq_length ):
         q = qual_array[i]
         if( q >= qual_size ):
            raise ValueError, "Too large quality value encountered."
         count_array[ i, q ] += 1
         
      return None
   
   cpdef SequenceWithQualities trim_left_end_with_quals( SequenceWithQualities self, 
         Sequence pattern, int max_mm_qual = 5 ):
      cdef int seqlen = len( self.seq )
      cdef int patlen = len( pattern.seq )
      cdef int minlen
      if seqlen < patlen:
         minlen = seqlen
      else: 
         minlen = patlen
      cdef char * seq_cstr = self.seq
      cdef char * pat_cstr = pattern.seq
      cdef int match = 0
      cdef int i, j
      cdef int sum_mm_qual
      if self._qualarr is None:
         self._fill_qual_arr()   
      cdef numpy.ndarray[ numpy.int_t, ndim=1 ] qual_array = self._qualarr
      for i in xrange( 1, minlen+1 ):
         num_mismatches = 0
         for j in xrange( i ):
            if seq_cstr[ j ] != pat_cstr[ patlen - i + j ]:
               sum_mm_qual += qual_array[ j ]
               if sum_mm_qual > max_mm_qual:
                  break
         else:
            match = i
      return self[ match : seqlen ]

   cpdef SequenceWithQualities trim_right_end_with_quals( SequenceWithQualities self, 
          Sequence pattern, int max_mm_qual = 5 ):
      cdef int seqlen = len( self.seq )
      cdef int patlen = len( pattern.seq )
      cdef int minlen
      if seqlen < patlen:
         minlen = seqlen
      else: 
         minlen = patlen
      cdef char * seq_cstr = self.seq
      cdef char * pat_cstr = pattern.seq
      cdef int match = 0
      cdef int i, j
      cdef int sum_mm_qual
      if self._qualarr is None:
         self._fill_qual_arr()   
      cdef numpy.ndarray[ numpy.int_t, ndim=1 ] qual_array = self._qualarr
      for i in xrange( 1, minlen+1 ):
         sum_mm_qual = 0
         for j in xrange( i ):
            if seq_cstr[ seqlen - i + j ] != pat_cstr[ j ]:
               sum_mm_qual += qual_array[ seqlen - i + j ]
               if sum_mm_qual > max_mm_qual:
                  break
         else:
            match = i
      return self[ 0 : seqlen-match ]

   
###########################
##   Alignment
###########################

cdef class Alignment( object ):

   """Alignment base type:

   An alignment object can be defined in different ways but will always
   provide these attributes:
     read:      a SequenceWithQualities object with the read
     aligned:   whether the read is aligned
     iv:        a GenomicInterval object with the alignment position 
   """
   
   def __init__( self, read, iv ):
      self._read = read
      self.iv = iv
      
   @property
   def read( self ):
      return self._read
   
   def __repr__( self ):
      cdef str s
      if self.paired_end:
         s = "Paired-end Read"
      else:
         s = "Read"
      if self.aligned:
         return "<%s object: %s '%s' aligned to %s>" % (
            self.__class__.__name__, s, self.read.name, str(self.iv) )
      else:
         return "<%s object: %s '%s', not aligned>" % (
            self.__class__.__name__, s, self.read.name )
            
   @property
   def paired_end( self ):
      return False
         
   @property
   def aligned( self ):
      """Returns True unless self.iv is None. The latter indicates that
      this record decribes a read for which no alignment was found.
      """
      return self.iv is not None

cdef class AlignmentWithSequenceReversal( Alignment ):

   """Many aligners report the read's sequence in reverse-complemented form
   when it was mapped to the reverse strand. For such alignments, a 
   daughter class of this one should be used.
   
   Then, the read is stored as aligned in the 'read_as_aligned' field,
   and get reverse-complemented back to the sequenced form when the 'read'
   attribute is sequenced.
   """

   def __init__( self, SequenceWithQualities read_as_aligned, GenomicInterval iv ):
      self.read_as_aligned = read_as_aligned      
      self._read_as_sequenced = None
      self.iv = iv
      
   property read:
      def __get__( self ):
         if self._read_as_sequenced is None:
            if (not self.aligned) or self.iv.strand != "-":
               self._read_as_sequenced = self.read_as_aligned
            else:
               self._read_as_sequenced = self.read_as_aligned.get_reverse_complement()
               self._read_as_sequenced.name = self.read_as_aligned.name
         return self._read_as_sequenced      
      #def __set__( self, read ):
      #   self.read_as_aligned = read
      #   self._read_as_sequenced = None
         
         
cdef class BowtieAlignment( AlignmentWithSequenceReversal ):

   """When reading in a Bowtie file, objects of the class BowtieAlignment
   are returned. In addition to the 'read' and 'iv' fields (see Alignment
   class), the fields 'reserved' and 'substitutions' are provided. These 
   contain the content of the respective columns of the Bowtie output 
   
   [A parser for the substitutions field will be added soon.]
   """
   
   cdef public str reserved
   cdef public str substitutions
   
   def __init__( self, bowtie_line ):
      cdef str readId, strand, chrom, position, read, qual
      cdef int positionint
      (readId, strand, chrom, position, read, qual, 
         self.reserved, self.substitutions) = bowtie_line.split( '\t' )
      positionint = int( position )
      AlignmentWithSequenceReversal.__init__( self, 
         SequenceWithQualities( read, readId, qual ),
         GenomicInterval( chrom, positionint, positionint + len(read), strand ) )
      
      
cigar_operation_names = {
   'M': 'matched',
   'I': 'inserted',
   'D': 'deleted',
   'N': 'skipped',
   'S': 'soft-clipped',
   'H': 'hard-clipped',
   'P': 'padded' }

cdef class CigarOperation( object ):
   
   cdef public str type
   cdef public int size
   cdef public GenomicInterval ref_iv
   cdef public int query_from, query_to
   
   def __init__( self, str type_, int size, int rfrom, int rto, int qfrom, 
         int qto, str chrom, str strand, bint check=True ):
      self.type = type_
      self.size = size
      self.ref_iv = GenomicInterval( chrom, rfrom, rto, strand )
      self.query_from = qfrom
      self.query_to = qto
      if check and not self.check():
         raise ValueError, "Inconsistent CIGAR operation."
               
   def __repr__( self ):
      return "< %s: %d base(s) %s on ref iv %s, query iv [%d,%d) >" % (
         self.__class__.__name__, self.size, cigar_operation_names[ self.type ],
         str( self.ref_iv ), self.query_from, self.query_to )
         
   def check( CigarOperation self ):
      cdef int qlen = self.query_to - self.query_from
      cdef int rlen = self.ref_iv.length
      if self.type == 'M':
         if not ( qlen == self.size and rlen == self.size ):
            return False
      elif self.type == 'I' or self.type == 'S':
         if not ( qlen == self.size and rlen == 0 ):
            return False
      elif self.type == 'D' or self.type == 'N':
         if not ( qlen == 0 and rlen == self.size ):
            return False
      elif self.type == 'H' or self.type == 'P':
         if not ( qlen == 0 and rlen == 0 ):
            return False
      else:
         return False
      return True

_re_cigar_codes = re.compile( '([A-Z])' )

cpdef list parse_cigar( str cigar_string, int ref_left = 0, str chrom = "", str strand = "." ):
   cdef list split_cigar, cl
   cdef int size
   cdef str code
   split_cigar = _re_cigar_codes.split( cigar_string )
   if split_cigar[-1] != '' or len(split_cigar) % 2 != 1:
      raise ValueError, "Illegal CIGAR string '%s'" % cigar_string
   cl = []
   for i in xrange( len(split_cigar) // 2 ):
      try:
         size = int( split_cigar[2*i] )
      except ValueError:
         raise ValueError, "Illegal CIGAR string '%s'" % cigar_string
      code  = split_cigar[2*i+1]
      cl.append( ( code, size ) )
   return build_cigar_list( cl, ref_left, chrom, strand  )

cpdef list build_cigar_list( list cigar_pairs, int ref_left = 0, str chrom = "", str strand = "." ):
   cdef list split_cigar, res
   cdef int rpos, qpos, size
   cdef str code
   rpos = ref_left
   qpos = 0
   res = []
   for code, size in cigar_pairs:
      if code == 'M':
         res.append( CigarOperation ( 
            'M', size, rpos, rpos + size, qpos, qpos + size, chrom, strand ) )
         rpos += size
         qpos += size
      elif code == 'I':
         res.append( CigarOperation ( 
            'I', size, rpos, rpos, qpos, qpos + size, chrom, strand ) )
         qpos += size
      elif code == 'D':
         res.append( CigarOperation ( 
            'D', size, rpos, rpos + size, qpos, qpos, chrom, strand ) )
         rpos += size
      elif code == 'N':
         res.append( CigarOperation ( 
            'N', size, rpos, rpos + size, qpos, qpos, chrom, strand ) )
         rpos += size
      elif code == 'S':
         res.append( CigarOperation ( 
            'S', size, rpos, rpos, qpos, qpos + size, chrom, strand ) )
         qpos += size
      elif code == 'H':
         res.append( CigarOperation ( 
            'H', size, rpos, rpos, qpos, qpos, chrom, strand ) )
      elif code == 'P':
         res.append( CigarOperation ( 
            'P', size, rpos, rpos, qpos, qpos, chrom, strand ) )
      else:
         raise ValueError, "Unknown CIGAR code '%s' encountered." % code
   return res
      
cdef _parse_SAM_optional_field_value( str field ):
   if len(field) < 5 or field[2] != ':' or field[4] != ':':
      raise ValueError, "Malformatted SAM optional field '%'" % field
   if field[3] == 'A':
      return field[5]
   elif field[3] == 'i':
      return int( field[5:] )
   elif field[3] == 'f':
      return float( field[5:] )
   elif field[3] == 'Z':
      return field[5:]
   elif field[3] == 'H':
      return int( field[5:], 16 )
   else:
      raise ValueError, "SAM optional field with illegal type letter '%s'" % field[2]

cigar_operation_codes = [ 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
      
      
cdef class SAM_Alignment( AlignmentWithSequenceReversal ):

   """When reading in a SAM file, objects of the class SAM_Alignment
   are returned. In addition to the 'read', 'iv' and 'aligned' fields (see 
   Alignment class), the following fields are provided:
    - aQual: the alignment quality score
    - cigar: a list of CigarOperatio objects, describing the alignment
    - tags: the extra information tags [not yet implemented]
   """
   
   @classmethod
   def from_pysam_AlignedRead( cls, read, samfile ):
      strand = "-" if read.is_reverse else "+"
      if read.tid != -1:
          chrom = samfile.getrname(read.tid)
          iv = GenomicInterval( chrom, read.pos, read.aend, strand )
      else:
          iv = None
      
      seq = SequenceWithQualities( read.seq, read.qname, read.qual )
      a = SAM_Alignment( seq, iv )
      a.cigar = build_cigar_list( [ (cigar_operation_codes[code], length) for (code, length) in read.cigar ] , read.pos, chrom, strand ) if iv != None else []
      a.inferred_insert_size = read.isize
      a.aQual = read.mapq
      a.proper_pair = read.is_proper_pair
      a.not_primary_alignment = read.is_secondary
      a.failed_platform_qc = read.is_qcfail
      a.pcr_or_optical_duplicate = read.is_duplicate
      a.original_sam_line = ""
      if read.is_paired:
         if read.is_proper_pair:
            strand = "-" if read.mate_is_reverse else "+"
            a.mate_start = GenomicPosition( samfile.getrname(read.mrnm), read.mpos, strand )
            a.pe_which = "first" if read.is_read1 else "second" #TODO:check wheter that actually works as expected, what about 'unknown'?
      return a
         
   @classmethod
   def from_SAM_line( cls, line ):
      cdef str qname, flag, rname, pos, mapq, cigar, 
      cdef str mrnm, mpos, isize, seq, qual
      cdef list optional_fields
      cdef int posint, flagint
      cdef str strand
      cdef list cigarlist
            
      fields = line.rstrip().split( "\t" )
      if len( fields ) < 10:
         raise ValueError, "SAM line does not contain at least 11 tab-delimited fields."
      (qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, 
         seq, qual) = fields[ 0:11 ]
      optional_fields = fields[ 11: ]      
      
      if seq.count( "=" ) > 0:
         raise ValueError, "Sequence in SAM file contains '=', which is not supported."
      if seq.count( "." ) > 0:
         raise ValueError, "Sequence in SAM file contains '.', which is not supported."
      flagint = int( flag )
        
      if flagint & 0x0004:     # flag "query sequence is unmapped" 
         iv = None
         cigar = None
         if rname != "*":     # flag "query sequence is unmapped"      
            warnings.warn( "Malformed SAM line: RNAME != '*' although flag bit &0x0004 set" )
      else:
         if rname == "*":
            raise ValueError, "Malformed SAM line: RNAME == '*' although flag bit &0x0004 cleared"
         posint = int( pos ) - 1   # SAM is one-based, but HTSeq is zero-based!
         if flagint & 0x0010:      # flag "strand of the query"
            strand = "-"
         else:
            strand = "+"
         cigarlist = parse_cigar( cigar, posint, rname, strand )
         iv = GenomicInterval( rname, posint, cigarlist[-1].ref_iv.end, strand )   
            
      alnmt = SAM_Alignment( SequenceWithQualities( seq.upper(), qname, qual ), iv )
         
      alnmt.cigar = cigarlist
      alnmt._optional_fields = optional_fields
      alnmt.aQual = int( mapq )
      alnmt.inferred_insert_size = int( isize )
      alnmt.original_sam_line = line
      
      if flagint & 0x0001:         # flag "read is paired in sequencing"
         if flagint & 0x0008:      # flag "mate is unmapped"
            if mrnm != "*":
               warnings.warn( "Malformed SAM line: MRNM != '*' although flag bit &0x0008 set" )
            alnmt.mate_start = None
         else:
            if mrnm == "*":
               raise ValueError, "Malformed SAM line: MRNM == '*' although flag bit &0x0008 cleared"               
            posint = int( mpos ) - 1
            if flagint & 0x0020:   # flag "strand of the mate"
               strand = "-"
            else:
               strand = "+"           
            alnmt.mate_start = GenomicPosition( mrnm, posint, strand )   
            if alnmt.mate_start.chrom == "=":
               alnmt.mate_start.chrom = alnmt.iv.chrom
         if flagint & 0x0040:
            alnmt.pe_which = intern( "first" )
         elif flagint & 0x0080:
            alnmt.pe_which = intern( "second" )
         else:
            alnmt.pe_which = intern( "unknown" )
      else:
         alnmt.mate_start = None
         alnmt.pe_which = intern( "not_paired_end" )
        
      alnmt.proper_pair = flagint & 0x0002 > 0
      alnmt.not_primary_alignment = flagint & 0x0100 > 0
      alnmt.failed_platform_qc = flagint & 0x0200 > 0
      alnmt.pcr_or_optical_duplicate = flagint & 0x0400 > 0

      return alnmt

   @property
   def paired_end( self ):
      return self.pe_which != "not_paired_end"

   @property
   def mate_aligned( self ):
      return self.mate_start is not None
         
   def get_sam_line( self ):
       
      cdef str cigar = ""
      cdef int flag = 0
      cdef GenomicInterval query_start, mate_start
      cdef CigarOperation cop
       
      if self.pe_which != "not_paired_end":
         self.flag |= 0x0001
         if self.proper_pair:
            self.flag |= 0x0002
         if self.pe_which == "first":
            self.flag |= 0x0040
         elif self.pe_which == "second":
            self.flag |= 0x0080
         if self.pe_which != "unknown":
            raise ValueError, "Illegal value in field 'pe_which'"
       
      if self.aligned:
         query_start = self.iv
         if self.iv.strand == "-":
            flag |= 0x0010
      else:
         query_start = GenomicPosition( "*", -1 )
         flag |= 0x0004
          
      if self.mate_start is not None:
         mate_start = self.mate_start
         if self.mate_start.strand == "-":
            flag |= 0x0020
      else:
         mate_start = GenomicPosition( "*", -1 )
         if self.paired_end:
            flag |= 0x0008
          
      if self.proper_pair:
         flag |= 0x0002
      if self.not_primary_alignment:
         flag |= 0x0100
      if self.failed_platform_qc: 
         flag |= 0x0200
      if self.pcr_or_optical_duplicate:
         flag |= 0x0400
          
      if self.cigar is not None:
         for cop in self.cigar:
            cigar += str(cop.size) + cop.type
      else:
         cigar = "*"
       
      return '\t'.join( ( self.read.name, str(flag), query_start.chrom, 
          str(query_start.start+1), str(self.aQual), cigar, mate_start.chrom, 
          str(mate_start.pos+1), str(self.inferred_insert_size), 
           self.read_as_aligned.seq, self.read_as_aligned.qualstr,
           ' '.join( self._optional_fields ) ) )      

   def optional_field( SAM_Alignment self, str tag ):
      cdef str field
      for field in self._optional_fields:
         if field[:2] == tag:
            return _parse_SAM_optional_field_value( field )
      raise KeyError, "SAM optional field tag %s not found" % tag

   def optional_fields( SAM_Alignment self ):
      cdef str field
      return dict( [ ( field[:2], _parse_SAM_optional_field_value( field ) ) 
         for field in self._optional_fields ] )

