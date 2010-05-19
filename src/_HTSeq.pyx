import sys
import math
import re
import csv
import gzip
import urllib
import itertools 
import collections
import cStringIO

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
         str strand ):
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
         if not( intern(strand) is strand_plus or intern(strand) is strand_minus or 
               intern(strand) is strand_nostrand ):
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
      return "<%s object '%s', [%d,%d), strand '%s'>" % \
         ( self.__class__.__name__, self.chrom, self.start, self.end, self.strand )
         
   def __str__( GenomicInterval self ):
         return "%s:[%d,%d)/%s" % \
            ( self.chrom, self.start, self.end, self.strand )

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
   
   
cdef class GenomicArray( object ):
   
   cdef public dict step_vectors
   cdef readonly bool stranded
   cdef readonly str typecode
   cdef public bool auto_add_chroms
   
   def __init__( self, object chroms, bool stranded=True, str typecode='d' ):
      cdef str chrom
      self.step_vectors = {}
      self.stranded = stranded
      self.typecode = typecode
      self.auto_add_chroms = chroms == "auto"
      if self.auto_add_chroms:
         chroms = []
      elif isinstance( chroms, list ):
         chroms = dict( [ ( c, sys.maxint ) for c in chroms ] )
      elif not isinstance( chroms, dict ):
         raise TypeError, "'chroms' must be a list or a dict or 'auto'."
      for chrom in chroms:
         self.add_chrom( chrom, chroms[chrom] )
            
   def __getitem__( self, index ):
      if isinstance( index, GenomicInterval ):
         if self.stranded and index.strand not in ( strand_plus, strand_minus ):
            raise KeyError, "Non-stranded index used for stranded GenomicArray."
         if isinstance( index, GenomicPosition ):
            if self.stranded:
               return self.step_vectors[ index.chrom ][ index.strand ][ index.pos ]
            else:
               return self.step_vectors[ index.chrom ][ index.pos ]
         else:
            if self.stranded:
               return self.step_vectors[ index.chrom ][ index.strand ][ index.start : index.end ]
            else:
               return self.step_vectors[ index.chrom ][ index.start : index.end ]
      else:
         return self.step_vectors[ index ]

   def __setitem__( self, index, value ):
      if isinstance( index, GenomicInterval ):
         if self.stranded and index.strand not in ( strand_plus, strand_minus ):
            raise KeyError, "Non-stranded index used for stranded GenomicArray."
         if self.auto_add_chroms and index.chrom not in self.step_vectors:
            self.add_chrom( index.chrom )
         if isinstance( index, GenomicPosition ):
            if self.stranded:
               self.step_vectors[ index.chrom ][ index.strand ][ index.pos ] = value
            else:
               self.step_vectors[ index.chrom ][ index.pos ] = value
         else:
            if self.stranded:
               self.step_vectors[ index.chrom ][ index.strand ][ index.start : index.end ] = value
            else:
               self.step_vectors[ index.chrom ][ index.start : index.end ] = value
      else:
         raise TypeError, "Only GenomicInterval and GenomicPosition objects " + \
            "are supported for element replacement in GenomicArray objects."
            
   def __richcmp__( self, GenomicArray other, int op ):
      if op == 2:
         return other is not None and \
            self.stranded == other.stranded and self.step_vectors == other.step_vectors
      elif op == 3:
         return other is None or \
            self.stranded != other.stranded or self.step_vectors != other.step_vectors
      else:
         raise NotImplemented
   
   def add_chrom( self, chrom, length = sys.maxint, start_index = 0 ):
      if self.stranded:
         self.step_vectors[ chrom ] = {
            strand_plus:  StepVector.StepVector( length, self.typecode, start_index ),
            strand_minus: StepVector.StepVector( length, self.typecode, start_index ) }
      else:   
         self.step_vectors[ chrom ] = StepVector.StepVector( length, self.typecode, start_index )
   
   
   def add_value( self, value, iv ):
      if self.stranded:
         self.step_vectors[ iv.chrom ][ iv.strand ].add_value( value, iv.start, iv.end )
      else:
         self.step_vectors[ iv.chrom ].add_value( value, iv.start, iv.end )
         
   def get_steps( self, GenomicInterval iv = None, bool values_only=False ):
      if iv is None:
         return _HTSeq_internal.Genomic_array_get_all_steps( self )
      if self.stranded:
         if iv.strand not in ( strand_plus, strand_minus ):
            raise KeyError, "Non-stranded index used for stranded GenomicArray."
         a = self.step_vectors[ iv.chrom ][ iv.strand ].get_steps( iv.start, iv.end, values_only )
      else:
         a = self.step_vectors[ iv.chrom ].get_steps( iv.start, iv.end, values_only )
      if values_only:
         return a
      else:
         if self.stranded:
            return _HTSeq_internal.GenomicArray_get_steps_convert_iv( a, iv.chrom, iv.strand )
         else:
            return _HTSeq_internal.GenomicArray_get_steps_convert_iv( a, iv.chrom, "." )
      
   def __reduce__( self ):
      return ( _GenomicArray_unpickle, ( self.stranded, self.typecode, self.step_vectors ) )
      
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
      for chrom in self.step_vectors:
         if self.stranded:
            sv =  self.step_vectors[ chrom ][ strand ]
         else:
            sv =  self.step_vectors[ chrom ]
         for start, stop, value in sv.get_steps():
            if start == -sys.maxint-1 or stop == sys.maxint:
               continue
            f.write( "%s\t%d\t%d\t%f\n" % (chrom, start, stop, value) )    
      if not hasattr( file_or_filename, "write" ):
         f.close()
    
   def apply( self, func, iv = None ):
      for siv, value in self.get_steps( iv ):
         if self.stranded:
             self.step_vectors[siv.chrom][siv.strand][siv.start:siv.end] = func( value )
         else:
             self.step_vectors[siv.chrom][siv.start:siv.end] = func( value )
    
   
def _GenomicArray_unpickle( stranded, typecode, step_vectors ):
   ga = GenomicArray( {}, stranded, typecode )
   ga.step_vectors = step_vectors
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
   
cdef str _translation_table_for_complementation = _make_translation_table_for_complementation( )

cpdef str reverse_complement( str seq ):
   """Returns the reverse complement of DNA sequence 'seq'. Does not yet
   work with extended IUPAC nucleotide letters or RNA."""

   return seq[ ::-1 ].translate( _translation_table_for_complementation )

base_to_column = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4 }

cdef class Sequence( object ):
   """A Sequence, typically of DNA, with a name.
   """
   
   def __init__( self, str seq, str name="unnamed" ):
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
      return "<%s.%s object '%s' (length %d)>" % ( self.__class__.__module__,
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
         elif b == 'N' or b == 'n':
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

   def __init__( self, str seq, str name, str qualstr, str qualscale="phred" ):
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
      self._qualstr_phred = ''

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
         self._qualstr = ""
         self._qualscale = "none"
         self._qualstr_phred = ""
              
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
            self._qualstr_phred = ' ' * seqlen
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
      for i in xrange( ( len( self ) + 1 ) // 70 + 1):
         fastq_file.write( self.seq[ i*70 : (i+1)*70 ] + "\n" )
      fastq_file.write( "+\n" )
      for i in xrange( ( len( self ) + 1 ) // 70 + 1):
         fastq_file.write( self.qualstr[ i*70 : (i+1)*70 ] + "\n" )

   def get_fastq_str( self, bool convert_to_phred=False ):
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
               sum_mm_qual += qual_array[ seqlen - 1 - i + j ]
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
               sum_mm_qual += qual_array[ seqlen - 1 - i + j ]
               if sum_mm_qual > max_mm_qual:
                  break
         else:
            match = i
      return self[ 0 : seqlen-match ]

   
###########################
##   Alignment
###########################

cdef class Alignment( object ):

   """Alignment abstract base type:

   An alignment object can be defined in different ways but will always
   provide these attributes:
     read:      a SequenceWithQualities object with the read
     aligned:   whether the read is aligned
     iv:        a GenomicInterval object with the alignment position 
   """
   
   def __init__( self ):
      raise NotImplemented, "Alignment is an abstract base class"
   
   def __repr__( self ):
      cdef str s
      if self.paired_end:
         s = "Paired-end read"
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
   
   def __init__( self, type_, size, rfrom, rto, qfrom, qto, chrom, strand ):
      self.type = type_
      self.size = size
      self.ref_iv = GenomicInterval( chrom, rfrom, rto, strand )
      self.query_from = qfrom
      self.query_to = qto
      
   def __repr__( self ):
      return "< %s: %d base(s) %s on ref iv %s, query iv [%d,%d) >" % (
         self.__class__.__name__, self.size, cigar_operation_names[ self.type ],
         str( self.ref_iv ), self.query_from, self.query_to )

cpdef list parse_cigar( str cigar_string, int ref_left = 0, str chrom = "", str strand = "." ):
   cdef list split_cigar, res
   cdef int rpos, qpos, size
   cdef str code
   cigar_codes = re.compile( '([A-Z])' )
   split_cigar = cigar_codes.split( cigar_string )
   if split_cigar[-1] != '' or len(split_cigar) % 2 != 1:
      raise ValueError, "Illegal CIGAR string '%s'" % cigar_string
   rpos = ref_left
   qpos = 0
   res = []
   for i in xrange( len(split_cigar) // 2 ):
      try:
         size = int( split_cigar[2*i] )
      except ValueError:
         raise ValueError, "Illegal CIGAR string '%s'" % cigar_string
      code  = split_cigar[2*i+1]
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
            'N', size, rpos, rpos, qpos, qpos + size, chrom, strand ) )
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
      
cdef class SAM_Alignment( AlignmentWithSequenceReversal ):

   """When reading in a SAM file, objects of the class SAM_Alignment
   are returned. In addition to the 'read', 'iv' and 'aligned' fields (see 
   Alignment class), the following fields are provided:
    - aQual: the alignment quality score
    - cigar: a list of CigarOperatio objects, describing the alignment
    - tags: the extra information tags [not yet implemented]
   """
     
   def __init__( self, line ):
      cdef str qname, flag, rname, pos, mapq, cigar, 
      cdef str mrnm, mpos, isize, seq, qual
      cdef list tags
      cdef int posint, flagint
      cdef str strand
      
      fields = line.rstrip().split( "\t" )
      if len( fields ) < 10:
         raise ValueError, "SAM line does not contain at least 11 tab-delimited fields."
      (qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, 
         seq, qual) = fields[ 0:11 ]
      tags = fields[ 11: ]      
      
      if seq.count( "=" ) > 0:
         raise ValueError, "Sequence in SAM file contains '=', which is not supported."
      if seq.count( "." ) > 0:
         raise ValueError, "Sequence in SAM file contains '.', which is not supported."
      flagint = int( flag )
        
      if flagint & 0x0004:         # flag "query sequence is unmapped"
         if rname != "*":
            raise ValueError, "Malformed SAM line: RNAME != '*' although flag bit 0x0004 set"
         iv = None
         self.cigar = None
      else:
         if rname == "*":
            raise ValueError, "Malformed SAM line: RNAME == '*' although flag bit &0x0004 cleared"
         posint = int( pos ) - 1   # SAM is one-based, but HTSeq is zero-based!
         if flagint & 0x0010:      # flag "strand of the query"
            strand = "-"
         else:
            strand = "+"
         self.cigar = parse_cigar( cigar, posint, rname, strand )
         iv = GenomicInterval( rname, posint, self.cigar[-1].ref_iv.end, strand )   
            
      AlignmentWithSequenceReversal.__init__( self,
         SequenceWithQualities( seq.upper(), qname, qual ), iv )
         
      self._tags = tags
      self.flags = flagint
      self.aQual = int( mapq )
      self.inferred_insert_size = int( isize )
      
      if flagint & 0x0001:         # flag "read is paired in sequencing"
         if flagint & 0x0008:      # flag "mate is unmapped"
            if mrnm != "*":
               raise ValueError, "Malformed SAM line: MRNM != '*' although flag bit &0x0008 set"
            self.mate_start = None
         else:
            if mrnm == "*":
               raise ValueError, "Malformed SAM line: MRNM == '*' although flag bit &0x0008 cleared"
            posint = int( mpos ) - 1
            if flagint & 0x0020:   # flag "strand of the mate"
               strand = "-"
            else:
               strand = "+"           
            self.mate_start = GenomicPosition( mrnm, posint, strand )   
            if self.mate_start.chrom == "=":
               self.mate_start.chrom = self.iv.chrom
         if flagint & 0x0040:
            self.pe_which = intern( "first" )
         elif flagint & 0x0080:
            self.pe_which = intern( "second" )
         else:
            self.pe_which = intern( "unknown" )
      else:
        self.mate_start = None
        self.pe_which = intern( "not_paired_end" )
        
   @property
   def paired_end( self ):
      return bool( self.flags & 0x0001 )

   @property
   def proper_pair( self ):
      return bool( self.flags & 0x0002 )
         
   @property
   def mate_aligned( self ):
      if not ( self.flags & 0x0001 ):
         raise ValueError, "Not a paired-end read"
      return not bool( self.flags & 0x0008 )

   @property
   def passed_filter( self ):
      return not bool( self.flags & 0x0200 )
      
   @property
   def primary( self ):
      return not bool( self.flags & 0x0100 )

   @property
   def duplicate( self ):
      return bool( self.flags & 0x0400 )
         
         

