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

cdef class Strand:

   """The class Strand is an encapsulation of an enumeration of three possible
   values, namely the plus strand ("+"), the minus strand ("-") and no strand
   ("."). These values should be used when specifying a genomic interval or
   position. Use "." when specifying a strand is inappropriate.
   """

   def __init__( Strand self, strand ):
      """To instatiate a Strand object, pass one of the strings '+', '-', or '.'
      to the constructor.
      """
      cdef char * strand_cstr
      if isinstance( strand, Strand ):
         self.se = (<Strand> strand).se
      else:
         if len( strand ) != 1:
            raise TypeError, "Invalid strand specification"   
         strand_cstr = <char*?> strand
         if strand_cstr[0] == ".":
            self.se = strand_nostrand
         elif strand_cstr[0] == "+":
            self.se = strand_plus
         elif strand_cstr[0] == "-":
            self.se = strand_minus
         else:
            raise TypeError, "Invalid strand specification"   
   
   def __str__( Strand self ):
      if self.se == strand_nostrand:
         return '.'
      elif self.se == strand_plus:
         return '+'
      elif self.se == strand_minus:
         return '-'
      
   def __repr__( Strand self ):
      return "<%s '%s'>" % (self.__class__.__name__, str(self) )
      
   def __hash__( Strand self ):
      return hash( self.se )                        
            
   def __richcmp__( Strand self, other, int op ):
      """Strands can be compared with the operators "==" and "!=".
      """
      if op == 2:
         if isinstance( other, Strand ):
            return self.se == (<Strand> other).se
         else:
            try:
               return self == Strand( str( other ) )
            except TypeError:
               return False
      elif op == 3:
         if isinstance( other, Strand ):
            return self.se != (<Strand> other).se
         else:
            try:
               return self != Strand( str( other ) )
            except TypeError:
               return True
      elif op == 0:
         return self.se < (<Strand> other).se
      elif op == 1:
         return self.se <= (<Strand> other).se
      elif op == 4:
         return self.se > (<Strand> other).se
      elif op == 5:
         return self.se >= (<Strand> other).se
      else:
         raise NotImplementedError
      
   def __reduce__( Strand self ):
      return( Strand, ( str( self ), ) )

   def __copy__( self ):
      constr, args = self.__reduce__()
      return constr( *args )         
   

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
      genome: A Genome object, specifying the reference (may be None)
      length: The length of the interval, i.e., end - start
      start_d: The "directional start" position. This is the position of the
        first base of the interval, taking the strand into account. Hence, 
        this is the same as 'start' except when strand == '-', in which 
        case it is end-1.
      end_d: The "directional end": Usually, the same as 'end', but for 
        strand=='-1', it is start+1.
   """
   
   def __init__( GenomicInterval self, str chrom, long start, long end, 
         object strand, object genome=None ):
      """See the class docstring for the meaning of the slots. Note that 
      there is also a factory function, 'from_directional', to be used if
      you wish to specify start_d and length.
      """
      self.chrom = intern( chrom )
      self.start = start
      self.end = end
      if isinstance( strand, Strand ):
         self.strand = <Strand> strand
      else:
         self.strand = Strand( strand )
      self.genome = genome
      if self.start > self.end:
         raise ValueError, "start is larger than end"
            
   def __reduce__( GenomicInterval self ):
      return GenomicInterval, ( self.chrom, self.start, self.end, 
         self.strand, self.genome )
         
   def __copy__( self ):
      constr, args = self.__reduce__()
      return constr( *args )   
   
   def __repr__( GenomicInterval self ):
      if self.genome:
         return "<%s on '%s' (%s genome), [%d,%d), strand '%s'>" % \
            ( self.__class__.__name__, self.chrom, str(self.genome), self.start, self.end, self.strand )
      else:
         return "<%s object '%s', [%d,%d), strand '%s'>" % \
            ( self.__class__.__name__, self.chrom, self.start, self.end, self.strand )
         
   def __str__( GenomicInterval self ):
         return "%s:[%d,%d)/%s" % \
            ( self.chrom, self.start, self.end, self.strand )

   property strand:
   
      def __get__( GenomicInterval self ):
         return self._strand
         
      def __set__( GenomicInterval self, object newStrand ):
         self._strand = Strand( newStrand )
         
   property length:

      """The length is calculated as end - start. If you set the length, 
      'start_d' will be preserved, i.e., 'end' is changed, unless the strand
      is '-', in which case 'start' is changed."""   
            
      def __get__( GenomicInterval self ):
         return self.end - self.start

      def __set__( GenomicInterval self, long newLength ):
         if self._strand.se != strand_minus:
            self.end = self.start + newLength
         else:
            self.start = self.end - newLength

   property start_d:
      """See the class docstring for the meaning of the 'directional start'.
      Note that if you set 'start_d', both the start and the end are changed, 
      such teh interval gets the requested new directional start and its
      length stays unchanged."""

      def __get__( GenomicInterval self ):
         if self._strand.se != strand_minus:
            return self.start
         else:
            return self.end - 1

      def __set__( GenomicInterval self, long newStartd ):
         if self._strand.se != strand_minus:
            self.end = newStartd + self.length
            self.start = newStartd
         else:
            self.start = newStartd + 1 - self.length
            self.end = newStartd + 1
         
   property end_d:
   
      def __get__( GenomicInterval self ):
         if self._strand.se != strand_minus:
            return self.end
         else:
            return self.start + 1
         
   def __richcmp__( GenomicInterval self, GenomicInterval other, int op ):
      if op == 2:  # ==
         if other == None:
            return False
         if self.genome != None and other.genome != None \
               and self.genome != other.genome:
            return False
         return self.strand == other._get_strand() and \
            self.start == other.start and self.end == other.end
      elif op == 3:  # !=
         return not ( self == other )
      else:
         raise NotImplementedError
         
   def __hash__( GenomicInterval self ):
      return hash( ( self.genome, self.chrom, self.start, self.end ) )
           
   cdef Strand _get_strand( GenomicInterval self ):
      return self._strand

   cpdef is_contained_in( GenomicInterval self, GenomicInterval iv ):
      """Returns a boolean value indicating whether the 'self' interval 
      is fully within the 'iv' interval.
      
      This is deemed the case if
        - both are on the same genome, or at least one of the has
           'None' as genome field, and
        - both are on the same chromosome, and    
        - both are on the same strand, or at least one of them is
           not stranded (i.e., has strand == '.'), and
        - self.start >= iv.start, and
        - self.end <= iv.end
      """
      if iv == None:
         return False
      if self.genome != None and iv.genome != None and not (self.genome is iv.genome):
         return False 
      if self.chrom != iv.chrom:
         return False
      if self._strand.se != strand_nostrand and iv._strand.se != strand_nostrand and \
            self._strand.se != iv._strand.se:
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
        - both are on the same genome, or at least one of the has
           'None' as genome field, and
        - both are on the same chromosome, and    
        - both are on the same strand, or at least one of them is
           not stranded (i.e., has strand == '.'), and
        - the actual intervals overlap
      """
      if iv == None:
         return False
      if self.genome != None and iv.genome != None and not (self.genome is iv.genome):
         return False 
      if self.chrom != iv.chrom:
         return False
      if self._strand.se != strand_nostrand and iv._strand.se != strand_nostrand and \
            self._strand.se != iv._strand.se:
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

   def xranged( GenomicInterval self, long int step = 1 ):
      """Generate an iterator over the GenomicPositions covered by the interval.
      running from start_d to end_d.
      """
      return _HTSeq_internal.GenomicInterval_xranged( self, step )
      
   cpdef extend_to_include( GenomicInterval self, GenomicInterval iv ):
      """Extend the interval such that it includes iv."""
      if iv is None:
         raise TypeError, "Cannot extend an interval to include None."
      if self.genome != None and iv.genome != None and not (self.genome is iv.genome):
         raise ValueError, "Cannot extend an interval to include an interval on another genome."
      if self.chrom != iv.chrom:
         raise ValueError, "Cannot extend an interval to include an interval on another chromosome."
      if self._strand.se != strand_nostrand and iv._strand.se != strand_nostrand and \
            self._strand.se != iv._strand.se:
         raise ValueError, "Cannot extend an interval to include an interval on another strand."
      self.start = min( self.start, iv.start )
      self.end = max( self.end, iv.end )


def GenomicInterval_from_directional( str chrom, long int start_d, long int length, object strand=".", genome=None ):
   cdef Strand strand2 = Strand( strand )
   if strand2.se != strand_minus:
      return GenomicInterval( chrom, start_d, start_d+length, strand2, genome )
   else:
      return GenomicInterval( chrom, start_d-length+1, start_d+1, strand2, genome )


cdef class GenomicPosition( GenomicInterval ):

   """A GenomicPosition specifies the position of a nucleotide or
   base pair on a reference genome.
   
   It has the following slots:
      chrom: The name of a sequence (i.e., chromosome, contig, or 
         the like). If 'genome' is given, it should know 'seqname'
      pos: The position on the sequence specified by seqname.
         The position should always be given as 0-based value!
      strand: The strand, as a single character, '+' or '-'. '.' indicates
         that the strand is irrelavant.
      genome: A Genome object, specifying the reference (may be None)

   The GenomicPosition class is derived from GenomicInterval. Hence,
   a GenomicPosition is always a GenomicInterval of length 1. Do not tinker
   with the exposed GenomeInterval slots.
   """

   def __init__( self, str chrom, long int pos, strand='.', genome=None ):
      GenomicInterval.__init__( self, chrom, pos, pos+1, strand, genome )
      
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
      if self.genome:
         return "<%s on '%s':%d (%s genome), strand '%s'>" % \
            ( self.__class__, self.chrom, self.pos, self.genome, self.strand )
      else:
         return "<%s object '%s':%d, strand '%s'>" % \
            ( self.__class__.__name__, self.chrom, self.pos, self.strand )               
            
   def __str__( self ):
      return "%s:%d/%s" % ( self.chrom, self.pos, self.strand )            

   def __reduce__( GenomicPosition self ):
      return GenomicPosition, ( self.chrom, self.pos, self.strand, self.genome )
   
   
cdef class GenomicArray( object ):
   
   cdef public dict step_vectors
   cdef readonly bool stranded
   cdef readonly str typecode
   
   def __init__( self, dict chrom_lengths, bool stranded=True, str typecode='d' ):
      cdef str chrom
      self.step_vectors = {}
      self.stranded = stranded
      self.typecode = typecode
      for chrom in chrom_lengths:
         if self.stranded:
            self.step_vectors[ chrom ] = {
               Strand("+"): StepVector.StepVector( chrom_lengths[chrom], typecode ),
               Strand("-"): StepVector.StepVector( chrom_lengths[chrom], typecode ) }
         else:   
            self.step_vectors[ chrom ] = StepVector.StepVector(
               chrom_lengths[chrom], typecode )
            
   def __getitem__( self, index ):
      if isinstance( index, GenomicInterval ):
         if self.stranded and index.strand not in ( Strand("+"), Strand("-") ):
            raise KeyError, "Non-stranded index used for stranded GenomicArray."
         if isinstance( index, GenomicPosition ):
            if self.stranded:
               return self.step_vectors[ index.chrom ][ index.strand ][ index.pos ]
            else:
               return self.step_vectors[ index.chrom ][ index.pos ]
         else:
            res = GenomicArray( {}, self.stranded, self.typecode )
            if self.stranded:
               res.step_vectors[index.chrom ] = self.step_vectors[ index.chrom ][ index.strand ][ index.start : index.end ]
            else:
               res.step_vectors[index.chrom ] = self.step_vectors[ index.chrom ][ index.start : index.end ]
            return res
      else:
         return self.step_vectors[ index ]

   def __setitem__( self, index, value ):
      if isinstance( index, GenomicInterval ):
         if self.stranded and index.strand not in ( Strand("+"), Strand("-") ):
            raise KeyError, "Non-stranded index used for stranded GenomicArray."
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
   
   def add_value( self, value, iv ):
      if self.stranded:
         self.step_vectors[ iv.chrom ][ iv.strand ].add_value( value, iv.start, iv.end )
      else:
         self.step_vectors[ iv.chrom ].add_value( value, iv.start, iv.end )
         
   def get_steps( self, GenomicInterval iv = None, bool values_only=False ):
      if iv is None:
         return _HTSeq_internal.Genomic_array_get_all_steps( self )
      if self.stranded:
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
      if self.stranded and strand not in ( "+", "-" ):
         raise ValueError, "Strand must be specified for stranded GenomicArray."
      if hasattr( file_or_filename, "write" ):
         f = file_or_filename
      else:
         f = open( file_or_filename, "w" )
      if track_options == "":
         f.write( "track type=bedGraph\n" )
      else:
         f.write( "track type=bedGraph %s\n" % track_options )
      strand = Strand( strand )
      for chrom in self.step_vectors:
         if self.stranded:
            sv =  self.step_vectors[ chrom ][ strand ]
         else:
            sv =  self.step_vectors[ chrom ]
         for start, stop, value in sv.get_steps():
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


cdef class Sequence( object ):
   """A Sequence, typically of DNA, with a name.
   """
   
   cdef public str seq
   cdef public str name
   cdef public str descr
   
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
      for i in xrange( ( len( self.seq ) + 1 ) // 70):
         fasta_file.write( self.seq[ i*70 : (i+1)*70 ] + "\n" )

cdef class SequenceWithQualities( Sequence ):
   """A Sequence with base-call quality scores.
   It now has property  'qual', an integer NumPy array of Sanger/Phred 
   quality scores of the  base calls.
   """

   cdef readonly str _qualstr
   cdef readonly str _qualscale
   cdef readonly object _qualarr

   def __init__( self, str seq, str name, str qualstr, str qualscale="phred" ):
      """ Construct a SequenceWithQuality object.
      
        seq       - The actual sequence.
        name      - The sequence name or ID
        qualstr   - The quality string. Must have the same length as seq
        qualscale - The encoding scale of the quality string. Must be one of
                      "phred", "solexa-old", "solexa-1.3" )
      """
      Sequence.__init__( self, seq, name )
      self._qualstr = qualstr
      self._qualscale = qualscale
      self._qualarr = None
      
   @property
   def qual( self ):
      if self._qualarr is not None:
         return self._qualarr
      else:
         assert len( self.seq ) == len( self._qualstr )
         if self._qualscale == "phred":
            self._qualarr = numpy.array( [ ord( c ) - 33 for c in self._qualstr ] )
         else:
            self._qualarr = numpy.array( [ ord( c ) - 64 for c in self._qualstr ], 'd' )
            if self._qualscale == "solexa-old":
               self._qualarr = 10 * numpy.log10(1 + 10 ** ( self._quastrl / 10.0 ) )
            else:
               if self._qualscale != "solexa-1.3":
                  raise ValueError, "Illegal quality scale '%s'." % self._qualscale
         return self._qualarr
              
   def __repr__( self ):
      return "<%s object '%s'>" % ( self.__class__.__name__, self.name )

   def __getitem__( self, item ):
      if self.name.endswith( "[part]" ):
         new_name = self.name
      else:
         new_name = self.name + "[part]"
      return SequenceWithQualities( 
         self.seq[ item ], new_name, self.qualstr[ item ], 
            self.qualstrIsSolexaScale )

   def write_to_fastq_file( self, fastq_file, bool convert_to_phred=False ):
      if convert_to_phred:
         raise NotImplemented
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

base_to_row = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4 }

def base_counts_by_position( sequences, seq_length=None ):
   """This function takes an iterable of Sequence objects and counts for each 
   position in the sequence the number of 'A's, 'C's, 'T's, 'G's and 'N's. It
   returns a two-dimensional integer numpy array with 'seq_length' rows and
   5 columns. The 5 columns correspond to the letters A, C, G, T, N, in this
   order. The 'base_to_row' dictionary allows to address the columns
   conveniently.
   
   If 'seq_length' is not provided, the length of the first sequence is used.
   If your sequences do not all have the same length, you should provide
   'seq_length'.
   
   
   Example:
   
   >>> reads = HTSeq.FastqReader( fastq_file )
   >>> counts = HTSeq.base_counts_by_position( reads )
   >>> print "There are", counts[ 17, 1 ], "reads with a 'C' at position 17."
   
   Instead of 'counts[ 17, 1 ]', you can write, more verbosely,
   'counts[ 17, HTSeq.base_to_row["C"] ]'.
   
   
   Notes:    
   
   - Make sure that all elements of the iterable are of the type Sequence
   (or a sub-class, e.g., SequenceWithQualities). For efficiency reasons, only
   the type of the first object in the sequence is checked, and a wrong type of
   a later object may crash the program.
   
   - Small letters are counted as well but letters other than A, C, G, T, N,
   a, c, g, t, n cause a ValueError.
   """
   cdef int seq_length_int
   cdef Sequence seq
   cdef numpy.ndarray[ numpy.int_t, ndim=2 ] counts
   first_seq, sequences = _HTSeq_internal.peek( sequences )
   if not isinstance( first_seq, Sequence ):
      raise TypeError, "First argument must be an iterable providing Sequence objects."
   if seq_length is None:
      seq_length_int = len( first_seq )
   else:
      seq_length_int = seq_length
   counts = numpy.zeros( ( seq_length_int, 5 ), dtype=numpy.int )
   for seq in sequences:
      if len( seq ) > seq_length_int:
         raise ValueError, "Sequence %s is longer than %d (the specified length)" % ( 
            str(seq), seq_length_int )
      for i in xrange( len( seq ) ):
         if seq.seq[i] == 'A' or seq.seq[i] == 'a':
            counts[ i, 0 ] += 1
         elif seq.seq[i] == 'C' or seq.seq[i] == 'c':
            counts[ i, 1 ] += 1
         elif seq.seq[i] == 'G' or seq.seq[i] == 'g':
            counts[ i, 2 ] += 1
         elif seq.seq[i] == 'T' or seq.seq[i] == 't':
            counts[ i, 3 ] += 1
         elif seq.seq[i] == 'N' or seq.seq[i] == 'n':
            counts[ i, 4 ] += 1
         else:
            raise ValueError, "Illegal letter in sequence " + str( seq )
   return counts
   


###########################
##   Alignment
###########################

cdef class Alignment( object ):

   """Alignment abstract base type:

   An alignment object can be defined in different ways but will always
   provide these attributes:
     read: a SequenceWithQualities object with the aligned read
     iv:   a GenomicInterval object with the alignment position
   """
   
   def __init__( self ):
      raise NotImplemented, "Alignment is an abstract base class"
   
   def __repr__( self ):
      return "<%s object: Read '%s' aligned to %s>" % (
         self.__class__.__name__, self.read.name, str(self.iv) )

cdef class AlignmentWithSequenceReversal( Alignment ):

   """Many aligners report the read's sequence in reverse-complemented form
   when it was mapped to the reverse strand. For such alignments, a 
   daughter class of this one should be used.
   
   Then, the read is stored as aligned in the 'read_as_aligned' field,
   and get reverse-complemented back to the sequenced form when the 'read'
   attribute is sequenced.
   """

   cdef public SequenceWithQualities read_as_aligned
   cdef public SequenceWithQualities _read_as_sequenced
   cdef public GenomicInterval iv

   def __init__( self, SequenceWithQualities read_as_aligned, GenomicInterval iv ):
      self.read_as_aligned = read_as_aligned      
      self._read_as_sequenced = None
      self.iv = iv
      
   property read:
      def __get__( self ):
         if self._read_as_sequenced is None:
            if self.iv.strand != "-":
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
      self.ref_iv = GenomicInterval( chrom, rfrom, rto, Strand( strand ) )
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
   are returned. In addition to the 'read' and 'iv' fields (see Alignment
   class), the fields [...] are provided. These 
   contain the content of the respective columns of the SAM file.
   """
   
   cdef public list cigar
   
   def __init__( self, line ):
      cdef str qname, flag, rname, pos, mapq, cigar, 
      cdef str mrnm, mpos, isize, seq, qual, tags
      cdef int posint, flagint
      cdef str strand
      
      (qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, 
         seq, qual, tags) = line.split( "\t", 11 )
      
      if seq.count( "=" ) > 0:
         raise ValueError, "Sequence in SAM file contains '=', which is not supported."
      if seq.count( "." ) > 0:
         raise ValueError, "Sequence in SAM file contains '.', which is not supported."
      if flagint & 0x0001:
         raise ValueError, "Paired-end data encountered; not yet supported."    
        
      posint = int( pos ) - 1   # SAM if one-based, but HTSeq is zero-based!
      flagint = int( flag )
      if flagint & 0x0010:
         strand = "-"
      else:
         strand = "+"
      self.cigar = parse_cigar( cigar, posint, rname, strand )
    
      AlignmentWithSequenceReversal.__init__( self,
         SequenceWithQualities( seq.upper(), qname, qual ),
         GenomicInterval( rname, posint, self.cigar[-1].ref_iv.end, strand ) )
      
         
