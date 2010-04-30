"""HTSeq is a package to process high-throughput sequencing data.

See http://www-huber.embl.de/users/anders/HTSeq for documentation.
"""

import itertools, warnings

from _HTSeq import *

from _HTSeq_internal import peek

from _version import __version__


#########################
## Utils
#########################


class FileOrSequence( object ):
   """ The construcutor takes one argument, which may either be a string,
   which is interpreted as a file name (possibly with path), or a
   connection, by which we mean a text file opened for reading, or 
   any other object that can provide an iterator over strings 
   (lines of the file).
   
   The advantage of passing a file name instead of an already opened file
   is that if an iterator is requested several times, the file will be
   re-opened each time. If the file is already open, its lines can be read
   only once, and then, the iterator stays exhausted.      
   
   Furthermore, if a file name is passed that end in ".gz" or ".gzip"
   (case insensitive), it is transparently gunzipped.
   """
   
   def __init__( self, filename_or_sequence ):      
      self.fos = filename_or_sequence
      
   def __iter__( self ):
      if isinstance( self.fos, str ):
         if self.fos.lower().endswith( ( ".gz" , ".gzip" ) ):
            lines = gzip.open( self.fos )
         else:
            lines = open( self.fos )
      else:
         lines = self.fos
      for line in lines:
         yield line
      if isinstance( self.fos, str ):
         lines.close()
         
   def __repr__( self ):
      if isinstance( self.fos, str ):
         return "<%s object, connected to file name '%s'>" % (
            self.__class__.__name__, self.fos )
      else:
         return "<%s object, connected to %s >" % (
            self.__class__.__name__, repr( self.fos ) )   

class GeneratedSequence( object ):
   """GeneratedSequence is a little helper class to get "respawnable"
   iterator generators. Pass its constructor a function that generates
   an iterator, and its arguments, and it will call this functionion 
   whenever an iterator is requested.
   """
   
   def __init__( self, generator, *args ):
      self.generator = generator
      self.args = args
      
   def __iter__( self ):
      return self.generator( *self.args )

_re_attr_main = re.compile( "\s*(\w+)[\s=]+(.*)" )


#########################
## Features
#########################

class GenomicFeature( object ):
   
   """A genomic feature, i.e., an interval on a genome with metadata.
   
   At minimum, the following information should be provided by slots:
   
     name: a string identifying the feature (e.g., a gene symbol)
     type: a string giving the feature type (e.g., "gene", "exon")
     iv: a GenomicInterval object specifying the feature locus     
   """
   
   def __init__( self, name, type_, interval ):
      self.name = name
      self.type = intern( type_ )
      self.iv = interval
      
   def __repr__( self ):
      return "<%s: %s '%s' at %s: %d -> %d (strand '%s')>" % \
         ( self.__class__.__name__, self.type, self.name, 
         self.iv.chrom, self.iv.start_d, self.iv.end_d, self.iv.strand )
         
   def __eq__( self, other ):
      if not isinstance( other, GenomicFeature ):
         return False
      return self.name == other.name and self.type == other.type and \
         self.iv == other.iv
         
   def __neq__( self, other ):
      if not isinstance( other, GenomicFeature ):
         return True
      return not self.__eq__( other )
   
   def get_gff_line( self, with_equal_sign=False ):
      try:
         source = self.source
      except AttributeError:
         source = "."
      try:
         source = self.score
      except AttributeError:
         score = "."
      try:
         frame = self.frame
      except AttributeError:
         frame = "."
      try:
         attr = self.attr
      except AttributeError:
         attr = { 'ID': self.name }
      if with_equal_sign:
         sep = "="
      else:
         sep = " "
      attr_str = '; '.join( [ '%s%s\"%s\"' % ( ak, sep, urllib.quote( attr[ak] ) ) for ak in attr ] )
      return "\t".join( str(a) for a in ( self.iv.chrom, source, 
         self.type, self.iv.start, self.iv.end, score, 
         self.iv.strand, frame, attr_str ) ) + "\n"
         


def parse_GFF_attribute_string( attrStr, extra_return_first_value=False ):
   """Parses a GFF attribute string and returns it as a dictionary.
   
   If 'extra_return_first_value' is set, a pair is returned: the dictionary
   and the value of the first attribute. This might be useful if this is the ID.
   """
   if attrStr.endswith( "\n" ):
      attrStr = attrStr[:-1]
   d = {}
   for (i, attr) in itertools.izip( itertools.count(), attrStr.split( ";" ) ):
      if attr == "":
         continue
      if attr.count( '"' ) not in ( 0, 2 ):
         raise ValueError, "The attribute string seems to contain mismatched quotes."
      mo = _re_attr_main.match( attr )
      if not mo:
         raise ValueError, "Failure parsing GFF attribute line"
      val = mo.group(2)
      if val.startswith( '"' ) and val.endswith( '"' ):
         val = val[1:-1]
      val = urllib.unquote( val )
      d[ intern(mo.group(1)) ] = intern(val)
      if extra_return_first_value and i == 0:
         first_val = val         
   if extra_return_first_value:
      return ( d, first_val )
   else:
      return d


class GFF_Reader( FileOrSequence ):

   """Parser a GFF file
   
   Pass the constructor either a file name or an iterator of lines of a 
   GFF files. If a file name is specified, it may refer to a gzip compressed
   file.
   
   Iterating over the object then yields GenomicFeature objects.
   """
   
   def __iter__( self ):
      for line in FileOrSequence.__iter__( self ):
         if line.startswith( '#' ):
            continue
         ( seqname, source, feature, start, end, score, 
            strand, frame, attributeStr ) = line.split( "\t", 8 )   
         ( attr, name ) = parse_GFF_attribute_string( attributeStr, True )
         f = GenomicFeature( name, feature,
             GenomicInterval( seqname, int(start)-1, int(end), strand ) )
         if score != ".":
            score = int( score )
         if frame != ".":
            frame = int( frame )
         f.source = source
         f.score = score
         f.frame = frame
         f.attr = attr
         yield f

def make_feature_dict( feature_sequence ):
   """A feature dict is a convenient way to organize a sequence of Feature 
   object (which you have got, e.g., from parse_GFF).
   
   The function returns a dict with all the feature types as keys. Each value
   of this dict is again a dict, now of feature names. The values of this dict
   is a list of feature.
   
   An example makes this clear. Let's say you load the C. elegans GTF file
   from Ensemble and make a feature dict:
   
   >>> worm_features_dict = HTSeq.make_feature_dict( HTSeq.parse_GFF( 
   ...     "test_data/Caenorhabditis_elegans.WS200.55.gtf.gz" ) )

   (This command may take a few minutes to deal with the 430,000 features
   in the GTF file. Note that you may need a lot of RAM if you have millions
   of features.)

   Then, you can simply access, say, exon 0 of gene "F08E10.4" as follows:
   >>> worm_features_dict[ 'exon' ][ 'F08E10.4' ][ 0 ]
   <GenomicFeature: exon 'F08E10.4' at V: 17479353 -> 17479001 (strand '-')>
   """
   
   res = {}
   for f in feature_sequence:
      if f.type not in res:
         res[ f.type ] = {}
      res_ftype = res[ f.type ]
      if f.name not in res_ftype:
         res_ftype[ f.name ] = [ f ]
      else:
         res_ftype[ f.name ].append( f )
   return res



#########################
## GenomicArray
#########################

def read_chrom_lens( filename, delimiter="\t" ):
   return dict( ( ( chrom, int(len) ) 
      for chrom, len in csv.reader( open(filename), delimiter=delimiter ) ) )
            

#########################
## Sequence readers
#########################

_re_fasta_header_line = re.compile( r'>\s*(\S+)\s*(.*)' )      

class FastaReader( FileOrSequence ):
   """A Fasta_Reader is associated with a FASTA file or an open connection
   to a file-like object with content in FASTA format.
   It can generate an iterator over the sequences.
   """
   
   def __iter__( self ):
      seq = None
      for line in FileOrSequence.__iter__( self ):
         if line.startswith( ">" ):
            if seq:
               s = Sequence( seq, name )
               s.descr = descr
               yield s
            mo = _re_fasta_header_line.match( line )
            name = mo.group(1)
            descr = mo.group(2)
            seq = ""
         else: 
            assert seq is not None, "FASTA file does not start with '>'."
            seq += line[:-1]
      if seq:
         s = Sequence( seq, name )
         s.descr = descr
         yield s
         
   def get_sequence_lengths( self ):   
      seqname = None
      seqlengths = {}
      for line in FileOrSequence.__iter__( self ):
         if line.startswith( ">" ):
            if seqname is not None:
               seqlengths[ seqname ] = length
               print seqname, length
            mo = _re_fasta_header_line.match( line )
            seqname = mo.group(1)
            length = 0
         else: 
            assert seqname is not None, "FASTA file does not start with '>'."
            length += len( line.rstrip() )
      if seqname is not None:
         seqlengths[ seqname ] = length
      return seqlengths
            
class FastqReader( FileOrSequence ):
   """A Fastq object is associated with a FASTQ self.file. When an iterator
   is requested from the object, the FASTQ file is read.
   
   qual_scale is one of "phred", "solexa", "solexa-old".
   """

   def __init__( self, file_, qual_scale = "phred" ):
      FileOrSequence.__init__( self, file_ )
      self.qual_scale = qual_scale
      if qual_scale not in ( "phred", "solexa", "solexa-old" ):
         raise ValueError, "Illegal quality scale."
      
   def __iter__( self ):
      fin = FileOrSequence.__iter__( self )
      while True:
         id1  = fin.next()
         seq  = fin.next()
         id2  = fin.next()
         qual = fin.next()
         if qual == "":
            if id1 != "":
               warnings.warn( "Number of lines in FASTQ file is not "
                  "a multiple of 4. Discarding the last, "
                  "incomplete record" )
            break
            
         if not qual.endswith( "\n" ):
            qual += "\n"
         if not id1.startswith( "@" ):
            raise ValueError( "Primary ID line in FASTQ file does"
               "not start with '@'. Either this is not FASTQ data or the parser got out of sync." )
         if not id2.startswith( "+" ):
            raise ValueError( "Secondary ID line in FASTQ file does"
               "not start with '+'. Maybe got out of sync." )
         if len( id2 ) > 2 and id1[1:] != id2[1:]:
            raise ValueError( "Primary and secondary ID line in FASTQ"
               "disagree." )
               
         yield SequenceWithQualities( seq[:-1], id1[1:-1], qual[:-1], 
            self.qual_scale )

class BowtieReader( FileOrSequence ):
   """A BowtieFile object is associated with a Bowtie output file that 
   contains short read alignments. It can generate an iterator of Alignment
   objects."""

   def __iter__( self ):
      for line in FileOrSequence.__iter__( self ):
         try:
            algnt = BowtieAlignment( line )
         except ValueError:
            if line.startswith( "Reported " ):
               continue
            warnings.warn( "BowtieReader: Ignoring the following line, which could not be parsed:\n%s\n" % line,
               RuntimeWarning )
         yield algnt
   
def bundle_multiple_alignments( sequence_of_alignments ):      
   """Some alignment programs, e.g., Bowtie, can output multiple alignments,
   i.e., the same read is reported consecutively with different alignments.
   This function takes an iterator over alignments and bundles consecutive
   alignments regarding the same read to a list of Alignment objects and
   returns an iterator over these.
   """
   alignment_iter = iter( sequence_of_alignments )
   algnt = alignment_iter.next()
   ma = [ algnt ]
   for algnt in alignment_iter:
      if algnt.read.name != ma[0].read.name:
         yield ma
         ma = [ algnt ]
      else:
         ma.append( algnt )
   yield ma   
   
  
class SolexaExportAlignment( Alignment ):
   """Iterating over SolexaExportReader objects will yield SoelxaExportRecord
   objects. These have four fields:
      read          - a SequenceWithQualities object 
      aligned       - a boolean, indicating whether the object was aligned
      iv            - a GenomicInterval giving the alignment (or None, if not aligned)
      passed_filter - a boolean, indicating whether the object passed the filter
      nomatch_code  - a code indicating why no match was found (or None, if the
        read was aligned)
        
   As long as 'aligned' is True, a SolexaExportRecord can be treated as an 
   Alignment object.
   """
   
   def __init__( self ):
      # Data is filled in by SolexaExportRecord
      pass
   
   def __repr__( self ):
      if self.aligned:
         return "< %s object: Read '%s', aligned to %s >" % (
            self.__class__.__name__, self.read.name, self.iv )
      else:
         return "< %s object: Non-aligned read '%s' >" % ( 
            self.__class__.__name__, self.read.name )

class SolexaExportReader( FileOrSequence ):
   """Parser for *_export.txt files from the SolexaPipeline software.
   
   Iterating over a SolexaExportReader yields SolexaExportRecord objects.
   """

   def __init__( self, filename_or_sequence, solexa_old = False ):
      FileOrSequence.__init__( self, filename_or_sequence)
      if solexa_old:
         self.qualscale = "solexa-old"
      else:
         self.qualscale = "solexa-1.3"

   @classmethod
   def parse_line_bare( dummy, line ):
      if line[-1] == "\n":
         line = line[:-1]
      res = {}
      ( res['machine'], res['run_number'], res['lane'], res['tile'], res['x_coord'], 
         res['y_coord'], res['index_string'], res['read_nbr'], res['read_seq'], 
         res['qual_str'], res['chrom'], res['contig'], res['pos'], res['strand'], 
         res['match_descr'], res['single_read_algnt_score'], 
         res['paired_read_algnt_score'], res['partner_chrom'], res['partner_contig'], 
         res['partner_offset'], res['partner_strand'], res['passed_filtering'] ) \
         = line.split( "\t" )
      return res

   def __iter__( self ):   
      for line in FileOrSequence.__iter__( self ):
         record = SolexaExportAlignment()
         fields = SolexaExportReader.parse_line_bare( line )
         if fields['read_nbr'] != "1":
            warnings.warn( "Paired-end read encountered. PE is so far supported only for " +
               "SAM files, not yet for SolexaExport. All PE-related fields are ignored. " )
         record.read = SequenceWithQualities( 
            fields['read_seq'], 
            "%s:%s:%s:%s:%s#0" % (fields['machine'], fields['lane'], fields['tile'], 
               fields['x_coord'], fields['y_coord'] ),
            fields['qual_str'], self.qualscale )
         if fields['passed_filtering'] == 'Y':
            record.passed_filter = True
         elif fields['passed_filtering'] == 'N':
            record.passed_filter = False
         else:
            raise ValueError, "Illegal 'passed filter' value in Solexa export data: '%s'." % fields['passed_filtering']
         record.index_string = fields['index_string']
         if fields['pos'] == '':
            record.iv = None
            record.nomatch_code = fields['chrom']
         else:
            if fields['strand'] == 'F':
               strand = '+'
            elif fields['strand'] == 'R':
               strand = '-'
            else:
               raise ValueError, "Illegal strand value in Solexa export data."
            start = int( fields['pos'] )
            chrom = fields['chrom']
            if fields['chrom'] == "":
               chrom = fields['contig']
            record.iv = GenomicInterval( chrom, start, 
               start + len( fields['read_seq'] ), strand )
         yield record
 
class SAM_Reader( FileOrSequence ):
   """A SAM_Reader object is associated with a SAM file that 
   contains short read alignments. It can generate an iterator of Alignment
   objects."""

   def __iter__( self ):
      for line in FileOrSequence.__iter__( self ):
         if line.startswith( "@" ):
	    # do something with the header line
	    continue
         algnt = SAM_Alignment( line )
         yield algnt
       


class GenomicArrayOfSets( GenomicArray ):

   """A GenomicArrayOfSets is a specialization of GenomicArray that allows to store
   sets of objects. On construction, the step vectors are initialized with empty sets.
   By using the 'add_value' method, objects can be added to intervals. If an object
   is already present in the set(s) at this interval, an the new object is added to
   the present set, and the set is split if necessary.
   """

   def __init__( self, chroms, stranded=True ):
      GenomicArray.__init__( self, chroms, stranded, 'O' )

   def add_chrom( self, chrom, length = sys.maxint, start_index = 0 ):
      GenomicArray.add_chrom( self, chrom, length, start_index )
      if self.stranded:
         self.step_vectors[ chrom ][ "+" ][ : ] = set()
         self.step_vectors[ chrom ][ "-" ][ : ] = set()
      else:
         self.step_vectors[ chrom ][ : ] = set()
      
   def add_value( self, value, iv ):
      def _f( oldset ):
         newset = oldset.copy()
         newset.add( value )
         return newset 
      self.apply( _f, iv )
      
       
###########################
##   paired-end handling
###########################
         
def pair_SAM_alignments( alignments ):

   def check_is_pe( read ):
      if not read.paired_end:
         raise ValueError, "'pair_alignments' needs a sequence of paired-end alignments"

   def process_paired_reads( read1, read2 ):   
      if read1.pe_which == "second":
         aux = read1
         read1 = read2
         read2 = aux
      if not ( ( read1.pe_which == "first" and read2.pe_which == "second" ) or 
               ( read1.pe_which == "unknown" and read2.pe_which == "unknown" ) ):
         warnings.warn( "Incorrect first/second assignments in mate pairs " + 
            read1.read.name )
      if not ( read1.proper_pair and read2.proper_pair ):
         warnings.warn( "Incorrect 'proper_pair' flag value for read pair " + 
            read1.read.name )
      if not ( read1.mate_start == read2.iv.start_as_pos and 
            read2.mate_start == read1.iv.start_as_pos ):
         warnings.warn( "Read pair " + read1.read.name +
            " show inconsistency between 'iv' and 'mate_start' values" )
      return ( read1, read2 )

   def process_single_read( read ):
      if read.mate_aligned:         
         warnings.warn( "Read " + read.read.name + " claims to have an aligned mate " +
            "which could not be found. (Is the SAM file properly sorted?)" )
      if read.pe_which == "second":
         return ( None, read )
      else:
         return ( read, None )

   alignments = iter( alignments )
   read1 = None
   read2 = None
   while True:
      if read1 is None:
         read1 = alignments.next()
         check_is_pe( read1 )
      read2 = alignments.next()
      check_is_pe( read2 )      
      if read1.read.name == read2.read.name:
         yield process_paired_reads( read1, read2 )
         read1 = None
         read2 = None
      else:
         yield process_single_read( read1 )
         read1 = read2
         read2 = None
                
