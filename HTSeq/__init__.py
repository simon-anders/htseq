"""HTSeq is a package to process high-throughput sequencing data.

See http://www-huber.embl.de/users/anders/HTSeq for documentation.
"""

import itertools, warnings, os, shlex

try:
   from _HTSeq import *
except ImportError:
   if os.path.isfile( "setup.py" ):
      raise ImportError( "Cannot import 'HTSeq' when working directory is HTSeq's own build directory.")
   else:
      raise
      
from _version import __version__

#from vcf_reader import *

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
      self.line_no = None
      
   def __iter__( self ):
      self.line_no = 1
      if isinstance( self.fos, str ):
         if self.fos.lower().endswith( ( ".gz" , ".gzip" ) ):
            lines = gzip.open( self.fos )
         else:
            lines = open( self.fos )
      else:
         lines = self.fos
      for line in lines:
         yield line
         self.line_no += 1
      if isinstance( self.fos, str ):
         lines.close()
      self.line_no = None
         
   def __repr__( self ):
      if isinstance( self.fos, str ):
         return "<%s object, connected to file name '%s'>" % (
            self.__class__.__name__, self.fos )
      else:
         return "<%s object, connected to %s >" % (
            self.__class__.__name__, repr( self.fos ) )   
            
   def get_line_number_string( self ):
      if self.line_no is None:
         if isinstance( self.fos, str ):
            return "file %s closed" % self.fos
         else:
            return "file closed"
      if isinstance( self.fos, str ):
         return "line %d of file %s" % ( self.line_no, self.fos )
      else:
         return "line %d" % self.line_no
         

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
         score = self.score
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
      attr_str = '; '.join( [ '%s%s\"%s\"' % ( ak, sep, attr[ak] ) for ak in attr ] )
      return "\t".join( str(a) for a in ( self.iv.chrom, source, 
         self.type, self.iv.start+1, self.iv.end, score, 
         self.iv.strand, frame, attr_str ) ) + "\n"
         

_re_attr_main = re.compile( "\s*([^\s\=]+)[\s=]+(.*)" )
_re_attr_empty = re.compile( "^\s*$" )

def parse_GFF_attribute_string( attrStr, extra_return_first_value=False ):
   """Parses a GFF attribute string and returns it as a dictionary.
   
   If 'extra_return_first_value' is set, a pair is returned: the dictionary
   and the value of the first attribute. This might be useful if this is the ID.
   """
   if attrStr.endswith( "\n" ):
      attrStr = attrStr[:-1]
   d = {}
   first_val = "_unnamed_"
   for (i, attr) in itertools.izip( itertools.count(), _HTSeq.quotesafe_split( attrStr ) ):
      if _re_attr_empty.match( attr ):
         continue
      if attr.count( '"' ) not in ( 0, 2 ):
         raise ValueError, "The attribute string seems to contain mismatched quotes."
      mo = _re_attr_main.match( attr )
      if not mo:
         raise ValueError, "Failure parsing GFF attribute line"
      val = mo.group(2)
      if val.startswith( '"' ) and val.endswith( '"' ):
         val = val[1:-1]
      #val = urllib.unquote( val )
      d[ intern(mo.group(1)) ] = intern(val)
      if extra_return_first_value and i == 0:
         first_val = val         
   if extra_return_first_value:
      return ( d, first_val )
   else:
      return d

_re_gff_meta_comment = re.compile( "##\s*(\S+)\s+(\S*)" )

class GFF_Reader( FileOrSequence ):

   """Parse a GFF file
   
   Pass the constructor either a file name or an iterator of lines of a 
   GFF files. If a file name is specified, it may refer to a gzip compressed
   file.
   
   Iterating over the object then yields GenomicFeature objects.
   """
   
   def __init__( self, filename_or_sequence, end_included=True ):
      FileOrSequence.__init__( self, filename_or_sequence )
      self.end_included = end_included
      self.metadata = {}
   
   
   def __iter__( self ):
      for line in FileOrSequence.__iter__( self ):
         if line == "\n":
            continue
         if line.startswith( '#' ):
            if line.startswith( "##" ):
               mo = _re_gff_meta_comment.match( line )
               if mo:
                  self.metadata[ mo.group(1) ] = mo.group(2)
            continue
         ( seqname, source, feature, start, end, score, 
            strand, frame, attributeStr ) = line.split( "\t", 8 )   
         ( attr, name ) = parse_GFF_attribute_string( attributeStr, True )
         if self.end_included:
            iv = GenomicInterval( seqname, int(start)-1, int(end), strand )
         else:
            iv = GenomicInterval( seqname, int(start)-1, int(end)-1, strand )
         f = GenomicFeature( name, feature, iv )
         if score != ".":
            score = float( score )
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
      if seq is not None:
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
            mo = _re_fasta_header_line.match( line )
            seqname = mo.group(1)
            length = 0
         else: 
            assert seqname is not None, "FASTA file does not start with '>'."
            length += len( line.rstrip() )
      if seqname is not None:
         seqlengths[ seqname ] = length
      return seqlengths
      
   @staticmethod   
   def _import_pysam():
      global pysam
      try:
         import pysam
      except ImportError:
         sys.stderr.write( "Please install the 'pysam' package to be able to use the Fasta indexing functionality." )
         raise
      
   def build_index( self, force = False ):
      self._import_pysam()
      if not isinstance( self.fos, str ):
         raise TypeError, "This function only works with FastaReader objects " + \
            "connected to a fasta file via file name"
      index_filename = self.fos + ".fai"
      if os.access( index_filename, os.R_OK ):
         if (not force) and os.stat( self.filename_or_sequence ).st_mtime <= \
               os.stat( index_filename ).st_mtime:
            # index is up to date
            return
      pysam.faidx( self.fos )
      if not os.access( index_filename, os.R_OK ):
         raise SystemError, "Building of Fasta index failed due to unknown error."
      
   def __getitem__( self, iv ):
      if not isinstance( iv, GenomicInterval ):
         raise TypeError, "GenomicInterval expected as key."
      if not isinstance( self.fos, str ):
         raise TypeError, "This function only works with FastaReader objects " + \
            "connected to a fasta file via file name"
      self._import_pysam()
      fasta = pysam.faidx( self.fos, "%s:%d-%d" % ( iv.chrom, iv.start, iv.end-1 ) )
      ans = list( FastaReader( fasta ) )
      assert len( ans ) == 1
      ans[0].name = str(iv)
      if iv.strand != "-":
         return ans[0]
      else:
         return ans[0].get_reverse_complement()
            
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
         self.qualscale = "solexa"

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
         try:
            algnt = SAM_Alignment.from_SAM_line( line )
         except ValueError, e:
            e.args = e.args + ( self.get_line_number_string(), )
            raise
         yield algnt
       
class GenomicArrayOfSets( GenomicArray ):

   """A GenomicArrayOfSets is a specialization of GenomicArray that allows to store
   sets of objects. On construction, the step vectors are initialized with empty sets.
   By using the 'add_value' method, objects can be added to intervals. If an object
   is already present in the set(s) at this interval, an the new object is added to
   the present set, and the set is split if necessary.
   """

   def __init__( self, chroms, stranded=True, storage='step', memmap_dir = "" ):
      GenomicArray.__init__( self, chroms, stranded, 'O', storage, memmap_dir )

   def add_chrom( self, chrom, length = sys.maxint, start_index = 0 ):
      GenomicArray.add_chrom( self, chrom, length, start_index )
      for cv in self.chrom_vectors[ chrom ].values():
         cv[:] = set()
         cv.is_vector_of_sets = True
      
       
###########################
##   paired-end handling
###########################
                
                
def pair_SAM_alignments( alignments, bundle=False ):

   mate_missing_count = [0]

   def process_list( almnt_list ):
      while len( almnt_list ) > 0:
         a1 = almnt_list.pop( 0 )
         # Find its mate
         for a2 in almnt_list:
            if a1.pe_which == a2.pe_which:
               continue
            if a1.aligned != a2.mate_aligned or a1.mate_aligned != a2.aligned:
               continue
            if not (a1.aligned and a2.aligned):
               break
            if a1.iv.chrom == a2.mate_start.chrom and a1.iv.start == a2.mate_start.pos and \
                  a2.iv.chrom == a1.mate_start.chrom and a2.iv.start == a1.mate_start.pos:
               break
         else:
            if a1.mate_aligned:
               mate_missing_count[0] += 1
               if mate_missing_count[0] == 1:
                  warnings.warn( "Read " + a1.read.name + " claims to have an aligned mate " +
                     "which could not be found in an adjacent line." )
            a2 = None
         if a2 is not None:
            almnt_list.remove( a2 )
         if a1.pe_which == "first":
            yield ( a1, a2 )
         else:
            assert a1.pe_which == "second"
            yield ( a2, a1 )

   almnt_list = []
   current_name = None
   for almnt in alignments:
      if not almnt.paired_end:
         raise ValueError, "'pair_alignments' needs a sequence of paired-end alignments"
      if almnt.pe_which == "unknown":
         raise ValueError, "Paired-end read found with 'unknown' 'pe_which' status."
      if almnt.read.name == current_name:
         almnt_list.append( almnt )
      else:
         if bundle:
            yield list( process_list( almnt_list ) )
         else:
            for p in process_list( almnt_list ):
               yield p
         current_name = almnt.read.name
         almnt_list = [ almnt ]
   if bundle:
      yield list( process_list( almnt_list ) )
   else:
      for p in process_list( almnt_list ):
         yield p
   if mate_missing_count[0] > 1:
      warnings.warn( "%d reads with missing mate encountered." % mate_missing_count[0] )


def pair_SAM_alignments_with_buffer( alignments, max_buffer_size=3000000 ):

   almnt_buffer = {}
   ambiguous_pairing_counter = 0
   for almnt in alignments:

      if not almnt.paired_end:
         raise ValueError, "Sequence of paired-end alignments expected, but got single-end alignment."
      if almnt.pe_which == "unknown":
         raise ValueError, "Cannot process paired-end alignment found with 'unknown' 'pe_which' status."

      matekey = ( 
         almnt.read.name, 
         "second" if almnt.pe_which == "first" else "first",
         almnt.mate_start.chrom if almnt.mate_aligned else None, 
         almnt.mate_start.pos if almnt.mate_aligned else None, 
         almnt.iv.chrom if almnt.aligned else None, 
         almnt.iv.start if almnt.aligned else None, 
         -almnt.inferred_insert_size if almnt.aligned and almnt.mate_aligned else None )

      if matekey in almnt_buffer:
         if len( almnt_buffer[ matekey ] ) == 1:            
            mate = almnt_buffer[ matekey ][ 0 ]
            del almnt_buffer[ matekey ]
         else:
            mate = almnt_buffer[ matekey ].pop( 0 )
            if ambiguous_pairing_counter == 0:
               ambiguous_pairing_first_occurance = matekey
            ambiguous_pairing_counter += 1
         if almnt.pe_which == "first":
            yield ( almnt, mate )
         else:
            yield ( mate, almnt )
      else:
         almntkey = ( 
            almnt.read.name, almnt.pe_which, 
            almnt.iv.chrom if almnt.aligned else None, 
            almnt.iv.start if almnt.aligned else None, 
            almnt.mate_start.chrom if almnt.mate_aligned else None, 
            almnt.mate_start.pos if almnt.mate_aligned else None, 
            almnt.inferred_insert_size if almnt.aligned and almnt.mate_aligned else None )
         if almntkey not in almnt_buffer:
            almnt_buffer[ almntkey ] = [ almnt ]
         else:
            almnt_buffer[ almntkey ].append( almnt )
         if len(almnt_buffer) > max_buffer_size:
            raise ValueError, "Maximum alignment buffer size exceeded while pairing SAM alignments."

   if len(almnt_buffer) > 0:
      warnings.warn( "Mate records missing for %d records; first such record: %s." % 
         ( len(almnt_buffer), str( almnt_buffer.values()[0][0] ) ) )
      for almnt_list in almnt_buffer.values():
         for almnt in almnt_list:
            if almnt.pe_which == "first":
               yield ( almnt, None )
            else:
               yield ( None, almnt )

   if ambiguous_pairing_counter > 0:
      warnings.warn( "Mate pairing was ambiguous for %d records; mate key for first such record: %s." %
         ( ambiguous_pairing_counter, str( ambiguous_pairing_first_occurance ) ) )


###########################
##   variant calls
###########################


_re_vcf_meta_comment = re.compile( "^##([a-zA-Z]+)\=(.*)$" )

_re_vcf_meta_descr = re.compile('ID=[^,]+,?|Number=[^,]+,?|Type=[^,]+,?|Description="[^"]+",?')

_re_vcf_meta_types = re.compile( "[INFO|FILTER|FORMAT]" )

_vcf_typemap = {
    "Integer":int,
    "Float":float,
    "String":str,
    "Flag":bool
}

class VariantCall( object ):
    
    def __init__( self, chrom = None, pos = None, identifier = None, ref = None, alt = None, qual = None, filtr = None, info = None ):
        self.chrom  = chrom
        self.pos    = pos
        self.id     = identifier
        self.ref    = ref
        self.alt    = alt
        self.qual   = qual
        self.filter = filtr
        self.info   = info
        self._original_line = None
    
    @classmethod
    def fromdict( cls, dictionary ):
        ret = cls()
        ret.chrom   = dictionary["chrom"]
        ret.pos     = dictionary["pos"]
        ret.id      = dictionary["id"]
        ret.ref     = dictionary["ref"]
        ret.alt     = dictionary["alt"]
        ret.qual    = dictionary["qual"]
        ret.filter  = dictionary["filter"]
        ret.info    = dictionary["info"]
        ret._original_line = None

    @classmethod
    def fromline( cls, line, nsamples = 0, sampleids = [] ):
        ret = cls()
        if nsamples == 0:
            ret.format = None
            ret.chrom, ret.pos, ret.id, ret.ref, ret.alt, ret.qual, ret.filter, ret.info = line.rstrip("\n").split("\t", 7)
        else:
            lsplit = line.rstrip("\n").split("\t")
            ret.chrom, ret.pos, ret.id, ret.ref, ret.alt, ret.qual, ret.filter, ret.info = lsplit[:8]
            ret.format = lsplit[8].split(":")
            ret.samples = {}
            spos=9
            for sid in sampleids:
                ret.samples[ sid ] = dict( ( name, value ) for (name, value) in itertools.izip( ret.format, lsplit[spos].split(":") ) )
                spos += 1
        ret.pos = GenomicPosition( ret.chrom, int(ret.pos) )
        ret.alt = ret.alt.split(",")
        ret._original_line = line
        return ret
    
    def infoline( self ):
        if self.info.__class__ == dict:
            return ";".join(map((lambda key: str(key) + "=" + str(self.info[key])), self.info ))
        else:
            return self.info
    
    def get_original_line( self ):
       warnings.warn( "Original line is empty, probably this object was created from scratch and not from a line in a .vcf file!" )
       return self._original_line
    
    def sampleline( self ):
       if self.format == None:
          print >> sys.stderr, "No samples in this variant call!" 
          return ""
       keys = self.format
       ret = [ ":".join( keys ) ]
       for sid in self.samples:
          tmp = []
          for k in keys:
             if k in self.samples[sid]:
                tmp.append( self.samples[sid][k] )
          ret.append( ":".join(tmp) )
       return "\t".join( ret )
    
    def to_line( self ):
       if self.format == None:
          return "\t".join( map( str, [ self.pos.chrom, self.pos.pos, self.id, self.ref, ",".join( self.alt ), self.qual, self.filter, self.infoline() ] ) ) + "\n"
       else:
          return "\t".join( map( str, [ self.pos.chrom, self.pos.pos, self.id, self.ref, ",".join( self.alt ), self.qual, self.filter, self.infoline(), self.sampleline() ] ) ) + "\n"
    
    def __descr__( self ):
        return "<VariantCall at %s, ref '%s', alt %s >" % (str(self.pos).rstrip("/."), self.ref, str(self.alt).strip("[]"))
    
    def __str__( self ):
        return "%s:'%s'->%s" % (str(self.pos).rstrip("/."), self.ref, str(self.alt).strip("[]"))
    
    def unpack_info( self, infodict ):
        tmp = {}
        for token in self.info.strip(";").split(";"):
            if re.compile("=").search(token):
                token = token.split("=")
                if infodict.has_key( token[0] ):
                    tmp[token[0]] = map( infodict[token[0]], token[1].split(",") )
                else:
                    tmp[token[0]] = token[1].split(",")
                if len( tmp[ token[0] ] ) == 1:
                    tmp[token[0]] = tmp[token[0]][0]
            else: #Flag attribute found
                tmp[token] = True
        diff = set( infodict.keys() ).difference( set( tmp.keys() ) )
        for key in diff:
            if infodict[key] == bool:
                tmp[key] = False
        self.info = tmp

class VCF_Reader( FileOrSequence ):

    def __init__( self, filename_or_sequence ):
        FileOrSequence.__init__( self, filename_or_sequence )
        self.metadata = {}
        self.info = {}
        self.filters = {}
        self.formats = {}
        self.nsamples = 0
        self.sampleids = []
        
    def make_info_dict( self ):
        self.infodict = dict( ( key, _vcf_typemap[self.info[key]["Type"]] ) for key in self.info.keys() )
    
    def parse_meta( self, header_filename = None ):
        if header_filename == None:
            the_iter = FileOrSequence.__iter__( self )
        else:
            the_iter = open( header_filename, "r" )
        
        for line in the_iter:
            if line.startswith( '#' ):
                if line.startswith( "##" ):
                    mo = _re_vcf_meta_comment.match( line )
                    if mo:
                        value = mo.group(2)
                        if mo.group(1) == "INFO":
                            value = dict( e.rstrip(",").split("=",1) for e in _re_vcf_meta_descr.findall(value) )
                            key = value["ID"]
                            del value["ID"]
                            self.info[ key ] = value
                        elif mo.group(1) == "FILTER":
                            value = dict( e.rstrip(",").split("=",1) for e in _re_vcf_meta_descr.findall(value) )
                            key = value["ID"]
                            del value["ID"]
                            self.filters[ key ] = value
                        elif mo.group(1) == "FORMAT":
                            value = dict( e.rstrip(",").split("=",1) for e in _re_vcf_meta_descr.findall(value) )
                            key = value["ID"]
                            del value["ID"]
                            self.formats[ key ] = value
                        else:
                            self.metadata[ mo.group(1) ] = mo.group(2)
                else:
                    self.sampleids = line.rstrip("\t\n").split("\t")[9:]
                    self.nsamples = len( self.sampleids )
                continue
            else:
                break
    
    def meta_info( self, header_filename = None ):
       ret = []
       if header_filename == None:
          the_iter = FileOrSequence.__iter__( self )
       else:
          the_iter = open( header_filename, "r" )
       
       for line in the_iter:
          if line.startswith( '#' ):
             ret.append( line )
          else:
             break
       return ret
    
    def __iter__( self ):
        for line in FileOrSequence.__iter__( self ):
            if line == "\n" or line.startswith( '#' ):
                continue
            vc = VariantCall.fromline( line, self.nsamples, self.sampleids )
            yield vc


class WiggleReader( FileOrSequence ):

    def __init__( self, filename_or_sequence, verbose = True ):
        FileOrSequence.__init__( self, filename_or_sequence )
        self.attributes = {}       
        self.stepType = 'none'
        self.verbose = verbose
        
    def __iter__( self ):
        span = 1
        pos = None
        step = None
        chrom = None
        for line in FileOrSequence.__iter__( self ):
            if line.startswith( 'track' ):
                fields = shlex.split(line)[1:]
                self.attributes = dict([(p[0], p[1].strip('"')) for p in [x.split("=") for x in fields]])
            elif line.startswith( 'fixedStep' ): # do fixed step stuff
                self.stepType = 'fixed'
                fields = shlex.split(line)[1:]
                declarations = dict([(p[0], p[1].strip('"')) for p in [x.split("=") for x in fields]])
                pos = int(declarations['start'])
                step = int(declarations['step'])
                chrom = declarations['chrom']
                if 'span' in declarations:
                    span = int(declarations['span'])
                else:
                    span = 1
            elif line.startswith( 'variableStep' ): # do variable step stuff
                self.stepType = 'variable'
                fields = shlex.split(line)[1:]
                declarations = dict([(p[0], p[1].strip('"')) for p in [x.split("=") for x in fields]])
                chrom = declarations['chrom']
                if 'span' in declarations:
                    span = int(declarations['span'])
                else:
                    span = 1
            elif line.startswith( 'browser' ) or line.startswith( '#' ): #Comment or ignored
                if self.verbose:
                    print "Ignored line:", line
                continue
            else:
                if self.stepType == 'fixed':
                    yield ( GenomicInterval( chrom, pos, pos + span, '.' ), float(line.strip()) )
                    pos += step
                elif self.stepType == 'variable':
                    tmp = line.strip().split(" ")
                    pos = int(tmp[0])
                    yield ( GenomicInterval( chrom, pos, pos + span, '.' ), float(tmp[1]) )
            
class BAM_Reader( object ):

    def __init__( self, filename ):
        global pysam
        self.filename = filename
        self.sf = None  # This one is only used by __getitem__
        self.record_no = -1
        try:
           import pysam
        except ImportError:
           sys.stderr.write( "Please Install PySam to use the BAM_Reader Class (http://code.google.com/p/pysam/)" )
           raise
    
    def __iter__( self ):
        sf = pysam.Samfile(self.filename, "rb")
        self.record_no = 0
        for pa in sf:
            yield SAM_Alignment.from_pysam_AlignedRead( pa, sf )
            self.record_no += 1
    
    def fetch( self, reference = None, start = None, end = None, region = None ):
        sf = pysam.Samfile(self.filename, "rb")
        self.record_no = 0
        try:
           for pa in sf.fetch( reference, start, end, region ):
              yield SAM_Alignment.from_pysam_AlignedRead( pa, sf )
              self.record_no += 1
        except ValueError as e:
           if e.message == "fetch called on bamfile without index":
              print "Error: ", e.message
              print "Your bam index file is missing or wrongly named, convention is that file 'x.bam' has index file 'x.bam.bai'!"
           else:
              raise
        except:
           raise

    def get_line_number_string( self ):
        if self.record_no == -1:
            return "unopened file %s" % ( self.filename )
        else:
            return "record #%d in file %s" % ( self.record_no, self.filename )
    
    def __getitem__( self, iv ):
        if not isinstance( iv, GenomicInterval ):
           raise TypeError, "Use a HTSeq.GenomicInterval to access regions within .bam-file!"        
        if self.sf is None:
           self.sf = pysam.Samfile( self.filename, "rb" )
           if not self.sf._hasIndex():
              raise ValueError, "The .bam-file has no index, random-access is disabled!"
        for pa in self.sf.fetch( iv.chrom, iv.start+1, iv.end ):
            yield SAM_Alignment.from_pysam_AlignedRead( pa, self.sf )
    
    def get_header_dict( self ):
       sf = pysam.Samfile(self.filename, "rb")
       return sf.header
    
               
class BAM_Writer( object ):
   def __init__( self, filename, template = None, referencenames = None, referencelengths = None, text = None, header = None ):
      try:
         import pysam
      except ImportError:
         sys.stderr.write( "Please Install PySam to use the BAM_Writer Class (http://code.google.com/p/pysam/)" )
         raise
      
      self.filename = filename
      self.template = template
      self.referencenames = referencenames
      self.referencelengths = referencelengths
      self.text = text
      self.header = header
      self.sf = pysam.Samfile( self.filename, mode="wb", template = self.template, referencenames = self.referencenames, referencelengths = self.referencelengths, text = self.text, header = self.header )
      
   @classmethod
   def from_BAM_Reader( cls, fn, br ):
      return BAM_Writer( filename = fn, header = br.get_header_dict() )
   
   def write( self, alnmt):
      self.sf.write( alnmt.to_pysam_AlignedRead( self.sf ) )
   
   def close( self ):
      self.sf.close()


class BED_Reader( FileOrSequence ):

   def __init__( self, filename_or_sequence ):
      FileOrSequence.__init__( self, filename_or_sequence )
        
   def __iter__( self ):
      for line in FileOrSequence.__iter__( self ):
         if line.startswith( "track" ):
            continue
         fields = line.split()
         if len(fields) < 3:
            raise ValueError, "BED file line contains less than 3 fields"
         if len(fields) > 9:
            raise ValueError, "BED file line contains more than 9 fields"
         iv = GenomicInterval( fields[0], int(fields[1]), int(fields[2]), fields[5] if len(fields) > 5 else "." )
         f = GenomicFeature( fields[3] if len(fields) > 3 else "unnamed", "BED line", iv )
         f.score = float( fields[4] ) if len(fields) > 4 else None
         f.thick = GenomicInterval( iv.chrom, int( fields[6] ), int( fields[7] ), iv.strand ) if len(fields) > 7 else None
         f.itemRgb = [ int(a) for a in fields[8].split(",") ]  if len(fields) > 8 else None
         yield(f)

