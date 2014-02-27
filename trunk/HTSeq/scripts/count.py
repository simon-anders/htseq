import sys, optparse, itertools, warnings, traceback, os.path

import HTSeq

class UnknownChrom( Exception ):
   pass
   
def invert_strand( iv ):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

def count_reads_in_features( sam_filename, gff_filename, samtype, order, stranded, 
      overlap_mode, feature_type, id_attribute, quiet, minaqual, samout ):
      
   def write_to_samout( r, assignment ):
      if samoutfile is None:
         return
      if not pe_mode:
         r = (r,)
      for read in r:
         if read is not None:
            samoutfile.write( read.original_sam_line.rstrip() + 
               "\tXF:Z:" + assignment + "\n" )
      

   if samout != "":
      samoutfile = open( samout, "w" )
   else:
      samoutfile = None
      
   features = HTSeq.GenomicArrayOfSets( "auto", stranded != "no" )     
   counts = {}

   # Try to open samfile to fail early in case it is not there
   if sam_filename != "-":
      open( sam_filename ).close()
      
   gff = HTSeq.GFF_Reader( gff_filename )   
   i = 0
   try:
      for f in gff:
         if f.type == feature_type:
            try:
               feature_id = f.attr[ id_attribute ]
            except KeyError:
               raise ValueError, ( "Feature %s does not contain a '%s' attribute" % 
                  ( f.name, id_attribute ) )
            if stranded != "no" and f.iv.strand == ".":
               raise ValueError, ( "Feature %s at %s does not have strand information but you are "
                  "running htseq-count in stranded mode. Use '--stranded=no'." % 
                  ( f.name, f.iv ) )
            features[ f.iv ] += feature_id
            counts[ f.attr[ id_attribute ] ] = 0
         i += 1
         if i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d GFF lines processed.\n" % i )
   except:
      sys.stderr.write( "Error occured when processing GFF file (%s):\n" % gff.get_line_number_string() )
      raise
      
   if not quiet:
      sys.stderr.write( "%d GFF lines processed.\n" % i )
      
   if len( counts ) == 0:
      sys.stderr.write( "Warning: No features of type '%s' found.\n" % feature_type )
   
   if samtype == "sam":
      SAM_or_BAM_Reader = HTSeq.SAM_Reader
   elif samtype == "bam":
      SAM_or_BAM_Reader = HTSeq.BAM_Reader
   else:
      raise ValueError, "Unknown input format %s specified." % samtype

   try:
      if sam_filename != "-":
         read_seq_file = SAM_or_BAM_Reader( sam_filename )
         read_seq = read_seq_file
         first_read = iter(read_seq).next()
      else:
         read_seq_file = SAM_or_BAM_Reader( sys.stdin )
         read_seq_iter = iter( read_seq_file )
         first_read = read_seq_iter.next()
         read_seq = itertools.chain( [ first_read ], read_seq_iter )
      pe_mode = first_read.paired_end
   except:
      sys.stderr.write( "Error occured when reading beginning of SAM/BAM file.\n" )
      raise

   try:
      if pe_mode:
         if order == "name":
            read_seq = HTSeq.pair_SAM_alignments( read_seq )
         elif order == "pos":
            read_seq = HTSeq.pair_SAM_alignments_with_buffer( read_seq )
         else:
            raise ValueError, "Illegal order specified."
      empty = 0
      ambiguous = 0
      notaligned = 0
      lowqual = 0
      nonunique = 0
      i = 0   
      for r in read_seq:
         if i > 0 and i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d SAM alignment record%s processed.\n" % ( i, "s" if not pe_mode else " pairs" ) )

         i += 1
         if not pe_mode:
            if not r.aligned:
               notaligned += 1
               write_to_samout( r, "__not_aligned" )
               continue
            try:
               if r.optional_field( "NH" ) > 1:
                  nonunique += 1
                  write_to_samout( r, "__alignment_not_unique" )
                  continue
            except KeyError:
               pass
            if r.aQual < minaqual:
               lowqual += 1
               write_to_samout( r, "__too_low_aQual" )
               continue
            if stranded != "reverse":
               iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0 )
            else:
               iv_seq = ( invert_strand( co.ref_iv ) for co in r.cigar if co.type == "M" and co.size > 0 )            
         else:
            if r[0] is not None and r[0].aligned:
               if stranded != "reverse":
                  iv_seq = ( co.ref_iv for co in r[0].cigar if co.type == "M" and co.size > 0 )
               else:
                  iv_seq = ( invert_strand( co.ref_iv ) for co in r[0].cigar if co.type == "M" and co.size > 0 )
            else:
               iv_seq = tuple()
            if r[1] is not None and r[1].aligned:            
               if stranded != "reverse":
                  iv_seq = itertools.chain( iv_seq, 
                     ( invert_strand( co.ref_iv ) for co in r[1].cigar if co.type == "M" and co.size > 0 ) )
               else:
                  iv_seq = itertools.chain( iv_seq, 
                     ( co.ref_iv for co in r[1].cigar if co.type == "M" and co.size > 0 ) )
            else:
               if ( r[0] is None ) or not ( r[0].aligned ):
                  write_to_samout( r, "__not_aligned" )
                  notaligned += 1
                  continue         
            try:
               if ( r[0] is not None and r[0].optional_field( "NH" ) > 1 ) or \
                     ( r[1] is not None and r[1].optional_field( "NH" ) > 1 ):
                  nonunique += 1
                  write_to_samout( r, "__alignment_not_unique" )
                  continue
            except KeyError:
               pass
            if ( r[0] and r[0].aQual < minaqual ) or ( r[1] and r[1].aQual < minaqual ):
               lowqual += 1
               write_to_samout( r, "__too_low_aQual" )
               continue         
         
         try:
            if overlap_mode == "union":
               fs = set()
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     fs = fs.union( fs2 )
            elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
               fs = None
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     if len(fs2) > 0 or overlap_mode == "intersection-strict":
                        if fs is None:
                           fs = fs2.copy()
                        else:
                           fs = fs.intersection( fs2 )
            else:
               sys.exit( "Illegal overlap mode." )
            if fs is None or len( fs ) == 0:
               write_to_samout( r, "__no_feature" )
               empty += 1
            elif len( fs ) > 1:
               write_to_samout( r, "__ambiguous[" + '+'.join( fs ) + "]" )
               ambiguous += 1
            else:
               write_to_samout( r, list(fs)[0] )
               counts[ list(fs)[0] ] += 1
         except UnknownChrom:
            write_to_samout( r, "__no_feature" )
            empty += 1

   except:
      sys.stderr.write( "Error occured when processing SAM input (%s):\n" % read_seq_file.get_line_number_string() )
      raise

   if not quiet:
      sys.stderr.write( "%d SAM %s processed.\n" % ( i, "alignments " if not pe_mode else "alignment pairs" ) )
         
   if samoutfile is not None:
      samoutfile.close()

   for fn in sorted( counts.keys() ):
      print "%s\t%d" % ( fn, counts[fn] )
   print "__no_feature\t%d" % empty
   print "__ambiguous\t%d" % ambiguous
   print "__too_low_aQual\t%d" % lowqual
   print "__not_aligned\t%d" % notaligned
   print "__alignment_not_unique\t%d" % nonunique

      
def main():
   
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] alignment_file gff_file",
      
      description=
         "This script takes an alignment file in SAM/BAM format and a " +
         "feature file in GFF format and calculates for each feature " +
         "the number of reads mapping to it. See " +
         "http://www-huber.embl.de/users/anders/HTSeq/doc/count.html for details." ,
         
      epilog = 
         "Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology " +
         "Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General " +
         "Public License v3. Part of the 'HTSeq' framework, version %s." % HTSeq.__version__ )
         
   optParser.add_option( "-f", "--format", type="choice", dest="samtype",
      choices = ( "sam", "bam" ), default = "sam",
      help = "type of <alignment_file> data, either 'sam' or 'bam' (default: sam)" )

   optParser.add_option( "-r", "--order", type="choice", dest="order",
     choices=("pos", "name"), default="name",
     help = "'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing " +
        "data must be sorted either by position or by read name, and the sorting order " +
        "must be specified. Ignored for single-end data." )

   optParser.add_option( "-s", "--stranded", type="choice", dest="stranded",
     choices = ( "yes", "no", "reverse" ), default = "yes",
     help = "whether the data is from a strand-specific assay. Specify 'yes', " +
        "'no', or 'reverse' (default: yes). " +
        "'reverse' means 'yes' with reversed strand interpretation" )
      
   optParser.add_option( "-a", "--minaqual", type="int", dest="minaqual",
      default = 10,
      help = "skip all reads with alignment quality lower than the given " +
         "minimum value (default: 10)" )
      
   optParser.add_option( "-t", "--type", type="string", dest="featuretype",
      default = "exon", help = "feature type (3rd column in GFF file) to be used, " +
         "all features of other type are ignored (default, suitable for Ensembl " +
         "GTF files: exon)" )
         
   optParser.add_option( "-i", "--idattr", type="string", dest="idattr",
      default = "gene_id", help = "GFF attribute to be used as feature ID (default, " +
      "suitable for Ensembl GTF files: gene_id)" )

   optParser.add_option( "-m", "--mode", type="choice", dest="mode",
      choices = ( "union", "intersection-strict", "intersection-nonempty" ), 
      default = "union", help = "mode to handle reads overlapping more than one feature " +
         "(choices: union, intersection-strict, intersection-nonempty; default: union)" )
         
   optParser.add_option( "-o", "--samout", type="string", dest="samout",
      default = "", help = "write out all SAM alignment records into an output " +
      "SAM file called SAMOUT, annotating each line with its feature assignment " +
      "(as an optional field with tag 'XF')" )

   optParser.add_option( "-q", "--quiet", action="store_true", dest="quiet",
      help = "suppress progress report" ) # and warnings" )

   if len( sys.argv ) == 1:
      optParser.print_help()
      sys.exit(1)

   (opts, args) = optParser.parse_args()

   if len( args ) != 2:
      sys.stderr.write( sys.argv[0] + ": Error: Please provide two arguments.\n" )
      sys.stderr.write( "  Call with '-h' to get usage information.\n" )
      sys.exit( 1 )
      
   warnings.showwarning = my_showwarning
   try:
      count_reads_in_features( args[0], args[1], opts.samtype, opts.order, opts.stranded, 
         opts.mode, opts.featuretype, opts.idattr, opts.quiet, opts.minaqual,
         opts.samout )
   except:
      sys.stderr.write( "  %s\n" % str( sys.exc_info()[1] ) )
      sys.stderr.write( "  [Exception type: %s, raised in %s:%d]\n" % 
         ( sys.exc_info()[1].__class__.__name__, 
           os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
           traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
      sys.exit( 1 )

def my_showwarning( message, category, filename, lineno = None, line = None ):
   sys.stderr.write( "Warning: %s\n" % message )

if __name__ == "__main__":
   main()

