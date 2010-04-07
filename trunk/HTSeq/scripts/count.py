import sys, optparse

import HTSeq

def count_reads_in_features( sam_filename, gff_filename, stranded, 
      overlap_mode, feature_type, id_attribute, quiet ):
      
   features = HTSeq.GenomicArrayOfSets( [], stranded )     
   counts = {}
      
   for f in HTSeq.GFF_Reader( gff_filename ):
      if f.iv.chrom not in features.step_vectors.keys():
	 features.add_chrom( f.iv.chrom )
      if f.type == feature_type:
         try:
            features.add_value( f.attr[ id_attribute ], f.iv )
	 except KeyError:
	    sys.exit( "Feature %s does not contain a '%s' attribute" % 
	       ( f.name, id_attribute ) )
	 counts[ f.attr[ id_attribute ] ] = 0
	 
   if len( counts ) == 0 and not quiet:
      sys.stderr.write( "Warning: No features of type '%s' found.\n" % feature_type )
   
   empty = 0
   ambiguous = 0
   i = 0
   for r in HTSeq.SAM_Reader( sam_filename ):
      if r.aligned:
      
         if r.iv.chrom not in features.step_vectors:
	    sys.stderr.write( ( "Warning: Skipping read '%s', aligned to %s, because " +
	       "chromosome '%s' did not appear in the GFF file.\n" ) % 
	       ( r.read.name, r.iv, r.iv.chrom ) )
	    continue

         if overlap_mode == "union":
            fs = set()
            for fs2 in features.get_steps( r.iv, values_only=True ):
               fs = fs.union( fs2 )
	 elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
	    fs = None
            for fs2 in features.get_steps( r.iv, values_only=True ):
	       if fs is None:
	          fs = fs2
	       else:
        	  if len(fs2) > 0 or overlap_mode == "intersection-strict":
	             fs = fs.intersection( fs2 )
         else:
	    sys.exit( "Illegal overlap mode." )
         if len( fs ) == 0:
            empty += 1
         elif len( fs ) > 1:
            ambiguous += 1
         else:
            counts[ list(fs)[0] ] += 1
      i += 1
      if i % 100000 == 0 and not quiet:
         sys.stderr.write( "%d reads processed.\n" % i )
	 
   for fn in sorted( counts.keys() ):
      print "%s\t%d" % ( fn, counts[fn] )
   print "no_feature\t%d" % empty
   print "ambiguous\t%d" % ambiguous

      
def main():
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] sam_file gff_file",
      
      description=
         "This script takes an alignment file in SAM format and a " +
	 "feature file in GFF format and calculates for each feature " +
	 "the number of reads mapping to it. See " +
	 "http://www-huber.embl.de/users/anders/HTSeq/doc/count.html for details.",
	 
      epilog = 
	 "Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology " +
	 "Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General " +
	 "Public License v3. Part of the 'HTSeq' framework." )
	 
   optParser.add_option( "-m", "--mode", type="choice", dest="mode",
      choices = ( "union", "intersection-strict", "intersection-nonempty" ), 
      default = "union", help = "mode to handle reads overlapping more than one feature" +
         "(choices: union, intersection-strict, intersection-nonempty; default: union)" )
	 
   optParser.add_option( "-t", "--type", type="string", dest="featuretype",
      default = "exon", help = "feature type (3rd column in GFF file) to be used, " +
         "all features of other type are ignored (default, suitable for Ensembl " +
	 "GTF files: exon)" )
	 
   optParser.add_option( "-i", "--idattr", type="string", dest="idattr",
      default = "gene_id", help = "GFF attribute to be used as feature ID (default, " +
      "suitable for Ensembl GTF files: gene_id)" )

   optParser.add_option( "-s", "--stranded", type="choice", dest="stranded",
      choices = ( "yes", "no" ), default = "yes",
      help = "whether the data is from a strand-specific assay (default: yes)" )
      
   optParser.add_option( "-q", "--quiet", action="store_true", dest="quiet",
      help = "suppress progress report" )

   if len( sys.argv ) == 1:
      optParser.print_help()
      sys.exit(1)

   (opts, args) = optParser.parse_args()

   if len( args ) != 2:
      sys.stderr.write( sys.argv[0] + ": Error: Please provide two arguments.\n" )
      sys.stderr.write( "  Call with '-h' to get usage information.\n" )
      sys.exit( 1 )
      
      
   count_reads_in_features( args[0], args[1], opts.stranded == "yes", 
      opts.mode, opts.featuretype, opts.idattr, opts.quiet )


if __name__ == "__main__":
   main()
