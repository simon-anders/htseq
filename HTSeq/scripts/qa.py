#!/usr/bin/env python

# HTSeq_QA.py
#
# (c) Simon Anders, European Molecular Biology Laboratory, 2010
# released under GNU General Public License

import sys, time, os.path, optparse
from itertools import *
import numpy
import HTSeq

def main():

   try:
      import matplotlib
   except ImportError:
      sys.stderr.write("This script needs the 'matplotlib' library, which ")
      sys.stderr.write("was not found. Please install it." )
   matplotlib.use('PDF')
   from matplotlib import pyplot


   # **** Parse command line ****

   optParser = optparse.OptionParser( usage = "%prog [options] read_file",
      description=
	 "This script take a file with high-throughput sequencing reads " +
	 "(supported formats: SAM, Solexa _export.txt, FASTQ, Solexa " +
	 "_sequence.txt) and performs a simply quality assessment by " +
	 "producing plots showing the distribution of called bases and " +
	 "base-call quality scores by position within the reads. The " +
	 "plots are output as a PDF file.",
      epilog = 
	 "Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology " +
	 " Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General " +
	 " Public License v3. Part of the 'HTSeq' framework, version %s." % HTSeq.__version__ )
   optParser.add_option( "-t", "--type", type="choice", dest="type",
      choices = ("sam", "bam", "solexa-export", "fastq", "solexa-fastq"),
      default = "sam", help="type of read_file (one of: sam [default], bam, " +
	 "solexa-export, fastq, solexa-fastq)" )
   optParser.add_option( "-o", "--outfile", type="string", dest="outfile",
      help="output filename (default is <read_file>.pdf)" )
   optParser.add_option( "-r", "--readlength", type="int", dest="readlen",
      help="the maximum read length (when not specified, the script guesses from the file" )
   optParser.add_option( "-g", "--gamma", type="float", dest="gamma", 
      default = 0.3,
      help="the gamma factor for the contrast adjustment of the quality score plot" )
   optParser.add_option( "-n", "--nosplit", action="store_true", dest="nosplit",
      help="do not split reads in unaligned and aligned ones" )
   optParser.add_option( "-m", "--maxqual", type="int", dest="maxqual", default=40,
      help="the maximum quality score that appears in the data (default: 40)" )

   if len( sys.argv ) == 1:
      optParser.print_help()
      sys.exit(1)

   (opts, args) = optParser.parse_args()

   if len( args ) != 1:
      sys.stderr.write( sys.argv[0] + ": Error: Please provide one argument (the read_file).\n" )
      sys.stderr.write( "  Call with '-h' to get usage information.\n" )
      sys.exit( 1 )

   readfilename = args[0]

   if opts.type == "sam":
      readfile = HTSeq.SAM_Reader( readfilename )
      isAlnmntFile = True
   elif opts.type == "bam":
      readfile = HTSeq.BAM_Reader( readfilename )
      isAlnmntFile = True
   elif opts.type == "solexa-export":
      readfile = HTSeq.SolexaExportReader( readfilename )
      isAlnmntFile = True
   elif opts.type == "fastq":
      readfile = HTSeq.FastqReader( readfilename )
      isAlnmntFile = False
   elif opts.type == "solexa-fastq":
      readfile = HTSeq.FastqReader( readfilename, "solexa" )
      isAlnmntFile = False
   else:
      sys.error( "Oops." )

   twoColumns = isAlnmntFile and not opts.nosplit

   if opts.outfile is None:
      outfilename = os.path.basename( readfilename ) + ".pdf"
   else:   
      outfilename = opts.outfile


   # **** Get read length ****

   if opts.readlen is not None:
      readlen = opts.readlen
   else:
      readlen = 0
      if isAlnmntFile:
	 reads = ( a.read for a in readfile )
      else:
	 reads = readfile
      for r in islice( reads, 10000 ):
	 if len( r ) > readlen:
            readlen = len( r )

   max_qual = opts.maxqual
   gamma = opts.gamma


   # **** Initialize count arrays ****

   base_arr_U = numpy.zeros( ( readlen, 5 ), numpy.int )
   qual_arr_U = numpy.zeros( ( readlen, max_qual+1 ), numpy.int )
   if twoColumns:
      base_arr_A = numpy.zeros( ( readlen, 5 ), numpy.int )
      qual_arr_A = numpy.zeros( ( readlen, max_qual+1 ), numpy.int )


   # **** Main counting loop ****

   i = 0
   try:
      for a in readfile:
	 if isAlnmntFile:
	    r = a.read
	 else:
	    r = a
	 if twoColumns and (isAlnmntFile and a.aligned):
	    r.add_bases_to_count_array( base_arr_A )
	    r.add_qual_to_count_array( qual_arr_A )
	 else:
	    r.add_bases_to_count_array( base_arr_U )
	    r.add_qual_to_count_array( qual_arr_U )   
	 i += 1
	 if i % 200000 == 0:
            print i, "reads processed"
   except:
      sys.stderr.write( "Error occured in: %s\n" %
         readfile.get_line_number_string() )
      raise
   print i, "reads processed"


   # **** Normalize result ****

   def norm_by_pos( arr ):
      arr = numpy.array( arr, numpy.float )
      arr_n = ( arr.T / arr.sum( 1 ) ).T
      arr_n[ arr == 0 ] = 0
      return arr_n

   def norm_by_start( arr ):
      arr = numpy.array( arr, numpy.float )
      arr_n = ( arr.T / arr.sum( 1 )[ 0 ] ).T
      arr_n[ arr == 0 ] = 0
      return arr_n


   base_arr_U_n = norm_by_pos( base_arr_U )
   qual_arr_U_n = norm_by_start( qual_arr_U )
   nreads_U = base_arr_U[0,:].sum()
   if twoColumns:
      base_arr_A_n = norm_by_pos( base_arr_A )
      qual_arr_A_n = norm_by_start( qual_arr_A )
      nreads_A = base_arr_A[0,:].sum()


   # **** Make plot ****

   def plot_bases( arr ):
      xg = numpy.arange( readlen )   
      pyplot.plot( xg, arr[ : , 0 ], marker='.', color='red')
      pyplot.plot( xg, arr[ : , 1 ], marker='.', color='darkgreen')
      pyplot.plot( xg, arr[ : , 2 ], marker='.',color='lightgreen')
      pyplot.plot( xg, arr[ : , 3 ], marker='.',color='orange')
      pyplot.plot( xg, arr[ : , 4 ], marker='.',color='grey')
      pyplot.axis( (0, readlen-1, 0, 1 ) )
      pyplot.text( readlen*.70, .9, "A", color="red" )
      pyplot.text( readlen*.75, .9, "C", color="darkgreen" )
      pyplot.text( readlen*.80, .9, "G", color="lightgreen" )
      pyplot.text( readlen*.85, .9, "T", color="orange" )
      pyplot.text( readlen*.90, .9, "N", color="grey" )

   pyplot.figure()
   pyplot.subplots_adjust( top=.85 )
   pyplot.suptitle( os.path.basename(readfilename), fontweight='bold' )

   if twoColumns:

      pyplot.subplot( 221 )
      plot_bases( base_arr_U_n )
      pyplot.ylabel( "proportion of base" )
      pyplot.title( "non-aligned reads\n%.0f%% (%.3f million)" % 
	 ( 100. * nreads_U / (nreads_U+nreads_A), nreads_U / 1e6 ) )

      pyplot.subplot( 222 )
      plot_bases( base_arr_A_n )
      pyplot.title( "aligned reads\n%.0f%% (%.3f million)" % 
	 ( 100. * nreads_A / (nreads_U+nreads_A), nreads_A / 1e6 ) )

      pyplot.subplot( 223 )
      pyplot.pcolor( qual_arr_U_n.T ** gamma, cmap=pyplot.cm.Greens,
	  norm=pyplot.normalize( 0, 1 ) )
      pyplot.axis( (0, readlen-1, 0, max_qual+1 ) )
      pyplot.xlabel( "position in read" )
      pyplot.ylabel( "base-call quality score" )

      pyplot.subplot( 224 )
      pyplot.pcolor( qual_arr_A_n.T ** gamma, cmap=pyplot.cm.Greens,
	   norm=pyplot.normalize( 0, 1 ) )
      pyplot.axis( (0, readlen-1, 0, max_qual+1 ) )
      pyplot.xlabel( "position in read" )

   else:

      pyplot.subplot( 211 )
      plot_bases( base_arr_U_n )
      pyplot.ylabel( "proportion of base" )
      pyplot.title( "%.3f million reads" % ( nreads_U / 1e6 ) )

      pyplot.subplot( 212 )
      pyplot.pcolor( qual_arr_U_n.T ** gamma, cmap=pyplot.cm.Greens,
	  norm=pyplot.normalize( 0, 1 ) )
      pyplot.axis( (0, readlen-1, 0, max_qual+1 ) )
      pyplot.xlabel( "position in read" )
      pyplot.ylabel( "base-call quality score" )


   pyplot.savefig( outfilename )

if __name__ == "__main__":
   main()
