quick = True

import os, sys
if os.system( "make" ) != 0:
   sys.exit()

import unittest
import cPickle
import copy
import pprint
import md5
from itertools import *
import numpy
import time

import HTSeq

worm_gtf_file = "test_data/Caenorhabditis_elegans.WS200.55.gtf.gz"
worm_chromlen_file = "test_data/Caenorhabditis_elegans.WS200.55.chromlengths" 
worm_fasta_file = "test_data/Caenorhabditis_elegans.WS200.55.dna.toplevel.fa.gz"
fastq_file = "test_data/SiBL1_27FEB09_s_5_sequence_head10000.txt"
hsapiens_chromlen_file = "test_data/Homo_sapiens.GRCh37.55.without_haplo_variants.chromlens"
yeasthybrid_rnaseq_bwt_file = "test_data/yeast-hybrid-RNASeq.bwtout"


class MyTestCase( unittest.TestCase ):
   
   def checkPickleAndCopy( self, obj ):
      self.assertEqual( obj, copy.copy( obj ), "Copying failed" )
      self.assertEqual( obj, cPickle.loads( cPickle.dumps( obj ) ), "Pickling failed" )
      
   def checkPickle( self, obj ):
      self.assertEqual( obj, cPickle.loads( cPickle.dumps( obj ) ), "Pickling failed" )
      

class Strand_Test( MyTestCase ):
   
   def setUp( self ):
      self.splus  = HTSeq.Strand( "+" )
      self.sminus = HTSeq.Strand( "-" )
      self.snone  = HTSeq.Strand( "." )
      
      self.checkPickleAndCopy( HTSeq.Strand( "+" ) )
      
   def test_str( self ):
      self.assertEqual( str( self.splus  ), "+" )
      self.assertEqual( str( self.sminus ), "-" )
      self.assertEqual( str( self.snone  ), "." )
      
   def test_cmp( self ):
      self.assert_( self.splus == HTSeq.Strand("+") )
      self.assert_( self.sminus == HTSeq.Strand("-") )
      self.assert_( self.snone == HTSeq.Strand(".") )
      self.assert_( self.sminus != HTSeq.Strand(".") )
      self.assert_( self.splus == "+" )
      self.assert_( self.sminus == "-" ) 
      self.assert_( self.snone == "." ) 
      self.assert_( self.sminus != "." ) 
      self.assert_( self.sminus != "ABC" ) 
      self.assert_( self.splus != self.sminus )
      self.assert_( self.sminus != self.snone )
      self.assert_( self.snone != self.splus )
      self.assert_( self.snone != None )
      self.assert_( self.snone != 3 )
      self.assert_( not self.snone == None )
      self.assert_( not self.snone == 3 )
      self.assertRaises( NotImplementedError, 
         lambda: self.splus > self.sminus )
                     

class GenomicInterval_Test( MyTestCase ):
   
   def test_various( self ):
      gi = HTSeq.GenomicInterval( "chr10", 1234567, 1234999, "-" )
      self.assertEqual( str(gi), "chr10:[1234567,1234999)/-" )
      gi2 = HTSeq.GenomicInterval_from_directional( "chr10", gi.start_d, gi.length, "-" )
      self.assert_( gi == gi2 )
      gi2.strand = HTSeq.Strand( "+" )
      self.assert_( gi != gi2 )
      
      self.assertEqual( gi2.start_d, gi2.start )
      self.assertEqual( gi2.end_d, gi2.end )
      self.assertEqual( gi.start_d, gi.end-1 )
      self.assertEqual( gi.end_d, gi.start+1 )
      
      self.failIf( gi.contains( gi2 ) )
      gi2.strand = "-"
      self.assert_( gi.contains( gi2 ) )
      gi2.start += 10
      self.assert_( gi.contains( gi2 ) )
      self.failIf( gi2.contains( gi ) )
      self.assert_( gi2.overlaps( gi ) )
      
      self.checkPickleAndCopy( gi )
      
      
class GenomicPosition_Test( MyTestCase ):
   
   def test_various( self ):
      gp = HTSeq.GenomicPosition( "chr10", 1234567, "-" )
      self.assertEqual( str(gp), "chr10:1234567/-" )
      gp2 = copy.copy( gp )
      gp2.pos += 1
      self.assertEqual( str(gp2), "chr10:1234568/-" )
      self.failIf( gp == gp2 )
      self.checkPickleAndCopy( gp )


class GenomicFeature_Test( MyTestCase ):
   
   def test_worm_peek( self ):
      worm_features = HTSeq.GFF_Reader( worm_gtf_file ) 
      self.assertEqual( str( HTSeq.peek( worm_features )[ 0 ] ),
         "<GenomicFeature: exon 'cTel54X.1' at III: 2916 -> 2817 (strand '-')>" )
      self.assertEqual( set( ( f.type for f in islice( worm_features, 100 ) ) ),
         set( ['start_codon', 'exon', 'stop_codon', 'CDS'] ) )

class GenomicArray_Test( MyTestCase ):
   
   def test_worm_nonstranded( self ):
      ga = HTSeq.GenomicArray( HTSeq.read_chrom_lens( worm_chromlen_file ), False, 'O' )
      worm_features = HTSeq.GFF_Reader( worm_gtf_file )       
      for f in islice( worm_features, 30 ):
         ga[ f.iv ] = f
      self.assertEqual( ga[ HTSeq.GenomicPosition( "III", 9300, "." ) ].name,
          "H10E21.2" )
      self.checkPickle( ga )    
            
   def test_worm_stranded( self ):
      ga = HTSeq.GenomicArray( HTSeq.read_chrom_lens( worm_chromlen_file ), True, 'O' )
      worm_features = HTSeq.GFF_Reader( worm_gtf_file )       
      for f in islice( worm_features, 30 ):
         ga[ f.iv ] = f
      self.assertRaises( KeyError, lambda: ga[ HTSeq.GenomicPosition( "III", 9300, "." ) ] )
      self.assertEqual( ga[ HTSeq.GenomicPosition( "III", 9300, "+" ) ].name,
          "H10E21.2" )
      self.checkPickle( ga )    

class FastaReader_Test( MyTestCase ):
   
   def test_chromlens( self ):
      if quick:
         return
      HTSeq.read_chrom_lens( worm_chromlen_file )
      self.assertEqual(
         dict( ( ( seq.name, len(seq) ) for seq in HTSeq.FastaReader( worm_fasta_file ) ) ),
         HTSeq.read_chrom_lens( worm_chromlen_file ) )
      

   
class FastqReader_Test( MyTestCase ):
   
   def test_namehash( self ):
      reads = HTSeq.FastqReader( fastq_file )
      m = md5.new()
      for r in reads:
         m.update( r.name )
      self.assertEqual( m.hexdigest(), "b3f4a3a81a2b850fdc7c619025773366" )
   
   def test_basestats( self ):
      reads = HTSeq.FastqReader( fastq_file )
      first, read = HTSeq.peek( reads )
      read_length = len( first )
      stats1 = numpy.zeros( ( 36, 5 ), 'i' )
      rows = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4 }
      for read in reads:
         for i in xrange( read_length ):
            stats1[ i, rows[read.seq[i]] ] += 1
      stats2 = HTSeq.base_counts_by_position( reads )
      self.assert_( ( stats1 == stats2 ).all() )
      
   def no_test_basestats_with_plot( self ):
      from matplotlib import pyplot
      stats = HTSeq.base_counts_by_position( HTSeq.FastqReader( fastq_file ) )
      for base in HTSeq.base_to_row:
         pyplot.plot( stats[ : , HTSeq.base_to_row[ base ] ] )
      pyplot.legend( HTSeq.base_to_row.keys() )
      pyplot.show ()
   
def table( iterable ):
   res = {}
   for a in iterable:
      try:
         res[ a ] += 1
      except KeyError:
         res[ a ] = 1
   return res
   
   
class BowtieReader_Test( MyTestCase ):   
   
   def test_multialgn_yeasthybr( self ):
      bwtout = HTSeq.BowtieReader( yeasthybrid_rnaseq_bwt_file )
      self.assertEqual(
         table( ( len(ma) for ma in HTSeq.bundle_multiple_alignments( bwtout ) ) ),
          {1: 34, 2: 2896, 3: 20, 4: 868, 6: 107} )

   def test_coverage( self ):
      ga = HTSeq.GenomicArray( HTSeq.read_chrom_lens( hsapiens_chromlen_file ), False, typecode='i' )
      for r in HTSeq.BowtieReader( "test_data/human_RNASeq.bwtout" ):
         r.iv.chrom = r.iv.chrom.split(" ")[0]
         try:
            ga.add_value( 1, r.iv )
         except KeyError:
            pass
      self.assertEqual( ga[ HTSeq.GenomicPosition("2", 177066000) ], 4)
         
   def test_rna_contamin( self ):
      gff = HTSeq.GFF_Reader( "test_data/yeast_small_acs.gff" )
      bwtout = HTSeq.BowtieReader(yeasthybrid_rnaseq_bwt_file )
      chromlens =  dict( ( (f.name, f.iv.length) for f in gff if f.type == "chromosome" ) )
      ga = HTSeq.GenomicArray( chromlens, True, 'b' )
      for f in gff:
         if f.type == "rRNA":
            ga[ f.iv ] = True
      rRNA = 0
      total = 0
      for ma in HTSeq.bundle_multiple_alignments( bwtout ):
         found = False
         for a in ma:
            a.iv.chrom = a.iv.chrom.split(".")[1]
            if a.iv.chrom == "scplasm1":
               continue
            if any( v for v in ga.get_steps( a.iv, values_only=True ) ):
               rRNA += 1
               found = True
               break
         total += 1
      self.assertEqual( ( rRNA, total ), (506, 3925) )
   

class SAM_Reader_Test( MyTestCase ):   
   
   def test( self ):
      samfile = HTSeq.SAM_Reader( "test_data/tophat_test.sam" )
      for algnt in samfile:
         if algnt.iv.length != len( algnt.read ):
            print algnt, algnt.iv.length, len( algnt.read )

   
execfile("StepVector_test.py")
   
unittest.main()