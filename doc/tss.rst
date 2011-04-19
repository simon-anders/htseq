.. _tss:

******************************
A detailed use case: TSS plots
******************************

.. currentmodule:: HTSeq

A common task in ChIP-Seq analysis is to get profiles of marks with respect to
nearby features. For example, when analysing histone marks, one is often interested
in the position and extend of such marks in the vicinity of transcription start
sites (TSSs). To this end, one commonly calculates the coverage of reads or fragments
across the whole genome, then marks out fixed-size windows centered  around the 
TSSs of all genes, takes the coverages in these windows and adds them up to a
"profile" that has the size of the window. This is a simple operation, which, however,
can become challenging, when working with large genomes and very many reads. Here,
we demonstrate various ways of how data flow can be organized in HTSeq by means of
different solutions to this task.

As example data, we use a short excerpt from
the data set by Barski et al. (Cell, 2007, 129:823). We downloaded the FASTQ file
for one of the H3K4me3 samples (Short Read Archive accession number SRR001432),
aligned it against the Homo sapiens genome build GRCh37 with BWA, and provide the 
start of this file (actually only containing reads aligned to chromosome 1) as
file ``SRR001432_head.sam`` with the HTSeq example files. As annotation, we use
the file ``Homo_sapiens.GRCh37.56_chrom1.gtf``, which is the part of the Ensembl
GTF file for Homo sapiens for chromosome 1.

Using the full coverage
-----------------------

We start withy the straight-forward way of calculating the full coverage first
and then summing up the profile. This can be done as described in the Tour::

   >>> import HTSeq
   >>> coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
   >>> bamfile = HTSeq.BAM_Reader( "SRR001432_head.sam" )
   >>> gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.56_chrom1.gtf" )
   >>> for almnt in bamfile:
   ...    if almnt.aligned:
   ...       coverage[ almnt.iv ] += 1

To find the location of all transcription start sites, we can look in the GTF
file for exons with exon number 1 (as indicated by the ``exon_number``
attribute in Ensembl GTF files) and ask for their directional start (``start_d``).
The following loop extracts and prints this information (using ``itertools.islice``
to go through only the first 100 features in the GTF file)::

   >>> import itertools
   >>> for feature in itertools.islice( gtffile, 100):
   ...    if feature.type == "exon" and feature.attr["exon_number"] == "1":
   ...       print feature.attr["gene_id"], feature.attr["transcript_id"], feature.iv.start_d_as_pos
   ENSG00000223972 ENST00000456328 1:11873/+
   ENSG00000223972 ENST00000450305 1:12009/+
   ENSG00000227232 ENST00000423562 1:29369/-
   ENSG00000227232 ENST00000438504 1:29369/-
   ENSG00000227232 ENST00000488147 1:29569/-
   ENSG00000227232 ENST00000430492 1:29342/-
   ENSG00000243485 ENST00000473358 1:29553/+
   ENSG00000243485 ENST00000469289 1:30266/+
   ENSG00000221311 ENST00000408384 1:30365/+
   ENSG00000237613 ENST00000417324 1:36080/-
   ENSG00000237613 ENST00000461467 1:36072/-
   ENSG00000233004 ENST00000421949 1:53048/+
   ENSG00000240361 ENST00000492842 1:62947/+
   ENSG00000177693 ENST00000326183 1:69054/+

As the GTF file contains several transcripts for each gene, one TSS may appear 
multiple times, giving undue weight to it. Hence, we collect them in a ``set``
as this data type enforces uniqueness.

   >>> tsspos = set()
   >>> for feature in gtffile:
   ...    if feature.type == "exon" and feature.attr["exon_number"] == "1":
   ...       tsspos.add( feature.iv.start_d_as_pos )


Let's take one of these starting positions. To get a nice one, we manually chose
this one here, just for demonstration purposes::

   >>> p = HTSeq.GenomicPosition( "1", 37945599, "+" )

This is really one of the TSSs in the set:

   >>> p is tsspos
   True

We can get a window centered on this TSS by just subtracting and adding a fixed
value (half of the desired window size, let's use 3 kb)::
   
   >>> halfwinwidth = 3000
   >>> window = HTSeq.GenomicInterval( p.chrom, p.start - halfwinwidth, p.end + halfwinwidth, "." )
   >>> window
   <GenomicInterval object '1', [37942599,37948600), strand '.'>

We can check the coverage in this window by just subsetting and transforming to a list::

   >>> list( coverage[window] )  #doctest:+ELLIPSIS
   [0, 0, 0, ..., 0, 0]

As we will work with numpy from now on, it may be better to get this as 
numpy array::

   >>> import numpy
   >>> wincvg = numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
   array([0, 0, 0, ..., 0, 0, 0], dtype=int32)

With matplotlib, we can see that this vector is, in effect, not all zero:

   >>> from matplotlib import pyplot
   >>> pyplot.plot( wincvg )
   >>> pyplot.show()

We now initialize a numpy vector of the size of our window with zeroes:

   >>> profile = numpy.zeros( 2*halfwinwidth, dtype='i' )

Now, we can go through the TSS positions and add the coverage in their windows 
to get the profile:

   >>> for p in tsspos:
   ...    window = HTSeq.GenomicInterval( p.chrom, p.start - halfwinwidth, p.end + halfwinwidth, "." )
   ...    wincvg = numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
   ...    if p.strand == "+":
   ...       profile += wincvg
   ...    else:
   ...       profile += wincvg[::-1]

Note that we add the window coverage reversed ("`[::-1]`") if the gene was on the minus
strand.

Using matplotlib, we can plot this:

   >>> pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )
   >>> pyplot.show()

We can see clearly that the reads concentrate around the TSS, with a prominent peak 
a bit downstream (if you use matplotlib's interactive zoom, you can easily see that
the peak is at 153 bp) and a dip upstream, at -79 bp.

Going back to the beginning, there is one possible improvement: When calculating the
coverage, we just added one to all the base pairs that the read covered. However, the
fragment extends beyond the read, to a length of about 200 bp (the fragment size for
which Barski et al. selected). maybe we get a better picture. For this, we just add one
line, to extend the read to 200 bp. Using this, we now put the whole script together::

   import HTSeq
   import numpy
   from matplotlib import pyplot

   coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
   bamfile = HTSeq.BAM_Reader( "SRR001432_head.bam" )
   gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.56_chrom1.gtf" )
   halfwinwidth = 3000
   fragmentsize = 200

   for almnt in bamfile:
      if almnt.aligned:
         almnt.iv.length = fragmentsize
         coverage[ almnt.iv ] += 1

   tsspos = set()
   for feature in gtffile:
      if feature.type == "exon" and feature.attr["exon_number"] == "1":
         tsspos.add( feature.iv.start_d_as_pos )

   profile = numpy.zeros( 2*halfwinwidth, dtype='i' )      
   for p in tsspos:
      window = HTSeq.GenomicInterval( p.chrom, p.start - halfwinwidth, p.end + halfwinwidth, "." )
      wincvg = numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
      if p.strand == "+":
         profile += wincvg
      else:
         profile += wincvg[::-1]

   pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )
   pyplot.show()

The plot looks much smoother with the extended fragments.

The coverage vector can be held in memory, even for a very large genome, because large
parts of it are zero and even where there are reads, the values tend to stay constant 
for stretches of several bases. Hence, GenomicArray's step storage mode is useful.
If, however, extremely large amounts of reads are processed, the coverage vector can become "rough" and
change value at every position. Then, the step storage mode becomes inefficient
and we might be better off with an ordinary dense vector such as provided by numpy.
As this numpy vector becomes very large, it may not fit in memory, and the 'memmap'
storage (using numpy's memmap facility) then uses temporary files on disk. We
mention these possibilities as they may be useful when working with the ful coverage vector
is required. Here, however, we can do otherwise.

Using indexed BAM files
=======================

We do not need the coverage everythere. We only need it close to the TSSs. As we
have the data  
