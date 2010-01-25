HTSeq: A Python package to analyse high-throughput sequencing data
==================================================================

HTSeq is a python package that provides infrastructure to process data
from high-throughput sequencing assays. This tutorial demonstrates the
functionality of HTSeq by performing a number of common analysis tasks,
such as

- Getting statistical summaries about the base-call quality scores to
  study the data quality.
- Calculating a coverage vector and exporting it for visualization in
  a genome browser.
- Reading in annotation data from a GFF file.
- Assigning aligned reads from an RNA-Seq experiments to exons and
  genes.
  
This tutorial assumes that the reader is familiar with Python and with HTS
data.
  
Prequisites and installation
============================

To use HTSeq, you need at least version 2.5 of Python_ (Python 3 does not work yet), 
together with NumPy_,
a commonly used Python package for numerical calculations. Mac and Linux users 
will find that this is often pre-installed. (To check whether you have a working
NumPy installation, start Python and simply type ``import numpy``. If you do not
get an error, NumPy is available.) Especially Windows users might find the
`Enthought Python Distribution`_ convenient, which contains
Python, NyumPy and other add-ons and is available for all common platforms.

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _`Enthought Python Distribution`: http://www.enthought.com/products/epd.php

To install HTSeq, get a package tarball from [no web page yet, ask me by e-mail 
for one] and unpack it.

If you have a binary package (not containing a subdirectory called ``src``), open
a terminal and type, in the unpacked directory,
::

   python setup.py install
   
This will not work if you do not have write permission for Python's ``site-package``
directory. In that case, just find the directory ``HTSeq`` inside the unpacked
tarball and add the full path of the directory that contains the ``HTSeq`` directory
to Python's search path. This is done by either adding it to the ``PYTHONPATH``
environment variable or by putting, at the beginning of all your scripts or interactive
sessions::

   import sys
   sys.path.append ('/path/to/the/directory')

If you have a source package (conatining a subdirectory called ``src``), you need to 
first build the package before installing it. For this, make sure you have the GNU tool chain installed. 
(On a Mac, you need to install XCode (available from the Apple web site). On Ubuntu 
Linux, install the ``build-essential`` package. ) Then, open a terminal, ``cd`` 
into the unpacked directory, and type
::

   python setup.py build
  
Afterwards, proceed as with the binary package.

To test your installation, simply start Python and type ``import HTSeq``. No error 
message should appear.


Reading in reads
================

In the example data, a FASTQ file is provided with example reads from a yeast RNA-Seq
experiment. The file ``yeast_RNASeq_excerpt_sequence.txt`` is an excerpt of the
``_sequence.txt`` file produced by the SolexaPipeline software. We can access it from
HTSeq with
::

   >>> import HTSeq
   >>> fastq_file = HTSeq.FastqReader( "yeast_RNASeq_excerpt_sequence.txt", "solexa" )
  
The first argument is the file name, the optional second argument indicates that 
the quality values are encoded according to Solexa's specification. If you omit it,
the default "phred" is assumed, which means the encoding originally suggested
by the Sanger Institute. (A third option is "solexa_old", for data from the SolexaPipeline
prior to version 1.3.)

The variable ``fastq_file`` now refers to the file::

   >>> fastq_file
   <FastqReader object, connected to file name 'yeast_RNASeq_excerpt_sequence.txt'>
  
When used in a ``for`` loop, it generates an iterator of objects representing the
reads. Here, we use the ``islice`` function from ``itertools`` to cut after 10
reads::

   >>> import itertools
   >>> for read in itertools.islice( fastq_file, 10 ):
   ...    print read

   CTTACGTTTTCTGTATCAATACTCGATTTATCATCT
   AATTGGTTTCCCCGCCGAGACCGTACACTACCAGCC
   TTTGGACTTGATTGTTGACGCTATCAAGGCTGCTGG
   ATCTCATATACAATGTCTATCCCAGAAACTCAAAAA
   AAAGTTCGAATTAGGCCGTCAACCAGCCAACACCAA
   GGAGCAAATTGCCAACAAGGAAAGGCAATATAACGA
   AGACAAGCTGCTGCTTCTGTTGTTCCATCTGCTTCC
   AAGAGGTTTGAGATCTTTGACCACCGTCTGGGCTGA
   GTCATCACTATCAGAGAAGGTAGAACATTGGAAGAT
   ACTTTTAAAGATTGGCCAAGAATTGGGGATTGAAGA
   
Of course, there is more to a read than its sequence. The variable ``read`` still
contains the tenth read, and we may examine it::

   >>> read
   <SequenceWithQualities object 'HWI-EAS225:1:10:1284:142#0/1'>

A ``Sequence`` object has two slots, called ``seq`` and ``name``. This here is
a ``SequenceWithQualities``, and it also has a slot ``qual``::

   >>> read.name
   'HWI-EAS225:1:10:1284:142#0/1'
   >>> read.seq
   'ACTTTTAAAGATTGGCCAAGAATTGGGGATTGAAGA'
   >>> read.qual
   >>> read.qual
   array([ 33.,  33.,  33.,  33.,  33.,  33.,  29.,  27.,  29.,  32.,  29.,
           30.,  30.,  21.,  22.,  25.,  25.,  25.,  23.,  28.,  24.,  24.,
           29.,  29.,  29.,  25.,  28.,  24.,  24.,  26.,  25.,  25.,  24.,
           24.,  24.,  24.])

The values in the quality array are, for each base in the sequence, the Phred
score for the correctness of the base.

As a first simple example for the use of HTSeq, we now calculate the average
quality score for each position in the reads by adding up the ``qual`` arrays 
from all reads and the dividing by the number of reads. We sum everything up in
the variable ``qualsum``, a ``numpy`` array of integers::

   >>> import numpy
   >>> len( read )
   36
   >>> qualsum = numpy.zeros( len(read), numpy.int )

Then we loop through the fastq file, adding up the quality scores and
counting the reads::

   >>> nreads = 0
   >>> for read in fastq_file:
   ...    qualsum += read.qual
   ...    nreads += 1

The average qualities are hence::

   >>> qualsum / float(nreads)
   array([ 31.56970279,  30.08384335,  29.43881755,  29.00560022,
           28.55406216,  28.26937077,  28.46781871,  27.59194368,
           27.34221369,  27.57438298,  27.11888476,  27.19552782,
           26.84143366,  26.76363055,  26.450098  ,  26.7925917 ,
           26.43025721,  26.49973999,  26.13728549,  25.95939838,
           25.55042202,  26.20548822,  25.42457698,  25.72422897,
           25.04280171,  24.7525101 ,  24.48661946,  24.27181087,
           24.10840434,  23.68142726,  23.52150086,  23.49549982,
           23.11188448,  22.55858234,  22.43665747,  22.62470499])

If you have ``matplotlib`` installed, you can plot this::

   >>> from matplotlib import pyplot
   >>> pyplot.plot( qualsum / nreads )
   [<matplotlib.lines.Line2D object at 0x2ea8450>]
   >>> pyplot.show()

.. image:: qualplot.png

This is a very simple way of looking at the quality scores. For more sophisticated 
quality-control techniques, see [to be filled in].


What if you did not get the ``_sequence.txt`` file from your core facility but 
instead the ``export.txt`` file? While the former contains only the reads and
their quality, the latter also contains the alignment of the reads to a reference
as found by Eland. If we are only interested in the quality values, we simply
rewrite the commands as follows

   >>> solexa_export_file = HTSeq.SolexaExportReader( "yeast_RNASeq_excerpt_export.txt" )
   >>> nreads = 0
   >>> for aln in solexa_export_file:
   ...    qualsum += aln.read.qual
   ...    nreads += 1

We have simple replaced the ``FastqReader`` with a ``SolexaExportReader``, which 
iterates, when used in a ``for`` loop, over ``SolexaAlignment``objects. Each of
these contain a field ``read`` that contains the ``SequenceWithQualities``, as
before. There are more parses, for example the ``SAM_Parser`` that can read SAM
files, and generates ``SAM_Alignment`` objects. As all ``Alignment`` object
contain a ``read`` slot with the ``SequenceWithQualities``, we can use the same
code with any alignment filw for which a parser has been provided, and all we have
to change is the name of the reader class in the first line.

The other fields that all ``Alignment`` objects contain, is a Boolean called ``aligned``
that tells us whether the read has been aligned at all, and a field called ``iv``
(for "interval") that shows where teh read was aligned to. We use this information in
the next section.



Calculating coverage vectors
============================

By a "coverage vector", we mean a vector (one-dimensional array) of the length of
a chromosome, where each element counts how many reads cover the correspoding
base pair in their alignment. As chromosomes can be very long, it would be very 
inefficient to hold a coverage vector in memory by reserving space for each base
pair. Rather, we take advantage of the fact that the value of the coverage vector
usually stays constant (often it is just zero) over stretches of varying length,
which we call steps. A ``StepVector`` is a data structure defined for this purpose.

It works as follows: Let's define a ``Stepvector`` of length 30::

   >>> sv = HTSeq.StepVector.StepVector( 30 )
   
Initially, it has value 0 everywhere. We set the positions 7 to 15 to the value 120::

   >>> sv[ 7:15 ] = 120

Internally, ``sv`` now does not hold 30 numbers, but 3 steps, as follows::

   >>> list( sv.get_steps() )
   [(0, 7, 0.0), (7, 15, 120.0), (15, 30, 0.0)]

Each step is a triple, giving start, end and value of the step. If we now add the
value 100 to the positions 10 to 20, the steps get split accordingly::

   >>> sv.add_value( 100, 10, 20 )
   >>> list( sv.get_steps() )
   [(0, 7, 0.0), (7, 10, 120.0), (10, 15, 220.0), (15, 20, 100.0), (20, 30, 0.0)]
   
If you iterate over a ``StepVector``, it behaves like a list::

   >>> list( sv )
   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, 120.0, 120.0, 220.0, 220.0, 220.0, 220.0, 220.0, 100.0, 100.0, 
   100.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   
You can also take parts of a ``StepVector``, which produces a new, shorter, ``StepVector``::

   >>> sv[6:12]
   <StepVector object, type 'd', index range 6:12, 3 step(s)>
   >>> sv[6:12].get_steps()
   <generator object at 0x2ba95a8>
   >>> list( sv[6:12].get_steps() )
   [(6, 7, 0.0), (7, 10, 120.0), (10, 12, 220.0)]
   >>> list( sv[6:12] )
   [0.0, 120.0, 120.0, 120.0, 220.0, 220.0]


   
In practice, you will not work with ``StepVector``s directly, but rather with objects
of class ``GenomicArray``. These hold several step vectors, either one for each chromosome   
("non-stranded genomic array") or one for each strand, i.e., two per chromosome
("stranded genomic array"). To specify the locations of steps on a ``GenomicArray``, objects
of class ``GenomicInterval`` are used, which are instantiated by specifying chromsome
name, start, end, and position::

   >>> iv = HTSeq.GenomicInterval( "II", 100234, 100789, "+" )
   >>> iv
   <GenomicInterval object 'II', [100234,100789), strand '+'>
   >>> print iv
   II:[100234,100789)/+
   
A ``GenomicInterval`` has four slots which allow to access its data::
   
   >>> iv.chrom
   'II'
   >>> iv.start
   100234
   >>> iv.end
   100789
   >>> iv.strand
   <Strand '+'>
   
Two notes: ``chrom`` does not have to be chromosome, it could also be a contig name,
or any other identifier. ``strand`` can be ``+``, ``-``, or ``.``, where the latter
means "no strand", to be used whenever specifying a strand would be meaning-less.

A ``GenomicInterval`` has some more features, e.g., to calculate overlaps etc. See
[...] for these.


In order to calculate the coverage vectors for our yeast RNA-Seq data, we first need
to knwo the lengths of the chromosomes. One might take them from the lengths of
the reference FASTA files, but we simply specify the values here::

   >>> yeast_chrom_lengths = {
   ...    "2-micron": 6318,
   ...    "MT":    85779,
   ...    "I":    230208,
   ...    "II":   813178,
   ...    "III":  316617,
   ...    "IV":  1531918,
   ...    "V":    576869,
   ...    "VI":   270148,
   ...    "VII": 1090946,
   ...    "VIII": 562643,
   ...    "IX":   439885,
   ...    "X":    745745,
   ...    "XI":   666454,
   ...    "XII": 1078175,
   ...    "XIII": 924429,
   ...    "XIV":  784333,
   ...    "XV":  1091289,
   ...    "XVI":  948062
   ... }

Now, we define a ``GenomicArray``:

   >>> cvg = HTSeq.GenomicArray( yeast_chrom_lengths, stranded=True, typecode='i' )
   
As we specified ``stranded=True``, there are now two ``StepVector``s for each
chromosome, all holding integer values (``typecode='i'``):

   >>> import pprint
   >>> pprint.pprint( cvg.step_vectors )
   {'2-micron': {<Strand '+'>: <StepVector object, type 'i', index range 0:6318, 1 step(s)>,
                 <Strand '-'>: <StepVector object, type 'i', index range 0:6318, 1 step(s)>},
    'I': {<Strand '+'>: <StepVector object, type 'i', index range 0:230208, 1 step(s)>,
          <Strand '-'>: <StepVector object, type 'i', index range 0:230208, 1 step(s)>},
    'II': {<Strand '+'>: <StepVector object, type 'i', index range 0:813178, 1 step(s)>,
           <Strand '-'>: <StepVector object, type 'i', index range 0:813178, 1 step(s)>},
    'III': {<Strand '+'>: <StepVector object, type 'i', index range 0:316617, 1 step(s)>,
            <Strand '-'>: <StepVector object, type 'i', index range 0:316617, 1 step(s)>},
    'IV': {<Strand '+'>: <StepVector object, type 'i', index range 0:1531918, 1 step(s)>,
           <Strand '-'>: <StepVector object, type 'i', index range 0:1531918, 1 step(s)>},
    'IX': {<Strand '+'>: <StepVector object, type 'i', index range 0:439885, 1 step(s)>,
           <Strand '-'>: <StepVector object, type 'i', index range 0:439885, 1 step(s)>},
    'MT': {<Strand '+'>: <StepVector object, type 'i', index range 0:85779, 1 step(s)>,
           <Strand '-'>: <StepVector object, type 'i', index range 0:85779, 1 step(s)>},
    'V': {<Strand '+'>: <StepVector object, type 'i', index range 0:576869, 1 step(s)>,
          <Strand '-'>: <StepVector object, type 'i', index range 0:576869, 1 step(s)>},
    'VI': {<Strand '+'>: <StepVector object, type 'i', index range 0:270148, 1 step(s)>,
           <Strand '-'>: <StepVector object, type 'i', index range 0:270148, 1 step(s)>},
    'VII': {<Strand '+'>: <StepVector object, type 'i', index range 0:1090946, 1 step(s)>,
            <Strand '-'>: <StepVector object, type 'i', index range 0:1090946, 1 step(s)>},
    'VIII': {<Strand '+'>: <StepVector object, type 'i', index range 0:562643, 1 step(s)>,
             <Strand '-'>: <StepVector object, type 'i', index range 0:562643, 1 step(s)>},
    'X': {<Strand '+'>: <StepVector object, type 'i', index range 0:745745, 1 step(s)>,
          <Strand '-'>: <StepVector object, type 'i', index range 0:745745, 1 step(s)>},
    'XI': {<Strand '+'>: <StepVector object, type 'i', index range 0:666454, 1 step(s)>,
           <Strand '-'>: <StepVector object, type 'i', index range 0:666454, 1 step(s)>},
    'XII': {<Strand '+'>: <StepVector object, type 'i', index range 0:1078175, 1 step(s)>,
            <Strand '-'>: <StepVector object, type 'i', index range 0:1078175, 1 step(s)>},
    'XIII': {<Strand '+'>: <StepVector object, type 'i', index range 0:924429, 1 step(s)>,
             <Strand '-'>: <StepVector object, type 'i', index range 0:924429, 1 step(s)>},
    'XIV': {<Strand '+'>: <StepVector object, type 'i', index range 0:784333, 1 step(s)>,
            <Strand '-'>: <StepVector object, type 'i', index range 0:784333, 1 step(s)>},
    'XV': {<Strand '+'>: <StepVector object, type 'i', index range 0:1091289, 1 step(s)>,
           <Strand '-'>: <StepVector object, type 'i', index range 0:1091289, 1 step(s)>},
    'XVI': {<Strand '+'>: <StepVector object, type 'i', index range 0:948062, 1 step(s)>,
            <Strand '-'>: <StepVector object, type 'i', index range 0:948062, 1 step(s)>}}

The integer values are all initialized to 0. We may put them to a value, say 100,
at the genomic interval ``iv`` defined above:

    >>> cvg[ iv ] = 100
 
If we want to add a value, we use

   >>> cvg.add_value( 50, iv )
   
To see the effect, let's make the interval slightly longer and then look at the steps::

   >>> iv.start -= 30
   >>> iv.end += 100
   >>> pprint.pprint( list( cvg.get_steps( iv ) ) )
   [(<GenomicInterval object 'II', [100204,100234), strand '+'>, 0),
    (<GenomicInterval object 'II', [100234,100789), strand '+'>, 150),
    (<GenomicInterval object 'II', [100789,100889), strand '+'>, 0)]

With these tools, we can now calculate the coverage vector very easily. We just iterate
through all the reads and add the value 1 at the interval to which each read was aligned
to:

   >>> sam_file = HTSeq.SAM_Reader( "yeast_RNASeq_excerpt.sam" )
   >>> cvg = HTSeq.GenomicArray( yeast_chrom_lengths, stranded=True, typecode='i' )
   >>> for alngt in sam_file:
   ...    if alngt.aligned():
   ...       cvg.add_value( 1, alngt.iv )


