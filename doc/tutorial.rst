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

The read quality 
