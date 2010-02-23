************************************************************
HTSeq: Analysing high-throughput sequencing data with Python
************************************************************

HTSeq is a Python package that provides infrastructure to process data
from high-throughput sequencing assays. 

Please see the chapter :ref:`tour` first for an overview on the kind of analysis 
you can do with HTSeq and the design of the package, and the look at the reference
documentation. 

For installation instructions, see :ref:`install`.

Documentation overview
======================

* :ref:`tour`

  The Tour shows you how to get started. It explains how to install HTSeq, and then
  demonstrates typical analysis steps with explicite examples. Read this first, and 
  then see the Reference for details.
  
* Reference documentation

  The various classes of `HTSeq` are described here.

  * :ref:`sequences` 
  
    In order to represent sequences and reads (i.e., sequences with base-call quality 
    information), the classes ``Sequence`` and ``SequenceWithQualities`` are used.
    The classes ``FastaReader`` and ``FastqReader`` allow to parse FASTA and FASTQ
    files.
  
  * :ref:`genomic`
  
    The classes ``GenomicInterval`` and ``GenomicPosition`` represent intervals and
    positions in a genome. The class `GenomicArray` is an all-purpose container
    with easy access vcia a genomic interval or position, and ``GenomicArrayOfSets``
    is a special case useful to deal with genomi features (such as genes, exons,
    etc.)
    
  * :ref:`alignments`
  
    To process the output from short read aligners in various formats (e.g., SAM),
    the classes described here are used, to represent output files and alignments,
    i.e., reads with their alignment information.

Detailed table of contents
==========================

.. toctree::
   :maxdepth: 2
   
   install.rst
   tour.rst
   sequences.rst
   genomic.rst
   alignments.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

