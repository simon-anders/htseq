.. _overview:

.. currentmodule:: HTSeq

************************************************************
HTSeq: Analysing high-throughput sequencing data with Python
************************************************************

HTSeq is a Python package that provides infrastructure to process data
from high-throughput sequencing assays.

* Please see the chapter :ref:`tour` first for an overview on the kind of analysis 
  you can do with HTSeq and the design of the package, and then look at the reference
  documentation. 

..

* While the main purpose of HTSeq is to allow you to write your own analysis scripts,
  customized to your needs, there are also a couple of stand-alone scripts for
  common tasks that can be used without any Python knowledge. See the *Scripts* 
  section in the overview below for what is available.

..

* For downloads and installation instructions, see :ref:`install`.

Authors
=======

HTSeq is currently developed by `Fabio Zanini`_ at University of New South Wales in Sydney
(fabio *dot* zanini *at* unsw *dot* edu *dot* au) and `Simon Anders`_
(anders *at* embl *dot* de) at `EMBL Heidelberg`_ (`Genome Biology Unit`_).

.. _`Fabio Zanini`: http://fabilab.org
.. _`Simon Anders`: http://www.embl.de/research/units/genome_biology/huber/members/index.php?s_personId=6001
.. _`EMBL Heidelberg`: http://www.embl.de/
.. _`Genome Biology Unit`: http://www.embl.de/research/units/genome_biology/index.html


Citation
========

HTSeq is described in the following publications:

  | Fabio Zanini et al.
  | * HTSeq 2.0 -- Efficient manipulation of high-throughput sequencing data with long genomes
  | in preparation (2020)

and:

  | Simon Anders, Paul Theodor Pyl, Wolfgang Huber
  | *HTSeq --- A Python framework to work with high-throughput sequencing data*
  | Bioinformatics (2014), in print, online at `doi:10.1093/bioinformatics/btu638`_

.. _`doi:10.1093/bioinformatics/btu638`: https://doi.org/10.1093/bioinformatics/btu638

If you use HTSeq in research, please cite this paper in your publication.   


Documentation overview
======================

* :ref:`install`

  Download links and installation instructions can be found here

* :ref:`tour`

  The Tour shows you how to get started. It explains how to install HTSeq, and then
  demonstrates typical analysis steps with explicit examples. Read this first, and 
  then see the Reference for details.
  
* :ref:`tss`

  This chapter explains typical usage patterns for HTSeq by explaining in detail 
  three different solutions to the same programming task.

* :ref:`counting`

  This chapter explorer in detail the use case of counting the overlap of reads
  with annotation features and explains how to implement custom logic by
  writing on's own customized counting scripts


* Reference documentation

  The various classes of `HTSeq` are described here.

  * :ref:`refoverview` 

    A brief overview over all classes.

  * :ref:`sequences` 
  
    In order to represent sequences and reads (i.e., sequences with base-call quality 
    information), the classes :class:`Sequence` and :class:`SequenceWithQualities` are used.
    The classes :class:`FastaReader` and :class:`FastqReader` allow to parse FASTA and FASTQ
    files.
  
  * :ref:`genomic`
  
    The classes :class:`GenomicInterval` and :class:`GenomicPosition` represent intervals and
    positions in a genome. The class :class:`GenomicArray` is an all-purpose container
    with easy access via a genomic interval or position, and :class:`GenomicArrayOfSets`
    is a special case useful to deal with genomic features (such as genes, exons,
    etc.)
    
  * :ref:`alignments`
  
    To process the output from short read aligners in various formats (e.g., SAM),
    the classes described here are used, to represent output files and alignments,
    i.e., reads with their alignment information.

  * :ref:`features`
  
    The classes :class:`GenomicFeature` and :class:`GFF_Reader` help to deal with genomic
    annotation data.
    
  * :ref:`otherparsers`
  
    This page describes classes to parse VCF, Wiggle and BED files.

  * :ref:`misc`


* Scripts

  The following scripts can be used without any Python knowledge.
  
  * :ref:`qa`
  
    Given a FASTQ or SAM file, this script produces a PDF file with plots depicting
    the base calls and base-call qualities by position in the read. This is useful to
    assess the technical quality of a sequencing run.
  
  * :ref:`count`
  
    Given a SAM file with alignments and a GFF file with genomic features, this script
    counts how many reads map to each feature. 


* Appendices

..

  * :ref:`history`
  
..  
  
  * :ref:`contrib`

..

  * :ref:`Table of Contents<tableofcontents>`

..

  * :ref:`genindex`


..
   * :ref:`modindex`
   * :ref:`search`


License
=======

HTSeq is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The full text of the GNU General Public License, version 3, can be found
here: http://www.gnu.org/licenses/gpl-3.0-standalone.html
