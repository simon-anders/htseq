.. _otherparsers:

*************
Other parsers
*************

.. currentmodule:: HTSeq

.. doctest:: 
   :hide:

   >>> import HTSeq


``VCF_Reader`` and ``VariantCall``
==================================

VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.

There is an option whether to contain genotype information on samples for each position or not.

See the definitions at

.. `1000genomes Project VCF <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40>_`
.. `1000genomes Project VCF for structural variants <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20%28Variant%20Call%20Format%29%20version%204.0/encoding-structural-variants>_'

As usual, there is a parser class, called **VCF_Reader**, that can generate an
iterator of objects describing the structural variant calls. These objects are of type :class:`VariantCall`
and each describes one line of a VCF file. See below for an example.

.. class:: VCF_Reader( filename_or_sequence )

   As a subclass of :class:`FileOrSequence`, VCF_Reader can be initialized either
   with a file name or with an open file or another sequence of lines.
   
   When requesting an iterator, it generates objects of type :class:`VariantCall`.
      
      .. attribute:: VCF_Reader.metadata
      
         VCF_Reader skips all lines starting with a single '#' as this marks
         a comment. However, lines starying with '##' contain meta data (Information about filters, and the fields in the 'info'-column).
      
      .. function:: parse_meta( header_filename = None )
      
         The VCF_Reader normally does not parse the meta-information and also the :class:`VariantCall` does not contain unpacked metainformation. The function parse_meta reads the header information either from the attached :class:`FileOrSequence` or from a file connection being opened to a provided 'header-filename'. This is important if you want to access sample-specific information for the :class`VariantCall`s in your .vcf-file.

      .. function:: make_info_dict( )
      
         This function will parse the info string and create the attribute :attr:`infodict` which contains a dict 
         with key:value-pairs containig the type-information for each entry of the :class:`VariantCall`'s info field.
  
.. class:: VariantCall( line, nsamples = 0, sampleids=[]  )

   A VariantCall object always contains the following attributes:
   
      .. attribute:: VariantCall.alt
         
         The alternative base(s) of the :class:`VariantCall`. This is a list containing all called alternatives.
         
      .. attribute:: VariantCall.chrom
         
         The Chromosome on which the :class:`VariantCall` was called.
         
      .. attribute:: VariantCall.filter
         
         This specifies if the :class:`VariantCall` passed all the filters given in the .vcf-header (value=PASS) or 
         contains a list of filters that failed (the filter-id's are specified in the header also).
         
      .. attribute:: VariantCall.format
         
         Contains the format string specifying which per-sample information is stored
         in :attr:`VariantCall.samples`.
         
      .. attribute:: VariantCall.id
         
         The id of the :class:`VariantCall`, if it has been found in any database, for unknown variants this will be 
         ".".
         
      .. attribute:: VariantCall.info
         
         This will contain either the string version of the info field for this :class:`VariantCall` or a dict with the 
         parsed and processed info-string.
         
      .. attribute:: VariantCall.pos
         
         A :class:`HTSeq.GenomicPosition` that specifies the position of the :class:`VariantCall`.
         
      .. attribute:: VariantCall.qual
         
         The quality of the :class:`VariantCall`.
         
      .. attribute:: VariantCall.ref
         
         The reference base(s) of the :class:`VariantCall`.

      .. attribute:: VariantCall.samples
         
         A dict mapping sample-id's to subdicts which use the :attr:`VariantCall.format` as 
         keys to store the per-sample information.
         
      .. function:: VariantCall.unpack_info( infodict )
         
         This function parses the info-string and replaces it with a dict rperesentation if the infodict of the 
         originating VCF_Reader is provided.

Example Workflow for reading the dbSNP in VCF-format (obtained from `dbSNP <ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/v4.0/00-All.vcf.gz>_`):

.. doctest::

    >>> vcfr = HTSeq.VCF_Reader( "00-All.vcf.gz" ) #doctest: +SKIP
    >>> vcfr.parse_meta() #doctest: +SKIP
    >>> vcfr.make_info_dict() #doctest: +SKIP
    >>> for vc in vcfr: #doctest: +SKIP
    ...    print vc,
    1:10327:'T'->'C'
    1:10433:'A'->'AC'
    1:10439:'AC'->'A'
    1:10440:'C'->'A'

*FIXME* The example above is not run, as the example file is still missing!


Wiggle Reader
=============

The `Wiggle format`_ (file extension often ``.wig``) is a format to describe numeric scores assigned to base-pair positions on a genome.
The class :class:`WiggleReader` is parser for such files. 

.. _`Wiggle format`: http://genome.ucsc.edu/goldenPath/help/wiggle.html

.. class:: WiggleReader( filename_or_sequence, verbose=True )

   The class is instatiated with the file name of a Wiggle file, or a sequence of lines in Wiggle format. A ``WiggleReader`` 
   object generates an iterator, which yields pairs of the form ``(iv, score)``, where ``iv`` is a :class:`GenomicInterval` 
   object and ``score`` is a ``float`` with the score that the file assigns to the specified interval. If ``verbose`` is set to
   True, the user is alerted to skipped lines (comments or ``browser`` lines) by a message printed to the standard output.


BED Reader
==========

The `BED format`_  is a format originally used to describe gene models but is also commonly used to describe other genomic features.

.. _`BED format`: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

.. class:: BED_Reader( filename_or_sequence )

   The class is instatiated with the file name of a BED file, or a sequence of lines in BED format. A ``BED_Reader`` 
   object generates an iterator, which yields a :class:`GenomicFeature` object for each line in the BED file (except for
   lines starting with ``track``, whcih are skipped).

   The attributes of the yielded ``GenomicFeature`` objects are as follows:

   ``iv``
      a :class:`GenomicInterval` object with the coordinates as given by the 1st, 2nd, 3rd, and 6th column of the BED file. If the
      BED file has less than 6 columns, the strand is set to "``.``".

   ``name``
      the name of feature as given in the 4th column, or ``unnamed``, if the file has only three columns

   ``type``
      always the string ``BED line``

   ``score``
      a float with the score as given by the 5th column (or ``None`` if the BED file has less 5 columns).

   ``thick``
      a :class:`GenomicInterval` object containg the "thick" part of the feature, as specified by the 6th and 7th column, with chromosome
      and strand copied from ``iv`` (or ``None`` if the BED file has less 7 columns).

   ``itemRgb``
      a list of three ``int`` values, taken from the 8th column (``None`` if the BED file has less 8 columns). In a BED file, this triple
      is meant to specify the colour in which the feature should be drawn in a browser.