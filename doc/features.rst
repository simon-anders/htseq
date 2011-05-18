.. _features:

********
Features
********


.. currentmodule:: HTSeq

.. doctest:: 
   :hide:

   >>> import HTSeq

The easiest way to work with annotation is to use :class:`GenomicArray` with ``typecode=='O'``
or :class:`GenomicArrayOfSets`. If you have your annotation in a flat file, with each
line describing a feature and giving its coordinates, you can read in the file line for line,
parse it (see the standard Python module ``csv``), use the information on chromosome, start,
end and strand to create a :class:`GenomicInterval` object and then store the data from the line
in the genomic array at the place indicated by the genomic interval.

For example, if you have data in a tab-separated file as follows:

.. doctest::

   >>> for line in open( "feature_list.txt" ):  #doctest:+NORMALIZE_WHITESPACE
   ...     print line,
   chr2  100	300	+	"gene A"
   chr2	200	400	-	"gene B"
   chr3	150	270	+	"gene C"

Then, you could load this information as follows::

   >>> import csv
   >>> genes = HTSeq.GenomicArray( [ "chr1", "chr2", "chr3" ], typecode='O' )
   >>> for (chrom, start, end, strand, name) in \
   ...        csv.reader( open("feature_list.txt"), delimiter="\t" ):
   ...     iv = HTSeq.GenomicInterval( chrom, int(start), int(end), strand )
   ...     genes[ iv ] = name

Now, to see whether there is a feature at a given :class:`GenomicPosition`, you just query the
genomic array::

   >>> print genes[ HTSeq.GenomicPosition( "chr3", 100, "+" ) ]
   None
   >>> print genes[ HTSeq.GenomicPosition( "chr3", 200, "+" ) ]
   gene C

See :class:`GenomicArray` and :class:`GenomicArrayOfSets` for more sophisticated use.


``GFF_Reader`` and ``GenomicFeature``
=====================================

One of the most common format for annotation data is GFF_ (which includes GTF_ as
a sub-type). Hence, a parse for GFF files is included in HTSeq.

.. _GFF: http://www.sanger.ac.uk/resources/software/gff/spec.html
.. _GTF: http://mblab.wustl.edu/GTF22.html

As usual, there is a parser class, called **GFF_Reader**, that can generate an
iterator of objects describing the features. These objects are of type :class`GenomicFeature`
and each describes one line of a GFF file. See Section :ref:`tour` for an example.

.. class:: GFF_Reader( filename_or_sequence, end_included=False )

   As a subclass of :class:`FileOrSequence`, GFF_Reader can be initialized either
   with a file name or with an open file or another sequence of lines.
   
   When requesting an iterator, it generates objects of type :class:`GenomicFeature`.
   
   The GFF specification is unclear on whether the end coordinate marks the last
   base-pair of the feature (closed intervals, ``end_included=True``) or the one
   after (half-open intervals, ``end_included=False``). The default, False, is
   correct for Ensembl GTF files. If in doubt, look for a CDS or stop_codon
   feature in you GFF file. Its length should be divisible by 3. If "end-start"
   is divisible by 3, you need ``end_included=False``. If "end-start+1" is
   divisible by 3, you need ``end_included=True``. 
   
   GFF_Reader will convert the coordinates from GFF standard (1-based, end
   maybe included) to HTSeq standard (0-base, end not included) by subtracting
   1 from the start position, and, for ``end_included=True``, also subtract 1 from
   the end position.
   
      .. attribute:: GFF_Reader.metadata
      
         GFF_Reader skips all lines starting with a single '#' as this marks
         a comment. However, lines starying with '##' contain meta data (at least
         accoring to the Sanger Institute's version of the GFF standard.) Such meta
         data has the format ``##key value``. When a metadata line is encountered,
         it is added to the ``metadata`` dictionary.
         
  
.. class:: GenomicFeature( name, type_, interval )

   A GenomicFeature object always contains the following attributes:
   
      .. attribute:: GenomicFeature.name
      
         A name of ID for the feature. As the GFF format does not have a dedicated
         field for this, the value of the first attribute in the *attributes* column is
         assumed to be the name of ID.
         
      .. attribute:: GenomicFeature.type
      
         The type of the feature, i.e., a string like ``"exon"`` or ``"gene"``. For GFF
         files, the 3rd column (*feature*) is taken as the type.
         
      .. attribute:: GenomicFeature.interval
      
         The interval that the feature covers on the genome. For GFF files, this information is taken
         from the first (*seqname*), the forth (*start*), the fifth (*end*), and the seventh (*strand*)
         column.
         
   When created by a :class:`GFF_Reader` object, the following attributes are also present, with the information
   from the remaining GFF columns:
   
      .. attribute:: GenomicFeature.source
      
         The 2nd column, denoted *source* in the specification, and intended to specify the
         data source.
     
      .. attribute:: GenomicFeature.frame
      
         The 8th column (*frame*), giving the reading frame in case of a coding feature. Its value
         is an integer (0, 1, or 2), or the string ``'.'`` in case that a frame is not specified or would not make sense.
   
      .. attribute:: GenomicFeature.score
      
         The 6th column (*score*), giving some numerical score for the feature. Its value
         is a float, or ``'.'`` in case that a score is not specified or would not make sense
      
      .. attribute:: GenomicFeature.attr
      
         The last (9th) column of a GFF file contains *attributes*, i.e. a list of name/value pairs.
         These are transformed into a dict, such that, e.g., ``gf.attr['gene_id']`` gives the value
         of the attribute ``gene_id`` in the feature described by ``GenomicFeature`` object ``gf``.
         The parser for the attribute field is reasonably flexible to deal with format variations
         (it was never clearly established whetehr name and value should be sperarated by a colon or an
         equal sign, and whether quotes need to be used) and also does a URL style decoding, as is often
         required.

   In order to write a GFF file from a sequence of features, this method is provided:
   
   .. method GenomicFeature.get_gff_line( with_equal_sign=False )
   
      Returns a line to describe the feature in the GFF format. This works even if the optional 
      attributes given above are missing. Call this for each of your ``GenomicFeature`` objects
      and write the lines into a file to get a GFF file.
      
.. function:: parse_GFF_attribute_string( attrStr, extra_return_first_value=False )      

   This is the function that :class:`GFF_Reader` uses to parse the attribute column. (See :attr:`GenomicFeature.attr`.)
   It returns a dict, or, if requested, a pair of the dict and the first value.

``VCF_Reader`` and ``VariantCall``
=====================================

VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.

There is an option whether to contain genotype information on samples for each position or not.

See the definitions at

.. `1000genomes Project VCF <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40>_`
.. `1000genomes Project VCF for structural variants <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20%28Variant%20Call%20Format%29%20version%204.0/encoding-structural-variants>_'

As usual, there is a parser class, called **VCF_Reader**, that can generate an
iterator of objects describing the structural variant calls. These objects are of type :class`VariantCall`
and each describes one line of a VCF file. See Section :ref:`tour` for an example.

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
      
         This function will parse the info string and create the attribute :attribute:`infodict` which contains a dict 
         with key:value-pairs containig the type-information for each entry of the :class:`VariantCall`'s info field.
  
.. class:: VariantCall( line, nsamples = 0, sampleids=[]  )

   A VariantCall object always contains the following attributes:
   
      .. attribute:: VariantCall.alt
         
         The alternative base(s) of the :class:`VariantCall`'. This is a list containing all called alternatives.
         
      .. attribute:: VariantCall.chrom
         
         The Chromosome on which the :class:`VariantCall`' was called.
         
      .. attribute:: VariantCall.filter
         
         This specifies if the :class:`VariantCall`' passed all the filters given in the .vcf-header (value=PASS) or 
         contains a list of filters that failed (the filter-id's are specified in the header also).
         
      .. attribute:: VariantCall.format
         
         Contains the format string specifying which per-sample information is stored
         in :attribute:`VariantCall.samples`.
         
      .. attribute:: VariantCall.id
         
         The id of the :class:`VariantCall`', if it has been found in any database, for unknown variants this will be 
         ".".
         
      .. attribute:: VariantCall.info
         
         This will contain either the string version of the info field for this :class:`VariantCall`' or a dict with the 
         parsed and processed info-string.
         
      .. attribute:: VariantCall.pos
         
         A :class:`HTSeq.GenomicPosition` that specifies the position of the :class:`VariantCall`'.
         
      .. attribute:: VariantCall.qual
         
         The quality of the :class:`VariantCall`'.
         
      .. attribute:: VariantCall.ref
         
         The reference base(s) of the :class:`VariantCall`'.

      .. attribute:: VariantCall.samples
         
         A dict mapping sample-id's to subdicts which use the :attribute:`VariantCall.format` as 
         keys to store the per-sample information.
         
      .. function:: VariantCall.unpack_info( infodict )
         
         This function parses the info-string and replaces it with a dict rperesentation if the infodict of the 
         originating VCF_Reader is provided.

Example Workflow for reading the dbSNP in VCF-format (obtained from `dbSNP <ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/v4.0/00-All.vcf.gz>_`):

.. doctest::

    >>> import HTSeq #doctest:+NORMALIZE_WHITESPACE
    vcfr = HTSeq.VCF_Reader( "00-All.vcf.gz" )
    vcfr.parse_meta()
    vcfr.make_info_dict()
    for vc in vcfr:
    ...    print vc,
    1:10327:'T'->'C'
    1:10433:'A'->'AC'
    1:10439:'AC'->'A'
    1:10440:'C'->'A'

