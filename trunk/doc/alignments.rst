.. _alignments:

***************
Read alignments
***************

Concepts
========

There are a large number of different tools to align short reads to a reference. Most of
them use their own output format, even though the `SAM format`_ seems to become the common 
standard now that many of the newer tools use it.

.. _`SAM format`: http://samtools.sourceforge.net/SAM1.pdf

HTSeq aims to offer a uniform way to analyse alignments from different tools. To this end,
for all supported alignment formats a parse class is offered that reads an alignment file
and generates an iterator over the individual alignment records. These are represented as
objects of a sub-class of :class:`Alignment` and hence all offer a common interface.

So, you can easily write code that should work for all aligner formats. As a simple example,
consider this function that counts the number of reads falling on each chromosome::

   >>> import collections
   ... def count_in_chroms( alignments ):
   ...     counts = collections.defaultdict( lambda: 0 )
   ...     for almnt in alignments:
   ...         if almnt.aligned:
   ...             counts[ almnt.iv.chrom ] += 1
   ...     return counts

If you have a SAM file (e.g., from BWA or BowTie), you can call it with::

   >>> count_in_chroms( HTSeq.SAM_Reader( "yeast_RNASeq_excerpt.sam" ) )
   defaultdict( ..., {'XVI': 1509, 'V': 999, ..., 'XV': 2133})

If, however, you have done your alignment with Eland from the SolexaPipeline, which
uses the "Solexa export" format, you can use the same function, only using :class:`SolexaExportReader` 
instead of :class:`SAM_Reader`:

.. doctest::

   >>> count_in_chroms( HTSeq.SolexaExportReader( "mydata_export.txt" ) ) #doctest:skip

Both class generate iterators of similar objects. On the other hand, some formats contain more information
and then the ``Alignment`` objects from these contain additional fields.


Parser classes
==============

Depending on the format of your alignment file, choose from the following parsers:

.. class:: HTSeq.BowtieReader( filename_or_sequence )
           HTSeq.SAM_Reader( filename_or_sequence )
           HTSeq.SolexaExportReader( filename_or_sequence )
          
   All of these are derived from :class:`FileOrSequence`. When asked for an iterator,
   they yield ``Alignment`` objects of types :class:`BowtieAlignment`, :class:`SAM_Alignment`,
   or :class:`SolexaExportAlignment`. See below for their properties.
   
   Adding support for a new format is very easy. Ask me if you need something and
   I can probably add it right-away.
   Alternatively, you can convert your format to the SAM format. The SAMtools_
   contain Perl skripts to convert nearly all common formats.
   
.. _SAMtools: http://samtools.sourceforge.net/

   Paired-end support will be added very soon. Contact me if you need it now.
   

             
``Alignment`` and ``AlignmentWithSequenceReversal``
===================================================

.. class:: HTSeq.Alignment( )

   This is the abstract abse class of all Alignment classes. Any class derived 
   from ``Alignment`` has at least the following attributes:
   
   .. attribute:: read  (:class:`SequenceWithQuality`)
   
      The read. See :class:`SequenceWithQuality` for the sub-attributes.
      
      Note that some aligners store the reverse complement of the read if it was
      aligned to the '-' strand. In this case, the parser revers-complements the read
      again, so that you can be sure that the read is always presented as it was sequenced
      (see also :class:`AlignmentWithSequenceReversal`).
   
   .. attribute:: aligned  (boolean)
   
      Some formats (e.g., those of Maq and Bowtie) contain only aligned
      reads (and the aligner collects the 
      unaligned reads in a seperate FASTQ file if requested). For these formats, ``aligned``
      is always ``True``. Other formats (e.g., SAM and Solexa Export) list all reads, including those which could
      not be aligned. In that case, check ``aligned`` to see whether the read has an
      alignment.
      
   .. attribute:: iv   (:class:`GenomicInterval` or ``None``)
   
      The genomic interval to which the read was aligned (or ``None`` if ``aligned=False``).
      See :class:`GenomicInterval` for the sub-attributes. Note that different formats
      have different conventions for genomic coordinates. The parser class takes care
      of normalizing this, so that ``iv`` always adheres to the conventions outlined
      in :class:GenomicInterval. Especially, all coordinates are counted from zero, not one.


.. class:: AlignmentWithSequenceReversal( read_as_aligned, iv )

      Some aligners store the reverse complement of the read if it was
      aligned to the '-' strand. For these aligners, the Alignment class is derived
      from ``AlignmentWithSequenceReversal``, which undoes the reverse-complement if necessary
      to ensure that the ``read`` attribute always presents the read in the ordder in which
      it was sequenced.
      
      To get better performance, this is done via lazy evaluation, i.e., the 
      reverse complement is only calculated when the ``read`` attribute is accessed 
      for the first time. The original read as read from the file is stored as well. You
      can access both with these attributes:
      
      .. attribute:: AlignmentWithSequenceReversal.read_as_aligned
      
         The read as it was found in the file.
         
      .. attribute:: AlignmentWithSequenceReversal.read_as_sequences
      
         The read as it was sequences, i.e., an alias for ``read``.
      


Format-specific Alignment classes
=================================

Note: All format-specific Alignment classes take a string as argument for their constructor. This
is a line from the alignment file describing the alignment and is passed in by the corresponding
``Reader`` object. As you do not create ``Alignment`` objects yourself but get them from the ``Reader``
object you typically never call the constructor yourself.

.. class:: BowtieAlignment( bowtie_line )

   ``BowtieAlignment`` objects contain all the attributes from :class:Alignment and 
   :class:AlignmentWithSequenceReversal, and, in addition, these:
   
   .. attribute:: BowtieAlignment.reserved
   
      The ``reserved`` field from the Bowtie output file as a string. See the Bowtie manual for its meaning.

   .. attribute:: BowtieAlignment.substitutions  (string)
   
      The substitutions string that describes mismatches in the format ``22:A>C, 25:C>T``
      to indicate a change from A to C in position 22 and from C to T in position 25.
      No further parsing for this is offered yet.
      
.. class:: SAM_Alignment( line )

   ``BowtieAlignment`` objects contain all the attributes from :class:Alignment and 
   :class:AlignmentWithSequenceReversal, and, in addition, these:
   
   .. attribute:: SAM_Alignment.aQual
   
      The alignment quality score (an int) in Phread style encoding.

   .. attribute:: SAM_Alignment.cigar
   
      A list of :class:CigarOperation objects, as parsed from the extended CIGAR string. See
      :class:CigarOperation for details.
      
   .. attribute:: SAM_Alignment._tags
   
      At the moment, the extra tags are just put into this string field. A parser will
      be added soon.
      
.. class:: SolexaExportAlignment( line )

   ``BowtieAlignment`` objects contain all the attributes from :class:Alignment and 
   :class:AlignmentWithSequenceReversal, and, in addition, these:
   
   .. attribute:: SolexaExportAlignment.passed_filter
   
      Whether the read passed the chastity filter. If ``passed_filter==False``, the
      ``aligned==False``.

   .. attribute:: SolexaExportAlignment.nomatch_code
   
      For ``aligned==False``, a code indicating why no match could be found. See the 
      description of the 11th column of the Solexa Export format in the SolexaPipeline
      manual for the meaning of the codes. For ``aligned==True``, ``nomatch_code==None``.
         


Multiple alignments
===================

.. function:: HTSeq.bundle_multiple_alignments( sequence_of_alignments )

Some alignment programs, e.g., Bowtie, can output multiple alignments,
i.e., the same read is reported consecutively with different alignments.
This function takes an iterator over alignments (as provided by one of the 
alignment Reader classes) and bundles consecutive alignments regarding the 
same read to a list of Alignment objects and returns an iterator over these.


CIGAR strings
=============

When reading in SAM files, the CIGAR string is parsed and stored as a list of
``CigarOperation`` objects. For example, assume, a 36 bp read has been aligned
to the '+' strand of chromosome 'chr3', extending to the right from position
1000, with the CIGAR string ``"20M6I10M"``. The function :function:parse_cigar
spells out what this means::

   >>> HTSeq.parse_cigar( "20M6I10M", 1000, "chr2", "+" )
   [< CigarOperation: 20 base(s) matched on ref iv chr2:[1000,1020)/+, query iv [0,20) >,
    < CigarOperation: 6 base(s) inserted on ref iv chr2:[1020,1020)/+, query iv [20,26) >,
    < CigarOperation: 10 base(s) matched on ref iv chr2:[1020,1030)/+, query iv [26,36) >]

We can see that the map includes an insert. Hence, the affected coordinates run from 1000
to 1030 on the reference (i.e., the chromosome) but only from 0 to 36 on the query (i.e., the read).

We can convenient access to the parsed data by looking at the attributes of the three ``CigarOperation``
objects in the list.

.. class:: HTSeq.CigarOperation( ... )

   The available attributes are:
   
   .. attribute:: CigarOperation.type
   
      The type of the operation. One of the letters M, I, D, N, S, H, or P. Use
      the dict :variable: to transform this to names::
      
          >>> HTSeq.cigar_operation_names
          {'D': 'deleted',
           'H': 'hard-clipped',
           'I': 'inserted',
           'M': 'matched',
           'N': 'skipped',
           'P': 'padded',
           'S': 'soft-clipped'}
           
   .. attribute:: CigarOperation.size
   
      The number of affected bases, an int.
      
   .. attribute:: CigarOperation.ref_iv
   
      A :class:`GenomicInterval` specifying the affected bases on the reference. In case
      of an insertion, this is a zero-length interval.
      
   .. attribute:: CigarOperation.query_from
                  CigarOperation.query_to
                  
      To ints, specifying the affected bases on the query (the read). In case of a
      deletion, ``query_from == query_to``.
   

