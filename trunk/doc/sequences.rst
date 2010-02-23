.. _sequences:

*******************************
Sequences and FASTA/FASTQ files
*******************************
.. doctest:: 
   :hide:

   >>> import HTSeq

``Sequence``
============

.. class:: Sequence

A ``Sequence`` object holds a DNA sequence. Besides the actual sequence, an object
may also hold a name.

Instantiation
   .. function:: Sequence( self, seq, name="unnamed" )
      :noindex:

   Pass the DNA sequence and, optionally, a name or ID to the constructor::
   
      >>> myseq = HTSeq.Sequence( "ACCGTTAC", "my_sequence" )
   
   (If the name is omitted, the default ``"unnamed"`` is used.)
   
Attributes
   .. attribute:: Sequence.seq
                  Sequence.name
                  Sequence.descr
   
   The information can be accessed via the attributes ``seq`` and ``name``, which are strings::   
   
      >>> myseq.seq
      'ACCGTTAC'
      >>> myseq.name
      'my_sequence'

   There is a third attribute, called ``descr``, which is by default ``None`` but may contain 
   a "description". See class ``FastaReader`` for more information.
   
Representation and string conversion
   
   The ``__repr__`` method gives name and length::
   
      >>> myseq
      <_HTSeq.Sequence object 'my_sequence' (length 8)>

   The ``__str__`` method returns just the sequence:

      >>> print myseq
      ACCGTTAC
      
   Note that the length of a sequence is the number of bases::
   
      >>> len( myseq )
      8

Subsetting
   Subsetting works as with strings::

      >>> myseq2 = myseq[3:5]
      >>> myseq2.name
      'my_sequence[part]'
      >>> myseq2.seq
      'GT'
   
   (Note that ``"[part]"`` is appended to the name of the subsetted copy.)
   
Reverse complement   

   .. method:: Sequence.get_reverse_complement( )

::   

      >>> print myseq.get_reverse_complement()
      GTAACGGT
      >>> rc = myseq.get_reverse_complement()
      >>> rc.name
      'revcomp_of_my_sequence'
      >>> rc.seq
      'GTAACGGT'

Writing to FASTA file
   
   .. method: Sequence.write_to_fasta_file( fasta_file )
   
   To write ``Sequence`` objects into a FASTA file, open a text file for writing,
   then call ``def write_to_fasta_file`` for each sequence, providing the open
   file handle as only argument, and close the file::
   
      >>> my_fasta_file = open( "test.fa", "w" )
      >>> myseq.write_to_fasta_file( my_fasta_file )
      >>> my_fasta_file.close()
   
   To read from a FASTA file, see class ``FastaReader``.
   
Extended UIPAC letters
   These are not (yet) supported. A sequence should only contain A, C, G, and T.   
   

``SequenceWithQuality``
=======================   

.. class:: SequenceWithQuality

The sequences obtained from high-throughput sequencing devices (in the following also
referred to as "reads") typically come with `base-call quality scores`, which indicate
how sure the software was that the right base was called. The class ``SequenceWithQuality`` represents such reads. 

``SequenceWithQualities`` is a daughter class of ``Sequence`` and inherits all its features.

Instantiation

   .. function:: SequenceWithQuality( seq, name qualstr, qualscale="phred" ) 
      :noindex:

   A ``SequenceWithQualities`` can be instantiated as a ``Sequence``, but now with
   a third argument, the quality string::

      >>> myread = HTSeq.SequenceWithQualities( "ACGACTGACC", "my_read", "IICGAB##(!" )
   
   The quality string is interpreted as Sanger-encoded string of Phred values, as
   defined in the `FASTQ format specification`_, i.e., each letter in the quality
   string corresponds to one base in the sequence and if the value 33 is subtracted
   from the quality characters ASCII value, the Phred score is obtained.
   
   The Phred scores can then be found in the slot ``qual``::

      >>> myread.qualstr
      'IICGAB##(!'
      >>> myread.qual
      array([40, 40, 34, 38, 32, 33,  2,  2,  7,  0])
      
   If the quality string follows the `Solexa FASTQ` specification, the value to be
   subtracted is not 33 but 64. If you pass a quality string in this format, set
   ``qualscale="solexa"``.
   
   Prior to version 1.3, the SolexaPipeline software used a yet another style of encoding
   quality string. If you want to use this one, specify ``qualscale="solexa-old"``
   
.. _`FASTQ format specification`: http://maq.sourceforge.net/fastq.shtml

Attributes

   As for ``Sequence`` objects, there are attributes ``name``, ``seq``, and ``descr``.
   
   Furthermore, we now have the attributes ``qual`` and ``qualstr``, already mentioned
   above.
   
   .. attribute:: SequenceWithQuality.qual
   
      ``qual`` is a ``numpy`` array of data type *integer*, with as many elements
      as there are bases. Each element is a `Phred score`. A Phred score *S* is
      defined to mean that the base caller estimates the probability *p* of the
      base call being wrong as *p* = -log10 ( *S*/10 ).
      
      Note that ``qual`` is always the probability, even if the ``solexa-old`` quality
      string format has been used, which encodes the odds *p* ( 1 - *p* ), i.e., in that case,
      the odds are converted to probabilities.
      
   .. attribute:: SequenceWithQuality.qualstr
   
      The quality string according to Sanger Phred encoding. In case the quality was
      originally given in ``solexa`` or ``solexa-old`` format, it is converted::
      
         >>> read2 = HTSeq.SequenceWithQualities( "ACGACTGACC", "my_read", "hhgddaZVFF", "solexa" )
         >>> read2.qual
         array([ 40.,  40.,  39.,  36.,  36.,  33.,  26.,  22.,   6.,   6.])
         >>> read2.qualstr
         "IIHEEB;7''"
     
Writing to FASTQ file
   
   .. method:: SequenceWithQuality.write_to_fastq_file( fasta_file )
   
   To write ``SequenceWithQualities`` objects into a FASTQ file, open a text file for writing,
   then call ``write_to_fastq_file`` for each sequence, providing the open
   file handle as only argument, and close the file::
   
      >>> my_fastq_file = open( "test.fq", "w" )
      >>> myread.write_to_fastq_file( my_fastq_file )
      >>> my_fastq_file.close()
   
   Note that the reads will always be written with quality strings in Sanger encoding.
   
   To read from a FASTQ file, see class ``FastqReader``.
   

``FastaReader`` and ``FastqReader``
===========================

.. class:: FastaReader

   

