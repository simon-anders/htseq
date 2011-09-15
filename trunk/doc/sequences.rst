.. _sequences:

*******************************
Sequences and FASTA/FASTQ files
*******************************

.. currentmodule:: HTSeq

.. doctest:: 
   :hide:

   >>> import HTSeq

``Sequence``
============

A **Sequence** object holds a DNA sequence. Besides the actual sequence, an object
may also hold a name.

Instantiation
   .. class:: Sequence( seq, name="unnamed" )

   Pass the DNA sequence and, optionally, a name or ID to the constructor::
   
      >>> myseq = HTSeq.Sequence( "ACCGTTAC", "my_sequence" )
   
   (If the name is omitted, the default ``"unnamed"`` is used.)
   
Attributes
   .. attribute:: Sequence.seq
                  Sequence.name
                  Sequence.descr
   
   The information can be accessed via the attributes **seq** and **name**, which are strings::   
   
      >>> myseq.seq
      'ACCGTTAC'
      >>> myseq.name
      'my_sequence'

   There is a third attribute, called **descr**, which is by default ``None`` but may contain 
   a "description". See class :class:`FastaReader` for more information.
   
Representation and string conversion
   
   The **__repr__** method gives name and length::
   
      >>> myseq
      <Sequence object 'my_sequence' (length 8)>

   The **__str__** method returns just the sequence:

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
   
   To write **Sequence** objects into a FASTA file, open a text file for writing,
   then call **write_to_fasta_file** for each sequence, providing the open
   file handle as only argument, and close the file::
   
      >>> my_fasta_file = open( "test.fa", "w" )
      >>> myseq.write_to_fasta_file( my_fasta_file )
      >>> my_fasta_file.close()
   
   To read from a FASTA file, see class :class:`FastaReader`.
   
Extended UIPAC letters
   These are not (yet) supported. A sequence should only contain A, C, G, T
   and N.    
   
Counting bases

   .. method: Sequence.add_bases_to_count_array( count_array )
   
   For read quality assessment, it is often helpful to count the proportions
   of called bases, stratified by position in the read. To obtain such counts,
   the following idiom is helpful:
   
   .. doctest::
   
      >>> import numpy
      >>> reads = HTSeq.FastqReader( "yeast_RNASeq_excerpt_sequence.txt" )  
      >>> counts = numpy.zeros( ( 36, 5 ), numpy.int )
      >>> for read in reads:
      ...     read.add_bases_to_count_array( counts )
      >>> counts        #doctest:+NORMALIZE_WHITESPACE,+ELLIPSIS
      array([[16194,  2048,  4017,  2683,    57],
             [10716,  3321,  4933,  6029,     0],
             [ 7816,  5024,  5946,  6213,     0],
             ...
             [ 8526,  4812,  5460,  6197,     4],
             [ 8088,  4915,  5531,  6464,     1]])      
   
   Here, a two-dimensional numpy array of integer zeroes is defined and then
   passed to the **add_bases_to_count_array** method of each Sequence object obtained
   from the Fastq file. The method *add_bases_to_count_array* adds, for each base,
   a one to one of the array elements such that, in the end, the 36 rows of the array
   correspond to the positions in the reads (all of length 36 bp in this example), and
   the 5 columns correspond to the base letters 'A', 'C', 'G', 'T', and 'N', as given by
   the constant **base_to_columns**
   
   .. data:: base_to_column = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4 }
   
   Hence, we can get the proportion of 'C's in each position as follows:
   
   .. doctest::
   
      >>> counts = numpy.array( counts, numpy.float ) #doctest:+ELLIPSIS,+NORMALIZE_WHITESPACE
      >>> #counts[ : , HTSeq.base_to_column['C'] ] / counts.sum(1)
      array([ 0.08192328,  0.13284531,  0.20096804,  0.16872675,  0.21200848,
              ...
              0.18560742,  0.19236769,  0.19088764,  0.17872715,  0.1924877 ,
              0.19660786])
   
   (Here, we first convert the count array to type ``float`` to allow to proper
   division, and then divide the second column (``HTSeq.base_to_column['C']``) by
   the row-wise sums (``counts.sum(1)``; the ``1`` requests summing along rows).)

Trimming reads

   .. method:: Sequence.trim_left_end( pattern, mismatch_prop = 0. )
               Sequence.trim_right_end( pattern, mismatch_prop = 0. )
               
   In high-throughput sequencing, reads are sometimes contaminated with adapters
   or sequencing primers. These function take a pattern and attempt to match either
   the right end of the pattern to the left end of the sequence (``trim_left_end``)
   or the left end of the pattern to the right end of the sequence (``trim_right_end``).
   The match is the trimmed off.
   
   Here is an example::
   
      >>> seq2 = HTSeq.Sequence( "ACGTAAAGCGGTACGGGGGG" )
      >>> left_seq = HTSeq.Sequence( "CCCACG" )
      >>> print seq2.trim_left_end( left_seq )
      TAAAGCGGTACGGGGGG
      
   The right end of the pattern ("ACG") matched the left end of the sequence, and
   has hence been trimmed off.
   
   The optional argument ``mismatch_prop`` is the number of allowed mismatches as
   proportion of the length of the match::

      >>> right_seq = HTSeq.Sequence( "GGGTGGG" )
      >>> print seq2.trim_right_end( right_seq )
      ACGTAAAGCGGTACGGG
      >>> print seq2.trim_right_end( right_seq, 1/6. )
      ACGTAAAGCGGTAC
      >>> print seq2.trim_right_end( right_seq, 1/7. )
      ACGTAAAGCGGTACGGG
      
   Here, if we allow at least one mismatch per six bases, the whole pattern gets cut off.
   
   If you have quality information, you can use this, too, to specify the allowed amount
   of mismatch. See :meth:`SequenceWithQualities.trim_left_end_with_quals` and 
   :meth:`SequenceWithQualities.trim_left_end_with_quals`.


``SequenceWithQualities``
=========================   

The sequences obtained from high-throughput sequencing devices (in the following also
referred to as "reads") typically come with `base-call quality scores`, which indicate
how sure the software was that the right base was called. The class ``SequenceWithQualities`` represents such reads. 

``SequenceWithQualities`` is a daughter class of :class:`Sequence` and inherits all its features.

Instantiation

   .. class:: SequenceWithQualities( seq, name qualstr, qualscale="phred" ) 

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
         array([40, 40, 39, 36, 36, 33, 26, 22,  6,  6])
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
   
   To read from a FASTQ file, see class :class:`FastqReader`.

Counting quality values

   .. method: Sequence.add_qual_to_count_array( count_array )
   
   Similar to :meth:`Sequence.add_bases_to_count_array`, this method counts the
   occuring quality values stratified by position. This then allows to calculate
   average qualities as well as histograms.
   
   Here is a usage example:
   
   .. doctest::
   
      >>> import numpy       #doctest:+NORMALIZE_WHITESPACE,+ELLIPSIS
      >>> reads = HTSeq.FastqReader( "yeast_RNASeq_excerpt_sequence.txt", "solexa" )
      >>> counts = numpy.zeros( ( 36, 41 ), numpy.int )
      >>> for read in reads:
      ...    read.add_qual_to_count_array( counts )
      >>> #counts
      array([[   0,    0,   64, ...,    0,    0,    0],
             [   0,    0,   93, ...,    0,    0,    0],
             ..., 
             [   0,    0, 2445, ...,    0,    0,    0],
             [   0,    0, 2920, ...,    0,    0,    0]])

   The value ``counts[i,j]`` is then the number of reads for which the base at
   position ``i`` hat the  quality scores ``j``. According to the Fastq standard,
   quality scores range from 0 to 40; hence, the array is initialized to have 
   41 columns.
   
Trimming reads

   .. method:: SequenceWithQualities.trim_left_end_with_quals( pattern, max_mm_qual = 5 )
               SequenceWithQualities.trim_right_end_with_quals( pattern, max_mm_qual = 5 )
               
   These methods work as :meth:`Sequence.trim_left_end` and :meth:`Sequence.trim_right_end`
   (which are, of course, avilable for ``SequenceWithQualities`` objects, too). The difference
   is, that for the ``_with_quals`` trimming methods, the maximum amount of allowed mismatch is
   specified as the maximum value that the sum of the quality scores of the mismatched bases
   may take.
   
   *TODO*: Add example

``FastaReader``
===============

The FastaReader class allows to read and parse a FASTA file. It can generates an
iterator of ``Sequence`` objects.

.. class:: FastaReader( filename_or_sequence )

   As daughter class of ``FileOrSequence``, ``FastaReader`` can be instantiated
   with either a filename, or with a sequence. See :class:`FileOrSequence` for details.
   
Example 1
   The typical use for FastaReader is to go through a FASTA file and do something with
   each sequence, e.g.::
   
      >>> for s in HTSeq.FastaReader( "fastaEx.fa" ):
      ...     print "Sequence '%s' has length %d." % ( s.name, len(s) )
      Sequence 'sequence1' has length 72.
      Sequence 'sequence2' has length 70.
      
Example 2
   Often, one might to read a whole Fasta file into memory to access it as a dict.
   This is a good idiom for this purpose::
   
      >>> sequences = dict( (s.name, s) for s in HTSeq.FastaReader("fastaEx.fa") )
      >>> sequences["sequence1"].seq
      'AGTACGTAGTCGCTGCTGCTACGGGCGCTAGCTAGTACGTCACGACGTAGATGCTAGCTGACTAAACGATGC'


``FastqReader``
===============

The **FastqReader** class works similar to :class:`FastaReader`. It reads a Fastq file
and generates an iterator over :class:`SequenceWithQualities` objects.

.. class:: FastqReader( filename_or_sequence, qual_scale="phred" )

   As daughter class of ``FileOrSequence``, ``FastaReader`` can be instantiated
   with either a filename, or with a sequence. See :class:`FileOrSequence` for details.
   
   By default, the quality strings are assumed to be encoded according to the
   Sanger/Phred standard. You may alternatively specify ``"solexa"`` or ``"solexa-old"``
   (see :class:`SequenceWithQuality`).
