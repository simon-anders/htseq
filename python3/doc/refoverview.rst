.. _refoverview:

.. currentmodule:: HTSeq

******************
Reference overview
******************

This page offers a brief overview over all classes and functions offered by HTSeq.

Parser and record classes
=========================

For all supported file formats, parser classes (called ``Reader``) are provided. These classes all
instatiated by giving a file name or an open file or stream and the function as iterator generators,
i.e., the parser objects can be used, e.g., in a ``for`` loop to yield a sequence of objects, each
desribing one record. The table gives the parse class and the record class yielded. For details,
see the linked documentation 

For most formats, functionality for writing files of the format is provided. See the detailed documentation 
as these methods and classes have varying semantics.

============================  ========================  ===========================  ===============================================  ===============================================
File format                   typical content           Parser class for reading     Record class yielded by parser                   Method/class method for writing
============================  ========================  ===========================  ===============================================  ===============================================
FASTA                         DNA sequences             :class:`FastaReader`         :class:`Sequence`                                :meth:`Sequence.write_to_fasta_file`
FASTQ                         sequenced reads           :class:`FastqReader`         :class:`SequenceWithQualities`                   :meth:`SequenceWithQuality.write_to_fastq_file`
GFF (incl. GFF3 and GTF)      genomic annotation        :class:`GFF_Reader`          :class:`GenomicFeature`                          :meth:`GenomicFeature.get_gff_line`
BED                           score data or annotation  :class:`BED_Reader`          :class:`GenomicFeature`                          
Wiggle (incl. BedGraph)       score data on a genome    :class:`WiggleReader`        pair: ``(``:class:`GenomicInterval`, ``float)``  :meth:`GenomicArray.write_bedgraph_file`
SAM                           aligned reads             :class:`SAM_Reader`          :class:`SAM_Alignment`                           :meth:`SAM_Alignment.get_sam_line`
BAM                           aligned reads             :class:`BAM_Reader`          :class:`SAM_Alignment`                           :class:`BAM_Writer`
VCF                           variant calls             :class:`VCF_Reader`          :class:`VariantCall`
Bowtie (legacy format)        aligned reads             :class:`BowtieReader`        :class:`BowtieAlignment`
SolexaExport (legacy format)  aligned reads             :class:`SolexaExportReader`  :class:`SolexaExportAlignment`
============================  ========================  ===========================  ===============================================  ===============================================

Most parser classes are sub-classes of class :class:`FileOrSequence`, which users will, however, rarely use directly.


Specifying genomic positions and intervals
==========================================

The class :class:`GenomicInterval` specifies an interval on a chromosome (or contig or the like). It is defined by specifying the chromosome (or contig) name,
the start and the end and the strand. Convenience methods are offered for different ways of accessing this information, and for tetsing
for overlap between intervals. A :class:`GenomicPosition`, technically a GenomicInterval of length 1, refers to a single nucleotide or base-pair
position.  

Objects of these classes are used internally wherever intervals or positions are represented, especially in record classes and as 
index keys to genomic array.

See page :ref:`genomic` for details.

 
Genomic arrays
==============

The classes :class:`GenomicArray` and :class:`GenomicArrayOfSets` are container classes to store data associated with genomic positions or intervals.

See page :ref:`genomic` for details.


Special features for SAM/BAM files
==================================

The class :class:`CigarOperation` offers a convenient way to handle the information encoded in the CIGAR field of SAM files.

The functions :func:`pair_SAM_alignments` and :func:`pair_SAM_alignments_with_buffer` help to ``pair up'' the records
in a SAM file that describe a pair of alignments for mated reads from the same DNA fragment.

Similarly, the function :func:`bundle_multiple_alignments` bundles multiple alignment record pertaining to the same read or read pair.

See page :ref:`alignments` for details.
