.. _count:

.. program:: htseq-count


***********************************************
Counting reads in features with ``htseq-count``
***********************************************

Given a file with aligned sequencing reads and a list of genomic
features, a common task is to count how many reads map to each feature.

A feature is here an interval (i.e., a range of positions) on a chromosome
or a union of such intervals.

In the case of RNA-Seq, the features are typically genes, where each gene
is considered here as the union of all its exons. One may also consider
each exon as a feature, e.g., in order to check for alternative splicing.
For comparative ChIP-Seq, the features might be binding region from a 
pre-determined list.

Special care must be taken to decide how to deal with reads that overlap more
than one feature. The ``htseq-count`` script allows to choose between three
modes. Of course, if none of these fits your needs, you can write your own
script with HTSeq. See the chapter :ref:`tour` for a step-by-step guide on 
how to do so.

The three overlap resolution modes of ``htseq-count`` work as follows. For 
each position `i` in the read, a set `S(i)` is defined as the set of all 
features overlapping position `i`. Then, consider the set `S`, which is 
(with `i` running through all position within the read)

* the union of all the sets `S(i)` for mode ``union``.

* the intersection of all the sets `S(i)` for mode ``intersection-strict``.

* the intersection of all non-empty sets `S(i)` for mode ``intersection-nonempty``.

If `S` contains precisely one feature, the read is counted for this feature. If
it contains more than one feature, the read is counted as ``ambiguous`` (and
not counted for any features), and if ``S`` is empty, the read is counted
as ``no_feature``. 

The following figure illustrates the effect of these three modes:

.. image:: count_modes.png


Usage
-----

After you have installed HTSeq (see :ref:`install`), you can run ``htseq-count`` from
the command line::

   htseq-count [options] <sam_file> <gff_file>
   
If the file ``htseq-qa`` is not in your path, you can, alternatively, call the script with

::
   
   python -m HTSeq.scripts.count [options] <sam_file> <gff_file>
   
The <sam_file> contains the aligned reads in the SAM format. (Note that the 
SAMtools_ contain Perl scripts to convert most alignment formats to SAM.)
Make sure to use a splicing-aware aligner such as TopHat. HTSeq-count makes 
full use of the information in the CIGAR field.

To read from standard input, use ``-`` as ``<sam_file>``.

If you have paired-end data, 
you have to sort the SAM file by read name first.  (If your sorting tool 
cannot handle big files, try e.g. Ruan Jue's *msort*, available from the 
SOAP_ web site.)
         
.. _SAMtools: http://samtools.sourceforge.net/
.. _SOAP: http://soap.genomics.org.cn

The <gff_file> contains the features in the `GFF format`_.

.. _`GFF format`: http://www.sanger.ac.uk/resources/software/gff/spec.html

The script outputs a table with counts for each feature, followed by
the special counters, which count reads that were not counted for any feature
for various reasons, namely:

* ``no_feature``: reads which could not be assigned to any feature 
  (set `S` as described above was empty).
   
* ``ambiguous``: reads which could have been assigned to more than 
  one feature and hence were not counted for any of these (set `S`
  had mroe than one element).
  
* ``too_low_aQual``: reads which were not counted due to the ``-a``
  option, see below
  
* ``not_aligned``: reads in the SAM file without alignment

* ``alignment_not_unique``: reads with more than one reported alignment.
  These reads are recognized from the ``NH`` optional SAM field tag. 
  (If the aligner does not set this field, multiply aligned reads will 
  be counted multiple times.)
  

*Important:* The default for strandedness is *yes*. If your RNA-Seq data has not been made
with a strand-specific protocol, this causes half of the reads to be lost.
Hence, make sure to set the option ``--stranded=no`` unless you have strand-specific
data!
  
      
Options
.......


.. cmdoption::  -m <mode>, --mode=<mode>  

   Mode to handle reads overlapping more than one feature. Possible values for
   `<mode>` are ``union``, ``intersection-strict`` and ``intersection-nonempty``
   (default: ``union``)

.. cmdoption:: -s <yes, no or reverse>, --stranded=<yes, no, or reverse>

   whether the data is from a strand-specific assay (default: ``yes``)
   
   For stranded=no, a read is considered overlapping with a feature regardless
   of whether it is mapped to the same or the opposite strand as the feature.
   For stranded=yes and single-end reads, the read has to be mapped to the same
   strand as the feature. For paired-end reads, the first
   read has to be on the same strand and the second read on the opposite strand.
   For stranded=reverse, these rules are reversed.

.. cmdoption:: -a <minaqual>, --a=<minaqual>

   skip all reads with alignment quality lower than the given
   minimum value (default: 0)


.. cmdoption:: -t <feature type>, --type=<feature type>

   feature type (3rd column in GFF file) to be used, all
   features of other type are ignored (default, suitable
   for RNA-Seq and `Ensembl GTF`_ files: ``exon``)
   
.. _`Ensembl GTF`: http://mblab.wustl.edu/GTF22.html

.. cmdoption:: -i <id attribute>, --idattr=<id attribute>

   GFF attribute to be used as feature ID. Several GFF lines with the same
   feature ID will be considered as parts of the same feature. The feature ID
   is used to identity the counts in the output table. The default, suitable 
   for RNA-SEq and Ensembl GTF files, is ``gene_id``. 

.. cmdoption:: -o <samout>, --samout=<samout>

   write out all SAM alignment records into an output SAM
   file called <samout>, annotating each line with its
   assignment to a feature or a special counter
   (as an optional field with tag 'XF')
   
.. cmdoption:: -q, --quiet           
   
   suppress progress report and warnings

.. cmdoption:: -h, --help

   Show a usage summary and exit  
