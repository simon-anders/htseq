.. _history:

***************
Version history
***************

Version 0.4.3
=============

2010-05-01

New argument to constructer of GFF_Reader: ``end_include``

* Version 0.4.3-p1 (2010-05-04): version number was messed up; fixed

* Version 0.4.3-p2 (2010-05-15): fixed '-q' option in htseq-count

* Version 0.4.3-p3 (2010-05-15): parse_GFF_attribute_string can now deal with
  empty fields; score treated as float, not int

* Version 0.4.3-p3 (2010-05-15): 
  - parse_GFF_attribute_string can now deal with empty fields; 
    score treated as float, not int
  - fixed bug in SAM_Reader: can now deal with SAM files with 11 columns
  - SAM_Alignment._tags is now a list of strings

* Version 0.4.3-p4 (2010-05-16):
  bumped version number again just to make sure

Version 0.4.2
=============

2010-04-19

Bug fixes to htseq-count and pair_SAM_alignments. Bumped version number to avoid
confusion.

* Version 0.4.2-p1 (2010-04-20): there was still a bug left in htseq-count, fixed.

* Version 0.4.2-p2 (2010-04-26): bug fix: adapter trimming failed if the adapter
  was completely included in the sequence

* Version 0.4.2-p3

* Version 0.4.2-p4 (2010-04-29): bug fix: error in warning when htseq-count
  encountered an unknown chromosome 

* Version 0.4.2-p5 (2010-04-30): bug fixes: error in warning when PE positions
  are mismatched, and misleading error when calling get_steps with unstranded
  interval in a stranded array  


Version 0.4.1
=============

2010-04-19

Bug fixes:

* Fixed bug in ``htseq-count``: CIGAR strings with gaps were not correctly handled

* Fixed bug in Tour (last section, on counting): An wrong indent, and accidental
change to the ``exons`` variable invalidated data.

* SolexaExportReader no longer complains about multiplexing (indexing) not being supported.

* Mention link to example data in Tour.

* Fix installation instructions. (``--user`` does not work for Python 2.5.)

Enhancements:

* Paired-end support for SAM_Alignment.

* "_as_pos" attributes for GenomicInterval


Version 0.4.0
=============

2010-04-07

First "official" release, i.e., uploaded to PyPI and announced at SeqAnswers

Version 0.3.7
=============

2010-03-12

First version that was uploaded to PyPI
