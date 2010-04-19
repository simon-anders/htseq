.. _history:

***************
Version history
***************

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
