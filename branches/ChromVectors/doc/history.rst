.. _history:

***************
Version history
***************

Version 0.4.7
=============

2010-12-22

- added new option ``-o`` (or ``--samout``) to htseq-count

Patch versions:

* Version 0.4.7p1 (2011-02-14)

  - bug fix: GFF files with empty attribute fiels are now read correctly

Version 0.4.6
=============

2010-12-09

- pair_SAM_alignments now handles multiple matches properly

- SAM_Alignments now allows access to optional fields via the new methods
  optional_field and optional_fields
  
- htseq-count now skips reads that are non-uniquely mapped according to the 'NH'
  optional field
  
- updated documentation    

Patch versions:

* Version 0.4.6p1 (2010-12-17)

  - updated htseq-count documentation page

  - htseq-count now accepts '-' as SAM file name

* Version 0.4.6p2 (2012-12-21)

  - corrected a bug in htseq-count regarding the handling of warnings and
    added SAM_Reader.peek.


Version 0.4.5
=============

2010-08-30

- correction to GenomicArray.get_steps() when called without arguments
- correction to FileOrSequence.get_line_number_string
- removed use of urllib's quote and unquote in GFF parsing/writing
- GFF_Reader now stores "meta information"
- qa.py now gives progress report
- auto add chrom now also works on read access
- refactored CIGAR parser
- added bool fields to SAM_Alignment for all flag bits

Patch versions:

* Version 0.4.5p1 (2010-10-08)

  - correction of a mistake in CIGAR checking, misreading symbol "N"

* Version 0.4.5p2 (2010-10-13)

  - Sequence.add_bases_to_count_array and hence htseq-qa now 
    accepts '.' instead of 'N' in a fastq file

* Version 0.4.5p3 (2010-10-20)

  - fixed error reporting for PE in htseq-count

* Version 0.4.5p4 (2010-10-21)

  - fixed another error reporting for PE in htseq-count

* Version 0.4.5p5 (2010-10-28)

  - Not only 'N' but also 'S' was read the wrong way. Fixed.
  
  - Cython had some odd way handling properties overloading attributes,
    which caused issues with 'Alignment.read'. Worked around.

* Version 0.4.5p6 (2010-11-02)

  - write_to_fastq should not break lines. Fixed.

* Version 0.4.5p7 (2010-11-16)

  - added fallback to distutils in case setuptools in unavailable
  
  - fixed documentation of '-a' option to htseq-count

Version 0.4.4
=============

2010-05-19

- StepVectors (and hence also GenomicArrays) now notice if, when setting the
  value of a step, this value is equal to an adjacent step and merge the steps.
  
- GenomicArray's constructor now allows the special value ``"auto"`` for its
  first arguments in order to start without chromosomes and automatically add
  them when first encountered.

Patch versions:

* Version 0.4.4p1 (2010-05-26):

  - minor change to make it run on Python 2.5 again
  - changed 'str' to 'bytes' at various places, now compiles with Cython 0.12
    (but no longer with Cython 0.11 and Python 2.5)

* Version 0.4.4p2 (2010-06-05):

  - change to SAM parser: if flag "query unmapped is set" but RNAME is not
    "*", a warning (rather than an error) is issued

* Version 0.4.4p3 (2010-06-25)

  - again removed an "except sth as e"

* Version 0.4.4p4 (2010-07-12)

  - dto.

* Version 0.4.4p5 (2010-07-13)

  - rebuilt with Cython 0.12.1 (previous one was accidently built with 
    Cython 0.11.1, causing it to fail with Python 2.5)

* Version 0.4.4p6 (2010-07-21)

  - fixed bug in error reporting in count.py
  - losened GFF attribute parsing
  - changed "mio" to "millions" in qa output
  - improved error reporting in GFF parser
  - made SAM parsing more tolerant


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
