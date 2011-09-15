.. _history:

***************
Version history
***************

Version 0.5.3
=============

2011-06-29

- added the '--stranded=reverse' option to htseq-count


Patch versions:

* v0.5.3p1 (2011-07-15):

  - fix a bug in pair_sam_Alignment (many thanks for Justin Powell for
    finding the bug and suggesting a patch)
    
* v0.5.3p2 (2011-09-15)

  - fixed a bug (and a documentation bug) in trim_left/right_end_with_quals

* v0.5.3p3 (2011-09-15)

  - p2 was built improperly


Version 0.5.2
=============

2011-06-24

- added the '--maxqual' option to htseq-qa


Version 0.5.1
=============

2011-05-03

- added steps method to GenomicArray

Patch versions:

* v0.5.1p1 (2011-05-11):

  - fixed a bug in step_vector.h causing linkage failure under GCC 4.2

* v0.5.1p2 (2011-05-12):

  - fixed pickling

* v0.5.1p3 (2011-05-22):

  - fixed quality plot in htseq-qa (top pixel row, for quality score 40, was cut off)

Version 0.5.0
=============

2011-04-21

- refactoring of GenomicArray class:

  - field ``step_vectors`` replaced with ``chrom_vector``. These now contain
    dicts of dicts of ``ChromVector`` objects rather than ``StepVector`` ones.
    
  - ``chrom_vectors`` is now always a dict of dict, even for unstranded GenomicArrays
    to make it easier to loop over them. (The inner dict has either keys ``"+"``
    and ``"-"``, or just one key, ``"."``.)
    
  - The new ``ChromVector`` class wraps the actual vector and supports three different
    storage modes: ``step``, ``ndarray`` and ``memmap``, the latter two being numpy
    arrays, without and with memory mapping.
    
  - The ``GenomicArray`` constructor now take two new arguments, one for the storage
    class, one for the memmap directory (if needed).
    
  - The ``add_value`` methods had been replaced with an ``__iadd__`` method, to
    enable the ``+=`` semantics.
    
  - Similarily, ``+=`` for ``GenomicArrayOfSets`` adds an element to the sets.
  
  - Instead of ``get_steps``, now use ``steps``.
  
  
- new parser class ``VCF_Reader`` and record class ``VariantCall``

- new parser class ``BAM_Reader``, to add BAM support (including indexed random access)
  (requires PySam to be installed)

- new documentation page :ref:`tss`

- ``Fasta_Reader`` now allows indexed access to Fasta files (requires Pysam to be 
  installed)
  
- peek function removed  

Patch Versions:

- v0.5.0p1  (2011-04-22):

  - build was incomplete; fixed

- v0.5.0p2 (2011-04-22):

  - build was still faulty; new try

- v0.5.0p3 (2011-04-26)

  - fixed regression bug in htseq-count

Version 0.4.7
=============

2010-12-22

- added new option ``-o`` (or ``--samout``) to htseq-count

Patch versions:

* Version 0.4.7p1 (2011-02-14)

  - bug fix: GFF files with empty attribute fiels are now read correctly

* Version 0.4.7p2 (2011-03-13)

  - fixed assertion error in pair_SAM_alignment, triggered by incorrect flags

* Version 0.4.7p3 (2011-03-15)

  - fixed problem due to SAM_Alignment.peek (by removing the method)

* Version 0.4.7p4 (2011-03-18)

  - removed left-over debugging print statement


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
