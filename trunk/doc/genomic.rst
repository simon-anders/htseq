.. _genomic:

************************************
Genomic intervals and genomic arrays
************************************

.. currentmodule:: HTSeq

.. doctest:: 
   :hide:

   >>> import HTSeq


``GenomicInterval``
===================

A genomic interval is a consecutive stretch on a genomic sequence such as a chromosome.
It is represented by a ``GenomicInterval`` object.

Instantiation
   .. class:: HTSeq.GenomicInterval( chrom, start, end, strand )

      ``chrom`` (string)
         The name of a sequence (i.e., chromosome, contig, or the like). 
         
      ``start``  (int)
         The start of the interval. Even on the reverse strand,
         this is always the smaller of the two values 'start' and 'end'.
         Note that all positions should be given and interpreted as 0-based value!
         
      ``end``  (int)
         The end of the interval. Following Python convention for 
         ranges, this in one more than the coordinate of the last base
         that is considered part of the sequence.
         
      ``strand``  (string)
         The strand, as a single character, ``'+'``, ``'-'``, or ``'.'``.
         ``'.'`` indicates that the strand is irrelevant.


Representation and string conversion
   The class's ``__str__`` method gives a spcae-saving description of the
   interval, the ``__repr__`` method is a bit more verbose::
   
      >>> iv = HTSeq.GenomicInterval( "chr3", 123203, 127245, "+" )
      >>> print iv
      chr3:[123203,127245)/+
      >>> iv
      <GenomicInterval object 'chr3', [123203,127245), strand '+'>

Attributes

   .. attribute:: GenomicInterval.chrom
                  GenomicInterval.start
                  GenomicInterval.end
                  GenomicInterval.strand
                 
      as above                 
                 
   .. attribute:: GenomicInterval.start_d

      The "directional start" position. This is the position of the
      first base of the interval, taking the strand into account. Hence, 
      this is the same as ``start`` except when ``strand == '-'``, in which 
      case it is ``end-1``.
      
      Note that if you set ``start_d``, both ``start`` and ``end`` are changed, 
      such that the interval gets the requested new directional start and its
      length stays unchanged.      
      
   .. attribute:: GenomicInterval.end_d

      The "directional end": The same as ``end``, unless ``strand=='-'``, 
      in which case it is ``start+1``. This convention allows to go from
      ``start_d`` to ``end_d`` (not including, as usual in Python, the last
      value) and get all bases in "reading" direction.
      
      ``end_d`` is not writable. 

   .. attribute:: GenomicInterval.length
   
      The length is calculated as end - start. If you set the length, 
      ``start_d`` will be preserved, i.e., ``end`` is changed, 
      unless the strand is ``-``, in which case ``start`` is changed.

Directional instantiation   
   .. function:: HTSeq.GenomicInterval_from_directional( chrom, start_d, length, strand="." )
   
      This function allows to create a new ``GenomicInterval`` object specifying
      directional start and length instead of start and end. 
   
Methods
   .. method:: GenomicInterval.is_contained_in( iv )
               GenomicInterval.contains( iv )
               GenomicInterval.overlaps( iv )
               
      These methods test whether the object is contained in, contains, or overlaps
      the second ``GenomicInterval`` object ``iv``. 
      
      For any of of these conditions
      to be true, the ``start`` and ``end`` values have to be appropriate, and furthermore,
      the ``chrom`` values have to be equal and the ``strand`` values consistent. The latter
      means that the strands have to be the same if both intervals have strand
      information. However, if at least one of the objects has ``strand == '.'``,
      the strand information of the other object is disregarded.
      
      Note that all three methods return ``True`` for identical intervals.
               
   .. method:: GenomicInterval.xrange( step = 1 )
               GenomicInterval.xrange_d( step = 1 )
               
      These methods yield iterators of :class:GenomicPosition objects from
      ``start`` to ``end`` (or, for ``xrange_d`` from ``start_d`` to ``end_d``).

   .. method:: GenomicInterval.extend_to_include( iv )
   
      Change the object's ``start`` end ``end`` values such that ``iv`` becomes contained.
         
Special methods      

   ``GenomicInterval`` implements the methods necessary for
   
   - obtaining a copy of the object (the ``copy`` method)
   - pickling the object
   - representing the object and converting it to a string (see above)
   - comparing two GenomicIntervals for equality and inequality
   - hashing the object


``GenomicPosition``
===================

A ``GenomicPosition`` represents the position of a single base or base pair, i.e., it is
an interval of length 1, and hence, the class is a subclass of :class:GenomicInterval.

.. class::  GenomicPosition( chrom, pos, strand='.' )

   The initialisation is as for a :class:GenomicInterval object, but no ``length`` argument is passed.
   
Attributes

   .. attribute:: pos
      ``pos`` is an alias for ``start_d``.
      
   All other attributes of ``GenomicInterval`` are still exposed. Refrain from 
   using them, unless you want to use the object as an interval, not as a position.
   Some of them are now read-only to prevent the length to be changed.
   
   
``StepVector``
==============


``GenomicArray``
================


``GenomicArrayOfSets``
======================
