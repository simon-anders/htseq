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

A ``StepVector`` is a container class to store a long one-dimensional array of values
which tend to be piece-wise constant. So instead of storing, say, 
``[ 3, 3, 3, 3, 3, 2, 2, 2 ]`` element for element, a step vector stores this
information in the form: one step from position 0 to 5 with the value 3, then another step 
from position 5 to 8 with the value 2. Each step consists of such a triple of start,
end, and value of the step (where, following Python convention, the end is not included
in the step), and a ``StepVector`` stores these steps in a way that is space efficient
and allows to quickly find all steps overlapping with a given query range.

A simple example: We create a ``StepVector`` of length 30, and then set 
positions 7 to 15 to the value 120.

   >>> sv = HTSeq.StepVector.StepVector( 30 )
   >>> sv[ 7:15 ] = 120
   >>> list( sv.get_steps() )
   [(0, 7, 0.0), (7, 15, 120.0), (15, 30, 0.0)]

Instantiation
   .. class:: HTSeq.StepVector.StepVector( length = sys.maxint, typecode = 'd', start_index = 0 )

      Create a ``StepVector`` of the given ``length``, with indices starting
      at the given ``start_index`` and counting up to (but not including)
      ``start_index + length``.
      
      Note that the ``StepVector`` class resides in its own sub-module called 
      ``Stepvector`` as well. Hence, the class's full name is ``HTSeq.StepVector.StepVector``.
      
      If no length is given, the ``StepVector`` extends indefinitely to the right,
      all the way to the largest representable `int` value (i.e., ``sys.maxint``).
      The ``start_index`` is zero by default. If you want your ``StepVector`` to
      extend indefinitely to the left as well, set ``start_index`` to the 
      ``-sys.maxint-1``, the smallest representable `int` value.
      
      The type codes are as in ``numpy``, i.e.:

      - ``'d'`` for float values (C type 'double'),
      - ``'i'`` for int values,
      - ``'b'`` for Boolean values,
      - ``'O'`` for arbitrary Python objects as value.
   
      The vector is initialized with the value zero (or, for typecode ``'O'``,
      with ``None``).

Attributes
   .. attribute:: StepVector.typecode
                  StepVector.start_index
                  
      see above

   .. attribute:: StepVector.num_steps
      
      the number of steps. (For a newly created object, this is 1.)
   
Getting values:
   Using a non-sliced array index returns a single value::
   
      >>> sv[12]
      120.0
      
   Using a slice produces a sub-vector::
   
      >>> sv[14:17]
      <StepVector object, type 'd', index range 14:17, 2 step(s)>
      
      Note that this sub-vector is a copy of the original vector. The copy is
      shallow, i.e., in case of type code ``'O'``, the elements still refer
      to the same objects.
      
      Slices with strides are not supported and raise an exception.
      
   Requesting an iterator from a StepVector generates an iterator that 
   goes through all elements:
   
   .. doctest::
   
      >>> iter(sv)              #doctest:+ELLIPSIS
      <generator object ...>
      >>> print list( sv )      #doctest:+NORMALIZE_WHITESPACE
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, 120.0, 120.0, 120.0, 120.0, 
      120.0, 120.0, 120.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      >>> for a in sv[5:8]:
      ...     print a
      0.0
      0.0
      120.0
      
   While the methods just describe access the data in its logical form, i.e.,
   as sequence of individual values, it is often desirable to access the steps, 
   which is done via the ``get_steps`` method.
   
   .. method:: StepVector.get_steps( start=None, stop=None, values_only = False )
   
      This method lists the steps as triples (start, end, value)::
      
         >>> list( sv.get_steps() )
         [(0, 7, 0.0), (7, 15, 120.0), (15, 30, 0.0)]
         >>> for (start, end, value) in sv.get_steps():
         ...     print "From %d to %d: %f" % ( start, end, value )
         From 0 to 7: 0.000000
         From 7 to 15: 120.000000
         From 15 to 30: 0.000000
         
      If ``start`` and/or ``stop`` is not specified, the start and end of the
      whole vector are assumed. If they are specified, the steps at the boundaries
      are reported as truncated to the extension of the query range::
   
         >>> list( sv.get_steps( 10, 20 ) )
         [(10, 15, 120.0), (15, 20, 0.0)]
         
      Often, the step boundaries are not of interest. Then, setting ``values_only=True``
      is helpful, as the method then returns an iterator over the covered values,
      instead of an iterator over the triples (start, end, value)::
      
         >>> list( sv.get_steps( 10, 20, True ) )
         [120.0, 0.0]
         
         
Setting values:
   To set a step, just use standard indexing:
   
      >>> sv[10] = 75
      >>> sv[13:15] = 20
      >>> list( sv[10:20] )
      [75.0, 120.0, 120.0, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            
Modifying values:
   In case of numeric data, the most common case is that on wishes to add a value.
   For technical reasons, the ``+=`` syntax is not supported. Instead, use the
   ``add_value`` method.
   
   .. method:: StepVector.add_value( value, start=None, stop=None )
   
   This method adds the value in the specified range, automatically introducing new step
   boundaries if necessary::
   
      >>> sv.add_value( 0.5, 10, 17 )
      >>> list( sv[10:20] )
      [75.5, 120.5, 120.5, 20.5, 20.5, 0.5, 0.5, 0.0, 0.0, 0.0]
      
   If ``start`` and ``stop`` are not specified, they default to start and end of the
   whole vector.
      
   If you need more complex manipulations, use ``apply``:
   
   .. method:: StepVector.apply( func, start = None, stop = None ):
   
   This method expects a function, which is called for the value of each step
   in the specified range::
   
      >>> list( sv[10:20] )
      [75.5, 120.5, 120.5, 20.5, 20.5, 0.5, 0.5, 0.0, 0.0, 0.0]
      >>> def times_seven( value ):
      ...     return value * 7
      >>> sv.apply( times_seven, 12, 20 )
      >>> list( sv[10:20] )
      [75.5, 120.5, 843.5, 143.5, 143.5, 3.5, 3.5, 0.0, 0.0, 0.0]

   This example could have been written easier with the ``lambda`` keyword. We undo the previous 
   operation in this example::
         
      >>> sv.apply( lambda x: x/7, 12, 20 )
      >>> list( sv[10:20] )
      [75.5, 120.5, 120.5, 20.5, 20.5, 0.5, 0.5, 0.0, 0.0, 0.0]
      
   Special care is needed when using ``apply`` for ``typecode`` ``'O'``. If the boundaries
   of the interval specified by ``start`` and ``stop`` do not both lie on
   step boundaries, one or both steps at the boundary are cut into two, both getting
   a reference to the same object as value. If the function then changes the object, this
   may affect the outer steps as well. The following example illustrates this::
   
      >>> sv2 = HTSeq.StepVector.StepVector( 6, 'O' )
      >>> sv2[ 2:5 ] = [ "foo", "bar" ]
      >>> list( sv2 )
      [None, None, ['foo', 'bar'], ['foo', 'bar'], ['foo', 'bar'], None]
      >>> def change_something( value ):
      ...     value[0] = "xxx"
      ...     return value
      >>> sv2.apply( change_something, 2, 4 )
      >>> list( sv2 )
      [None, None, ['xxx', 'bar'], ['xxx', 'bar'], ['xxx', 'bar'], None]
      
   Even though ``apply`` was supposed to only change the positions 2:4, the
   change affected positions 2:5. The mistake is that ``change_first`` did
   not make a copy of the object before modifying it. The following works as
   intended::
   
      >>> sv2 = HTSeq.StepVector.StepVector( 6, 'O' )
      >>> sv2[ 2:5 ] = [ "foo", "bar" ]
      >>> list( sv2 )
      [None, None, ['foo', 'bar'], ['foo', 'bar'], ['foo', 'bar'], None]
      >>> def change_something( value ):
      ...     new_value = value[:]   # make a copy of the list
      ...     new_value[0] = "xxx"
      ...     return new_value
      >>> sv2.apply( change_something, 2, 4 )
      >>> list( sv2 )
      [None, None, ['xxx', 'bar'], ['xxx', 'bar'], ['foo', 'bar'], None]
 
Special methods

- ``StepVector`` supports the protocols for copying and pickling.
- The length of a ``StepVector``, as reported by the ``len`` function, is 
  ``end - start``.
- ``StepVector``s can be compared for equality and inequality. The test checks
  whether the elements are logically equal, i.e., it will not get confused
  if, say, one ``StepVector`` has a step from 10 to 20 with value 5, and the
  other has to steps, one from 10 to 15, and the other from 15 to 20, but
  both with the value 5.
   

``GenomicArray``
================


``GenomicArrayOfSets``
======================
