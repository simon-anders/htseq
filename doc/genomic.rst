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
   .. class:: GenomicInterval( chrom, start, end, strand )

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
      
   .. attribute:: GenomicInterval.start_as_pos
                  GenomicInterval.end_as_pos
                  GenomicInterval.start_d_as_pos
                  GenomicInterval.end_d_as_pos
                  
      These attributes return :class:`GenomicPosition` objects referring to the
      respective positions.

Directional instantiation   
   .. function:: GenomicInterval_from_directional( chrom, start_d, length, strand="." )
   
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
   
      **pos** is an alias for :attr:`GenomicInterval.start_d`.
      
   All other attributes of :class:`GenomicInterval` are still exposed. Refrain from 
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
   .. class:: StepVector.StepVector( length = sys.maxint, typecode = 'd', start_index = 0 )

      Create a ``StepVector`` of the given ``length``, with indices starting
      at the given ``start_index`` and counting up to (but not including)
      ``start_index + length``.
      
      Note that the ``StepVector`` class resides in its own sub-module called 
      ``StepVector`` as well. Hence, the class's full name is ``HTSeq.StepVector.StepVector``.
      
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
- ``StepVector`` objects can be compared for equality and inequality. The test checks
  whether the elements are logically equal, i.e., it will not get confused
  if, say, one ``StepVector`` has a step from 10 to 20 with value 5, and the
  other has to steps, one from 10 to 15, and the other from 15 to 20, but
  both with the value 5.
   
Limitations

- ``StepVector`` does not realize when ``apply`` or ``add_value`` causes to
  adjacent steps to have the same value. If this situation arises when
  setting a value with indexing (``sv[3:5]``), it (usually but not always)
  merges the adjacent steps.
- If many steps are very narrow, a ``StepVector`` object may be less effective
  than an ordinary ``array``.


``GenomicArray``
================

A ``GenomicArray`` is a collection of :class:`StepVector` objects, typically either one or two 
for each chromosome of a genome. It allows to access the data in these StepVectors 
transparently via :class:`GenomicInterval` objects.

Instantiation
   .. class:: GenomicArray( chroms, stranded=True, typecode='d' )

   Creates a ``GenomicArray``. 
   
   If ``chroms`` is a list of chromosome names,  two (or one, see below) ``StepVector`` 
   objects for each chromosome are created, with start index 0 and indefinite
   length. If ``chroms`` is a ``dict``, the keys are used for the chromosome names
   and the values should be the lengths of the chromosome, i.e., the StepVectors
   index ranges are then from 0 to these lengths. (Note that the term chromosome
   is used only for convenience. Of course, you can as well specify contig IDs
   or the like.) Finally, if ``chroms`` is the string ``"auto"``, the GenomicArray
   is created without any chromosomes but whenever the user attempts to assign a 
   value to a yet unknown chromosome, a new one is automatically created with 
   :method:`GenomicArray.add_chrom`.
   
   If ``stranded`` is ``True``, two ``StepVector`` objects are created for each chromosome,
   one for the '+' and one for the '-' strand. For ``stranded == False``, only one
   ``StepVector`` per chromosome is used. In that case, the strand argument of
   all ``GenomicInterval`` objects that are used later to specify regions in the
   ``GenomicArray`` are ignored.
   
   The ``typecode`` is as described in :class:`StepVector`.
   
Attributes
   .. attribute: GenomicArray.stranded
                 GenomicArray.typecode
                 
      see above
      
   .. attribute: GenomicArray.step_vectors   
   
      a dict of (or a dict of dicts of) ``StepVector`` objects, using the chromosome
      names, and for stranded ``GenomicArray``s also ``'+'`` and ``'-'``, as keys::
      
      .. doctest::
      
         >>> ga = HTSeq.GenomicArray( [ "chr1", "chr2" ], stranded=False )
         >>> ga.step_vectors #doctest:+NORMALIZE_WHITESPACE
         {'chr2': <StepVector object, type 'd', index range 0:inf, 1 step(s)>,
          'chr1': <StepVector object, type 'd', index range 0:inf, 1 step(s)>}
         >>> ga = HTSeq.GenomicArray( [ "chr1", "chr2" ], stranded=True )
         >>> ga.step_vectors #doctest:+NORMALIZE_WHITESPACE
         {'chr2': {'+': <StepVector object, type 'd', index range 0:inf, 1 step(s)>,
                   '-': <StepVector object, type 'd', index range 0:inf, 1 step(s)>},
          'chr1': {'+': <StepVector object, type 'd', index range 0:inf, 1 step(s)>,
                   '-': <StepVector object, type 'd', index range 0:inf, 1 step(s)>}}

   .. attribute: GenomicArray.auto_add_chroms
   
      A boolean. This attribute is set to True if the GenomicArray was created with the ``"auto"``
      arguments for the ``chroms`` parameter. If it is true, an new chromosome
      will be added whenever needed.

                   
Data access
   One way to access your data is to access the ``step_vectors`` field directly.
   The more elegant way, however, is to use :class:GenomicInterval objects.
   
   To set an single position or an interval, use::
   
      >>> ga[ HTSeq.GenomicPosition( "chr1", 100, "+" ) ] = 7
      >>> ga[ HTSeq.GenomicInterval( "chr1", 250, 400, "+" ) ] = 20
      
   To read a single position::
   
      >>> ga[ HTSeq.GenomicPosition( "chr1", 300, "+" ) ]
      20.0
      
   To read an interval, use a ``GenomicInterval`` object as index, and
   obtain a ``StepVector`` with a sub-view::
   
      >>> sv = ga[ HTSeq.GenomicInterval( "chr1", 250, 450, "+" ) ]
      >>> sv
      <StepVector object, type 'd', index range 250:450, 2 step(s)>
      >>> list( sv.get_steps() )
      [(250, 400, 20.0), (400, 450, 0.0)]
      
   .. method:: GenomicArray.get_steps( iv = None, values_only = False )
      
   To get the steps, a method ``get_steps``, similar to :method:StepVector.get_steps,
   is supplied. It takes a ``GenomicInterval`` and returns an iterator generator of
   pairs. Each pair describes one step, first the range as ``GenomicInterval``,
   then the value::
   
      >>> list( ga.get_steps( HTSeq.GenomicInterval( "chr1", 0, 300, "+" ) ) ) #doctest:+NORMALIZE_WHITESPACE
      [(<GenomicInterval object 'chr1', [0,100), strand '+'>, 0.0),
       (<GenomicInterval object 'chr1', [100,101), strand '+'>, 7.0),
       (<GenomicInterval object 'chr1', [101,250), strand '+'>, 0.0),
       (<GenomicInterval object 'chr1', [250,300), strand '+'>, 20.0)]
         
   If ``iv`` is not specified, all steps are returned:
   
   .. doctest::
   
      >>> list( ga.get_steps(  ) )  #doctest:+ELLIPSIS,+NORMALIZE_WHITESPACE
      [..., [0,100), strand '+'>, 0.0), ..., 
       (<GenomicInterval object 'chr1', [0,9223372036854775807), strand '-'>, 0.0)]

   (Note that the long number here is ``sys.maxint``. This is because we have not
   specified an array length.)

   As in :class:StepVector.get_steps, setting ``values_only = True`` causes only the values being returned
         
Modifying values
   The two methods ``add_value`` and ``apply`` work as in :class:StepVector:
   
   .. method:: GenomicArray.add_value( value, iv )
   
      Add ``value`` to the elements specified by the ``GenomicInterval`` (or
      ``GenomicPosition``) object ``iv``. This only works for data type, for
      which the ``+`` operation is defined, i.e., numerical data.
      
      ::
      
         >>> ga.add_value( 10, HTSeq.GenomicInterval( "chr1", 300, 400, "+" ) )
   
   .. method:: GenomicArray.apply( func, iv )
   
      This method works the same way as :method:StepVector.apply, only now the range
      is specified as a ``GenomicInterval`` object. See :method:StepVector.apply for
      caveats.
   
Writing to a file
   .. method:: GenomicArray.write_bedgraph_file( file_or_filename, strand=".", track_options="" )
   
      Write out the data in the GenomicArray as a BedGraph_ track. This is a subtype of the Wiggle_ format
      (i.e., the file extension is usually ".wig") and such files can be conveniently viewed in
      a genome browser, e.g., with IGB_. 
      
      This works only for numerical data, i.e., ``datatype`` ``'i'`` or ``'d'``. As a bedgraph track
      cannot store strand information, you have to specify either ``'+'`` or ``'-'`` as the
      strand argument if your ``GenomicArray`` is stranded (``stranded==True``). Typically, you
      will write two wiggle files, one for each strand, and display them together.
      
Adding a chromosome
   .. method:: GenomicArray.add_chrom( chrom, length=sys.maxint, start_index=0 )
   
      Adds step vector(s) for a further chromosome. This is useful if you do not have a full
      list of chromosome names yet when instantiating the ``GenomicArray``.
      
.. _BedGraph: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
.. _Wiggle: http://genome.ucsc.edu/goldenPath/help/wiggle.html
.. _IGB: http://igb.bioviz.org/
   


``GenomicArrayOfSets``
======================

A ``GenomicArrayOfSets`` is a sub-class of :class:`GenomicArray` that deal with the common
special case of overlapping features. This is best explained by an example: Let's say, we
have two features, ``"geneA"`` and ``"geneB"``, that are at overlapping positions::

   >>> ivA = HTSeq.GenomicInterval( "chr1", 100, 300, "." )
   >>> ivB = HTSeq.GenomicInterval( "chr1", 200, 400, "." )
   
In a ``GenomicArrayOfSets``, the value of each step is a ``set`` and so can hold more than
one object:

.. doctest::

   >>> gas = HTSeq.GenomicArrayOfSets( ["chr1", "chr2"], stranded=False )
   >>> gas.add_value( "gene A", ivA )
   >>> gas.add_value( "gene B", ivB )
   >>> list( gas.get_steps( HTSeq.GenomicInterval( "chr1", 0, 500, "." ) ) ) #doctest:+NORMALIZE_WHITESPACE
   [(<GenomicInterval object 'chr1', [0,100), strand '.'>, set([])),
    (<GenomicInterval object 'chr1', [100,200), strand '.'>, set(['gene A'])),
    (<GenomicInterval object 'chr1', [200,300), strand '.'>, set(['gene A', 'gene B'])),
    (<GenomicInterval object 'chr1', [300,400), strand '.'>, set(['gene B'])),
    (<GenomicInterval object 'chr1', [400,500), strand '.'>, set([]))]


.. class:: GenomicArrayOfSets( chroms, stranded = True )

   Instantiation is as in :class:`GenomicArray`, only that ``datatype`` is always ``'O'``.
   
.. method:: GenomicArray.add_value( value, iv )

   The ``add_value`` method is changed to now add the value to the set of the step, performing
   proper copying, if necessary.
   


