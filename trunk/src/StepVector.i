%module StepVector

%include "exception.i"

%{ #define AUTOPYOBJPTR_EXTRAOPS %}
%include "AutoPyObjPtr.i"

%{
   #include "step_vector.h"

   class pystyle_stopiteration {};

   template< class T >
   class step_vector_pystyle_iterator
   {
      typename step_vector<T>::const_iterator current;
      typename step_vector<T>::const_iterator last;
     public: 
      step_vector_pystyle_iterator( typename step_vector<T>::const_iterator first,
         typename step_vector<T>::const_iterator last_ );
      std::pair< long int, T > next( );
      step_vector_pystyle_iterator<T> * __iter__( );
   };

   template< class T >
   class step_vector_for_python : public step_vector<T>
   {
     public:
      step_vector_for_python( long int length, 
         long int min_index_=0 );
      step_vector_pystyle_iterator<T> get_all_values_pystyle( ) const;
      step_vector_pystyle_iterator<T> get_values_pystyle( long int from ) const;
      int num_values( ) const;
   };

   template< class T >
   step_vector_pystyle_iterator<T>::step_vector_pystyle_iterator( 
         typename step_vector<T>::const_iterator first,
         typename step_vector<T>::const_iterator last_ )
    : current( first ), last( last_ ) 
   {
   }

   template< class T >
   step_vector_for_python<T>::step_vector_for_python( long int length, 
         long int min_index_ )
    : step_vector<T>( length, min_index_ )
   {}

   template< class T >
   step_vector_pystyle_iterator<T> step_vector_for_python<T>::get_all_values_pystyle( ) const
   {
      return step_vector_pystyle_iterator<T>( this->begin(), this->end() );
   }   

   template< class T >
   step_vector_pystyle_iterator<T> step_vector_for_python<T>::get_values_pystyle( long int from ) const
   {
      return step_vector_pystyle_iterator<T>( this->get_values( from ), this->end() );
   }   
   
   template< class T >
   int step_vector_for_python<T>::num_values( ) const
   {
      return this->m.size();
   }

   template< class T >
   std::pair< long int, T > step_vector_pystyle_iterator<T>::next( )
   {
      if( current == last )
         throw pystyle_stopiteration ();
      else {
         return *current++;
      }
   }   

   template< class T >   
   step_vector_pystyle_iterator<T> * step_vector_pystyle_iterator<T>::__iter__( )
   {
      return this;
   }   
               
%}

%exception next {
   try {
      $action
   } catch (pystyle_stopiteration &e) {
      PyErr_SetString( PyExc_StopIteration, "" );
      return NULL;
   }
}

template< class T1, class T2 >
class std::pair
{
  public:
   T1 first;
   T2 second;
   pair( T1 first_, T2 second_);
};

template< class T >
class step_vector_pystyle_iterator
{
  public: 
   step_vector_pystyle_iterator( typename step_vector<T>::const_iterator first,
      typename step_vector<T>::const_iterator last_ );
   std::pair< long int, T > next( );
   step_vector_pystyle_iterator<T> * __iter__( );   
};

#define TRY_CATCH_INDEX \
   try {                                                                   \
      $action                                                              \
   } catch (std::out_of_range &e) {                                        \
      SWIG_exception(SWIG_IndexError, e.what() );                          \
   } catch (type_error_non_arith &e) {                                               \
      SWIG_exception(SWIG_TypeError, "Illegal arithmetic operation" );     \
   }

%exception set_value { TRY_CATCH_INDEX }
%exception add_value { TRY_CATCH_INDEX }

template< class T >
class step_vector_for_python {
  public: 
   long int min_index;
   long int max_index;
   step_vector( long int length, long int min_index_=0 );
   void set_value( long int from, long int to, T value );
   void add_value( long int from, long int to, T value );
   step_vector_pystyle_iterator<T> get_all_values_pystyle( ) const;  
   step_vector_pystyle_iterator<T> get_values_pystyle( long int from ) const;
   int num_values( ) const;
};

%template( _Pair_int_float ) std::pair< long int, double >;
%template( _StepVector_Iterator_float ) step_vector_pystyle_iterator< double >; 
%template( _StepVector_float ) step_vector_for_python< double >; 

%template( _Pair_int_int ) std::pair< long int, int >;
%template( _StepVector_Iterator_int ) step_vector_pystyle_iterator< int >; 
%template( _StepVector_int ) step_vector_for_python< int >; 

%template( _Pair_int_bool ) std::pair< long int, bool >;
%template( _StepVector_Iterator_bool ) step_vector_pystyle_iterator< bool >; 
%template( _StepVector_bool ) step_vector_for_python< bool >; 

%template( _Pair_int_obj ) std::pair< long int, AutoPyObjPtr >;
%template( _StepVector_Iterator_obj ) step_vector_pystyle_iterator< AutoPyObjPtr >; 
%template( _StepVector_obj ) step_vector_for_python< AutoPyObjPtr >; 

%pythoncode %{

import sys

class StepVector( object ):

   """A step vector is a vector with integer indices that is able to store
   data efficiently if it is piece-wise constant, i.e., if the values change
   in "steps". So, if a number of adjacent vectort elements have the same
   value, this values will be stored only once.
   
   The data can be either one of a number of elementary types, or any object.
   
   Usage example:
   
   >>> sv = StepVector.StepVector( 20 )
   >>> sv[5:17] = 13
   >>> sv[12]
   13.0
   >>> list( sv )
   [0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 0.0, 0.0, 0.0]
   >>> list( sv.get_steps() )
   [(0, 5, 0.0), (5, 17, 13.0), (17, 20, 0.0)]

   """

   def __init__( self, length = sys.maxint, typecode = 'd', start_index = 0 ):
      """Construct a StepVector of the given length, with indices starting
      at the given start_index and counting up to (but not including)
      start_index + length.
      
      The typecode may be:
        'd' for float values (C type 'double'),
        'i' for int values,
        'b' for Boolean values,
        'O' for arbitrary Python objects as value.
   
      The vector is initialized with the value zero (or, for typecode 'O',
      with None).
      """
      if typecode == 'd':
         swigclass = _StepVector_float
      elif typecode == 'i':
         swigclass = _StepVector_int
      elif typecode == 'b':
         swigclass = _StepVector_bool
      elif typecode == 'O':
         swigclass = _StepVector_obj
      else:
         raise ValueError, "unsupported typecode"
      self._typecode = typecode
      self._swigobj = swigclass( length, start_index )    
   
   def __setitem__( self, index, value ):
      """To set element i of StepVector sv to the value v, write
         sv[i] = v
      If you want to set a whole step, say, all values from i to j (not
      including j), write
         sv[i:j] = v
      Note that the StepVector class will only notice that all the values
      from i to j are equal if you assign them in this fashion. Assigning each
      item individually in a loop from i to j will result in the value v being
      stored many times.
      """
      if isinstance( index, slice ):
         if index.step is not None and index.step != 1:
             raise ValueError, "Striding slices (i.e., step != 1) are not supported"
         start = index.start if index.start is not None else self.start_index()
         stop = index.stop if index.stop is not None else ( start + len(self) )
         self._swigobj.set_value( start, stop-1, value )
         # Note the "-1": The C++ object uses closed intervals, but we follow
         # Python convention here and use half-open ones.
      else:
         self._swigobj.set_value( index, index, value )
    
   def get_steps( self, start=None, stop=None, values_only = False ):
      """To get a succinct representation of the StepVector's content, call
      the 'get_steps' method. It returns an iterator that generates triples
      of values. Each triple contains one step, giving first the start index
      of the step, then the stop index (i.e., one more than the index of the 
      last element), and as third element the value of the step.
      
      If you want to see only a part of the StepVector, use the 'start' and
      'stop' parameters of 'get_steps' to specify a window.
      
      Sometimes, one might only be interested in the values, not the step
      boundaries. Then, set 'values_only' to true, and the iterator generates
      only the values insted of the triples.
      """
      if start is None:
         startvals = self._swigobj.get_all_values_pystyle( )
         pair = startvals.next()
         prevstart, prevval = pair.first, pair.second
      else:
         startvals = self._swigobj.get_values_pystyle( start )
         prevstart = start
         prevval = startvals.next().second
      for pair in startvals:
         stepstart, value = pair.first, pair.second
         if stop is not None and stepstart >= stop:
            if not values_only:
               yield prevstart, stop, prevval
            else:
               yield prevval
            return
         if not values_only:
            yield prevstart, stepstart, prevval
         else:
            yield prevval
         prevstart, prevval = stepstart, value
      else:
         if not values_only:
            if stop is None: 
               yield prevstart, self._swigobj.max_index+1, prevval
            else:
               yield prevstart, min( stop, self._swigobj.max_index+1), prevval
         else:
            yield prevval

   def __getitem__( self, index ):
      """Given a StepVector sv, writing sv[i] returns sv's element i (where i
      is an integer). 
      
      If you use a slice, i.e., 'sv[i:j]', you get a new StepVector which
      contains the elements from i to j. Beware that this performs a (shallow) 
      copy, i.e., if you need performance, using the 'get_steps' method may
      be preferable.
      
      Calling 'list' on a StepVector returned by 'sv[i:j]' is a convenient way
      to get a part of a Stepvector as an ordinary list, e.g., for inspection.      
      """
      if isinstance( index, slice ):
         if index.step is not None and index.step != 1:
             raise ValueError, "Striding slices (i.e., step != 1) are not supported"
         start = index.start if index.start is not None else self.start_index()
         stop = index.stop if index.stop is not None else ( start + len(self) )
         res = StepVector( stop - start, self._typecode, start )
         for stepstart, stepstop, value in self.get_steps( start, stop ):
            res[ stepstart : stepstop ] = value
         return res
      else:
         return self._swigobj.get_values_pystyle( index ).next().second
      
   def __iter__( self ):
      """When asked to provide an iterator, a StepVector will yield all its
      value, repeating each value according to the length of the step.
      Hence, calling, e.g., 'list( sv )' will transform the StepVector 'sv'
      into an ordinary list.
      """
      for start, stop, value in self.get_steps():
         for i in xrange( start, stop ):
            yield value
       
   def __repr__( self ):
      if self.start_index() == -sys.maxint - 1:
         start_s = "-inf"
      else:
         start_s = str( self.start_index() )
      if len(self) == sys.maxint:
         stop_s = "inf"
      else:
         stop_s = str( self.start_index() + len(self) )
      return "<%s object, type '%s', index range %s:%s, %d step(s)>" % (
         self.__class__.__name__, self.typecode(), start_s,
         stop_s, self.num_steps() )
       
   def typecode( self ):
      "Returns the typecode."
      return self._typecode
      
   def start_index( self ):
      "Returns the start index."
      return self._swigobj.min_index
      
   def __len__( self ):
      """The length of a StepVector is defined by its index range, not by
      the number of steps.
      """
      return self._swigobj.max_index - self._swigobj.min_index + 1
      
   def num_steps( self ):
      """Returns the number of steps, i.e., the number of triples that get_steps
      returns.
      """
      return self._swigobj.num_values()
      
   def __copy__( self ):
      return self[:]
      
   def __eq__( self, other ):
      """StepVectors can be compared for equality. This is conceptually done
      element for element, but, for performance, taking steps in one go.
      """
      if self.start_index() != other.start_index() or len(self) != len(other) or \
            self.typecode() != other.typecode():
         print "Mark A"
         return False
      selfsteps = self.get_steps()
      othrsteps = other.get_steps()
      selfstart, selfstop, selfval = selfsteps.next()
      othrstart, othrstop, othrval = othrsteps.next()
      while selfstop < self.start_index() + len(self) and \
            othrstop < other.start_index() + len(other):
         assert selfstart < othrstop and othrstart < selfstop
         if not( selfval == othrval ):
            return False
         if selfstop < othrstop:
            selfstart, selfstop, selfval = selfsteps.next()
         elif othrstop < selfstop:
            othrstart, othrstop, othrval = othrsteps.next()
         else:
            selfstart, selfstop, selfval = selfsteps.next()
            othrstart, othrstop, othrval = othrsteps.next()
      return True
      
   def __neq__( self, other ):
      return not ( self == other )
      
   def __reduce__( self ):
      return ( _StepVector_unpickle, ( len(self), self.typecode(),
         self.start_index(), list( self.get_steps() ) ) )
    
   def add_value( self, value, start=None, stop=None ):
      """Adds the given value to the vector elements given by start and stop.
      If they are missing, the beginning and/or end of the vector is taken.
      """
      self._swigobj.add_value( 
         start if start is not None else self._swigobj.min_index,
         stop-1 if stop is not None else self._swigobj.max_index,
         value )
    
   def apply( self, func, start = None, stop = None ):
      for stepstart, stepstop, value in self.get_steps( start, stop ):
         self[ stepstart : stepstop ] = func( value )
      
def _StepVector_unpickle( length, typecode, start_index, steps ):
   sv = StepVector( length, typecode, start_index )
   for start, stop, value in steps:
      sv[ start:stop ] = value
   return sv
      
%}
