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

#define TRY_CATCH_INDEX                                                    \
   try {                                                                   \
      $action                                                              \
   } catch (std::out_of_range &e) {                                        \
      SWIG_exception(SWIG_IndexError, e.what() );                          \
   } catch (type_error_non_arith &e) {                                     \
      SWIG_exception(SWIG_TypeError, "Illegal arithmetic operation" );     \
   }

%exception set_value { TRY_CATCH_INDEX }
%exception add_value { TRY_CATCH_INDEX }

template< class T >
class step_vector_for_python {
  public: 
   static const long int min_index;
   static const long int max_index;
   step_vector_for_python( );
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

   @classmethod
   def create( cls, length = sys.maxint, typecode = 'd', start_index = 0 ):
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
      obj = cls()
      obj._typecode = typecode
      obj._swigobj = swigclass( )    
      obj.start = start_index
      obj.stop = start_index + length
      return obj
     
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
      if isinstance( value, StepVector ):
         if self._swigobj is value._swigobj and \
               value.start == index.start and value.stop == index.stop:
            return
         else:
            raise NotImplemented, "Stepvector-to-Stepvector assignment still missing"
      if isinstance( index, slice ):
         if index.step is not None and index.step != 1:
             raise ValueError, "Striding slices (i.e., step != 1) are not supported"
         if index.start is None:
            start = self.start
         else:
            if index.start < self.start:
               raise IndexError, "start too small"
            start = index.start
         if index.stop is None:
            stop = self.stop
         else:
            if index.stop > self.stop:
               raise IndexError, "stop too large"
            stop = index.stop
         self._swigobj.set_value( start, stop-1, value )
         # Note the "-1": The C++ object uses closed intervals, but we follow
         # Python convention here and use half-open ones.
      else:
         self._swigobj.set_value( index, index, value )
    
   def get_steps( self, values_only = False, merge_steps=True ):
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
      startvals = self._swigobj.get_values_pystyle( self.start )
      prevstart = self.start
      prevval = startvals.next().second
      for pair in startvals:
         stepstart, value = pair.first, pair.second
         if merge_steps and value == prevval:
            continue
         if self.stop is not None and stepstart >= self.stop:
            if not values_only:
               yield prevstart, self.stop, prevval
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
            yield prevstart, min( self.stop, self._swigobj.max_index+1), prevval
         else:
            yield prevval

   def __getitem__( self, index ):
      """Given a StepVector sv, writing sv[i] returns sv's element i (where i
      is an integer). 
      
      If you use a slice, i.e., 'sv[i:j]', you get a view on the StepVector,
      i.e., the same data, but changed boundaries.
      """
      if isinstance( index, slice ):
         if index.step is not None and index.step != 1:
             raise ValueError, "Striding slices (i.e., step != 1) are not supported"
         if index.start is None:
            start = self.start
         else:
            if index.start < self.start:
               raise IndexError, "start too small"
            start = index.start
         if index.stop is None:
            stop = self.stop
         else:
            if index.stop > self.stop:
               raise IndexError, "stop too large"
            stop = index.stop
         res = self.__class__()
         res._typecode = self.typecode
         res._swigobj = self._swigobj
         res.start = start
         res.stop = stop
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
      if self.start == -sys.maxint - 1:
         start_s = "-inf"
      else:
         start_s = str( self.start )
      if self.stop == sys.maxint:
         stop_s = "inf"
      else:
         stop_s = str( self.stop )
      return "<%s object, type '%s', index range %s:%s, %d step(s)>" % (
         self.__class__.__name__, self.typecode(), start_s,
         stop_s, self.num_steps() )
       
   def typecode( self ):
      "Returns the typecode."
      return self._typecode
      
   def __len__( self ):
      """The length of a StepVector is defined by its index range, not by
      the number of steps.
      """
      return self.stop - self.start
      
   def num_steps( self ):
      """Returns the number of steps, i.e., the number of triples that get_steps
      returns.
      """
      return self._swigobj.num_values()
      
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
      if self.__class__ is not StepVector:
         raise NotImplemented, "Attempting to pickle a subclass of StepVector without redefined __reduce__."
      return ( 
         _StepVector_unpickle, 
         ( self.stop - self.start, self._typecode, self.start ),
         None,
         None,
         ( ( slice( start, stop ), val ) for start, stop, val in self.get_steps() ) )
    
   def __iadd__( self, value ):
      self._swigobj.add_value( self.start, self.stop-1, value )
      return self
    
   def apply( self, func, start = None, stop = None ):
      # TODO: check!
      for stepstart, stepstop, value in self.get_steps( start, stop ):
         self[ stepstart : stepstop ] = func( value )
      
def _StepVector_unpickle( length, typecode, start ):
   return StepVector.create( length, typecode, start )
    
%}
