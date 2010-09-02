# Attempt to wrap DynamicStepVector.hpp with Cython instead of SWIG

import sys

cdef extern from "Python.h":
   void Py_XINCREF( object o )

cdef extern from "<string>" namespace "std":
   cdef cppclass string:
      pass

cdef extern from "AutoPyObjPtr.h":
   cdef cppclass AutoPyObjPtr:
      AutoPyObjPtr( )
      AutoPyObjPtr( object )
      AutoPyObjPtr( AutoPyObjPtr )   
      object get_obj( )
   cdef void print_std_string( string s )

cdef extern from "DynamicStepVector.hpp":
   cdef cppclass DSV[ TKey, TValue ]:
      DSV( )
      void set( TKey, TKey, TValue )
      TValue get( TKey )
      string info( )
      
def test( ):
   cdef DSV[int,AutoPyObjPtr] dsv
   cdef int i
   print_std_string( dsv.info() )
   o = "A"
   #print ">", sys.getrefcount( o )
   dsv.set( 3, 5, AutoPyObjPtr( o ) )
   print_std_string( dsv.info() )
   for i in xrange( 0, 10 ):
      o2 = dsv.get( i ).get_obj()
      print i, o2
   print_std_string( dsv.info() )
   print "Now clearing ---------"
   dsv.set( 0, 100, AutoPyObjPtr( ) )
   print_std_string( dsv.info() )
   #print ">", sys.getrefcount( o )
   for i in xrange( 0, 10 ):
      o2 = dsv.get( i ).get_obj()
      print i, o2
