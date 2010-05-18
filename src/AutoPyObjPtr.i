/* 
This file defines a smart pointer "AutoAutoPyObjPtr" for PyObjects that deals 
with Python reference counting. This allows to store Python objects in a C++ 
container.

Typemaps are provided to transform to standard Python objects.

Simon Anders, 2009-08-28
*/

%{

   struct AutoPyObjPtr {
      PyObject * obj;
      AutoPyObjPtr( PyObject * o = Py_None );
      AutoPyObjPtr( const AutoPyObjPtr & op );
      AutoPyObjPtr & operator= ( const AutoPyObjPtr & po ); 
      bool operator== ( const AutoPyObjPtr & po ) const; 
      ~AutoPyObjPtr( );
    #ifdef AUTOPYOBJPTR_EXTRAOPS
      AutoPyObjPtr & operator+=( const AutoPyObjPtr & po );   
      AutoPyObjPtr & operator+( const AutoPyObjPtr & po );   
    #endif
   };
  
   AutoPyObjPtr::AutoPyObjPtr( PyObject * o )
    : obj( o )
   {
      Py_XINCREF( obj );
   }
   
   AutoPyObjPtr::AutoPyObjPtr( const AutoPyObjPtr & op )
    : obj( op.obj )
   {
      Py_XINCREF( obj );
   }
   
   AutoPyObjPtr & AutoPyObjPtr::operator= ( const AutoPyObjPtr & po )
   {
      Py_XDECREF( obj );
      obj = po.obj;
      Py_XINCREF( obj );
      return *this;
   }

   bool AutoPyObjPtr::operator== ( const AutoPyObjPtr & po ) const
   {
      int res = PyObject_RichCompareBool( obj, po.obj, Py_EQ );
      assert( res == 0 || res == 1 );
      return res;
   }

   AutoPyObjPtr::~AutoPyObjPtr( )
   {
      Py_XDECREF( obj );
   }

 #ifdef AUTOPYOBJPTR_EXTRAOPS

   class type_error_non_arith {};

   AutoPyObjPtr & AutoPyObjPtr::operator+= ( const AutoPyObjPtr & po )
   { throw type_error_non_arith(); }

   AutoPyObjPtr & AutoPyObjPtr::operator+ ( const AutoPyObjPtr & po )
   { throw type_error_non_arith(); }

 #endif

%}

%typemap(in) AutoPyObjPtr {
   $1 = AutoPyObjPtr( $input );
}

%typemap(out) AutoPyObjPtr {
   Py_XINCREF( $1.obj );
   $result = $1.obj;
}
