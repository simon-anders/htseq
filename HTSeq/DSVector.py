import sys, warnings

import DynamicStepVector


def isiterable( obj ):
   try:
      iter( obj )
   except TypeError:
      return False
   return True


class DSVector:

   def __init__( self, length = sys.maxint, typecode = 'd', start = 0, offset = 0, _dsv = None ):
      self.typecode = typecode
      self.start = start
      self.stop = start + length
      self.offset = offset
      if _dsv is None:
         if typecode == 'd':
            self._dsv = DynamicStepVector.floatDSV()
         elif typecode == 'i':
            self._dsv = DynamicStepVector.intDSV()
         elif typecode == 'b':
            raise NotImplementedError, "typecode 'b' does not work for now."
         elif typecode == 's':
            self._dsv = DynamicStepVector.strDSV()
         elif typecode == 'O':
            self._dsv = DynamicStepVector.pyDSV()
         else:
            raise ValueError, "Illegal typecode"
      else:
         if typecode not in ( 'd', 'i', 'b', 's', 'O' ):
            raise ValueError, "Illegal typecode"
         if (
               ( typecode == 'd' and not isinstance( _dsv, DynamicStepVector.floatDSV ) ) or
               ( typecode == 'i' and not isinstance( _dsv, DynamicStepVector.intDSV ) ) or
            #  ( typecode == 'b' and not isinstance( _dsv, DynamicStepVector.boolDSV ) ) or
               ( typecode == 's' and not isinstance( _dsv, DynamicStepVector.strDSV ) ) or
               ( typecode == 'O' and not isinstance( _dsv, DynamicStepVector.pyDSV ) ) ):
            raise ValueError, "typecode does not match _dsv"
         if typecode == 'b':
            raise NotImplementedError, "typecode 'b' does not work for now."
         self._dsv = _dsv
      assert isinstance( self.typecode, str )
      assert isinstance( self.start, int )
      assert isinstance( self.stop, int )
      assert isinstance( self.offset, int )
         
   def __getitem__( self, key ):

      if isinstance( key, int ):
         if key < self.start or key >= self.stop:
            raise IndexError, "Index out of bounds"
         return self._dsv.get( key + self.offset )

      elif isinstance( key, slice ):
         if key.step is not None and key.step != 1:
            raise IndexError, "Slices with steps are not supported"
         return self.slice( key.start, key.stop, False )
         
      else:
         raise TypeError, "Illegal index"

   def slice( self, start, stop, withOffset = True ):
      if start is None:
         start = self.start
      if stop is None:
         stop = self.stop
      if start < self.start or stop > self.stop:
         raise IndexError, "Index out of bounds"
      if withOffset:
         return DSVector( stop-start, self.typecode, 0, start + self.offset, self._dsv )
      else:
         return DSVector( stop-start, self.typecode, start, self.offset, self._dsv )
   

   def __setitem__( self, key, value ):      

      if isinstance( key, int ):
         if key < self.start or key >= self.stop:
            raise IndexError, "Index out of bounds"
         self._check_value_type( value )
         self._dsv.set( key + self.offset, key + self.offset + 1, value )

      elif isinstance( key, slice ):         
         start = key.start if key.start is not None else self.start
         stop  = key.stop  if key.stop  is not None else self.stop 
         if start < self.start or stop > self.stop:
            raise IndexError, "Index out of bounds"
         if key.step is not None and key.step != 1:
            raise IndexError, "Slices with steps are not supported"
         if stop - start > sys.maxint / 2 - 2:
            raise IndexError, "Attempt to assign to unlimited slice"

         if isinstance( value, DSVector ):
            if stop - start != len( value ):
               raise IndexError, "Attempt to assign a DSVector of wrong length"
            if self._dsv is value._dsv and start + self.offset == value.start + value.offset:
               # Copying is unnecessary:
               return   
            if self.typecode != value.typecode:
               raise TypeError, "Type mismatch in copying between DSVector objects."
            print "Debug notice: copying element-wise!"
            for i in xrange( start, stop ):
               self._dsv.set( i + self.offset, i + self.offset + 1, 
                  value._dsv.get( i - start + value.offset ) )
               # TODO: This is inefficient
         
         else:
            self._check_value_type( value )
            self._dsv.set( start + self.offset, stop + self.offset, value )
         
      else:
         raise TypeError, "Illegal index type"
         
   def __len__( self ):
      if self.start == -sys.maxint-1 or self.stop == sys.maxint:
         raise ValueError, "DSVector with unlimited range has undefined length."
      return self.stop - self.start
         
   def __iadd__( self, value ):
     self._dsv.add( self.start + self.offset, self.stop + self.offset, value )
     return self
     
   def __add__( self, value ):
      res = self.copy()
      res += value
      return res

   def copy( self ):
      # TODO: This is not compliant with the deepcopy protocol
      if self.typecode == "i":
         return DSVector( self.stop-self.start, self.typecode, self.start, self.offset, 
            DynamicStepVector.intDSV( self._dsv ) )
      else:
         raise NotImplemented
         
   def __copy__( self ):
      return self.copy( )
         
   def __iter__( self ):
      for i in xrange( self.start, self.stop ):
         yield self[i]
         
   def __len__( self ):
      return self.stop - self.start
      
   def __repr__( self ):
      return "<%s object, typecode '%s', range %s:%s%s>" % ( 
         self.__class__.__name__, self.typecode, 
         str(self.start) if self.start != -sys.maxint-1 else "-inf",
         str(self.stop) if self.stop != sys.maxint else "inf",
         (", offset %d" % self.offset) if self.offset != 0 else "")
   
   def values_iter( self ):
      return iter(self)
      
   def steps_iter( self ):
      # TO DO: Make this efficient
      valstart = self.start
      val = self[ self.start ]
      for i in xrange( self.start+1, self.stop ):
         if self[ i ] != val:
            yield( valstart, i, val )
            valstart = i
            val = self[ i ]
      yield( valstart, self.stop, val )

   def get_steps( self, start=None, stop=None, values_only = False ):
      warnings.warn( "'get_steps' is deprecated, use 'steps_iter' instead.",
         DeprecationWarning )
      if not values_only:
         return self[start:stop].steps_iter()
      else:
         return ( value for start_, stop_, value in self[start:stop].steps_iter() )
      
   def add_value( self, value, start=None, stop=None ):
      warnings.warn( "'add_value' is deprecated, use '+=' instead.",
         DeprecationWarning )
      self[start:stop] += value
         
   def _check_value_type( self, value ):
      if (
            ( self.typecode == 'd' and not ( type( value ) == float or type( value ) == int ) ) or
            ( self.typecode == 'i' and not type( value ) == float ) or
            ( self.typecode == 'b' and not type( value ) == bool ) or
            ( self.typecode == 's' and not type( value ) == str ) ):
         raise TypeError, "Type mismatch in assignment to step vector."
      
         
     
