import sys

import DynamicStepVector


def isiterable( obj ):
   try:
      iter( obj )
   except TypeError:
      return False
   return True


class DSVector:

   def __init__( self, typecode = 'i', start = 0, stop = sys.maxint, offset = 0, _dsv = None ):
      self.typecode = typecode
      self.start = start
      self.stop = stop
      self.offset = offset
      if _dsv is None:
         if typecode == 'i':
            self._dsv = DynamicStepVector.intDSV()
         else:
            raise NotImplementedError, "So far, only type code 'i' works."
      else:
         self._dsv = _dsv
      assert isinstance( self.typecode, str )
      assert isinstance( self.start, int )
      assert isinstance( self.stop, int )
      assert isinstance( self.offset, int )
      assert isinstance( self._dsv, DynamicStepVector.intDSV )
         
   def __getitem__( self, key ):

      if isinstance( key, int ):
         if key < self.start or key >= self.stop:
            raise IndexError, "Index out of bounds"
         return self._dsv.get( key + self.offset )

      elif isinstance( key, slice ):
         start = key.start if key.start is not None else self.start
         stop  = key.stop  if key.stop  is not None else self.stop 
         if start < self.start or stop > self.stop:
            raise IndexError, "Index out of bounds"
         if key.step is not None and key.step != 1:
            raise IndexError, "Slices with steps are not supported"
         return DSVector( self.typecode, 0, stop-start, start + self.offset, self._dsv )
         
      else:
         raise TypeError, "Illegal index"

   def __setitem__( self, key, value ):      

      if isinstance( key, int ):
         if key < self.start or key >= self.stop:
            raise IndexError, "Index out of bounds"
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

         if isinstance( value, self.__class__ ):
            if stop - start != len( value ):
               raise IndexError, "Attempt to assign a DSVector of wrong length"
            if self._dsv is value._dsv and start + d1.offset == value.start + value.offset:
               # Copying is unnecessarY:
               return   
            print "Debug notice: copying element-wise!"
            for i in xrange( start, stop ):
               self._dsv.set( i + self.offset, i + self.offset + 1, 
                  value._dsv.get( i - start + value.offset ) )
               # TODO: This is inefficient
         
         else:
            self._dsv.set( start + self.offset, stop + self.offset, value )
         
      else:
         raise TypeError, "Illegal index type"
         
   def __len__( self ):
      if self.start == -sys.maxint-1 or self.stop == sys.maxint:
         raise ValueError, "DSVector with unlimited range has undefined length."
      return self.stop - self.start
         
   def __iadd__( self, value ):
     self._dsv.add( self.start, self.stop, value )
     return self
     
   def __add__( self, value ):
      res = self.copy()
      res += value
      return res

   def copy( self ):
      # TODO: This is not compliant with the deepcopy protocol
      if self.typecode == "i":
         return DSVector( self.typecode, self.start, self.stop, self.offset, 
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
      valstart = self.start
      val = self[ self.start ]
      for i in xrange( self.start+1, self.stop ):
         if self[ i ] != val:
            yield( valstart, i, val )
            valstart = i
            val = self[ i ]
      yield( valstart, self.stop, val )

            
         
         
     
