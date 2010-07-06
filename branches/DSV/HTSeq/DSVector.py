import sys

import DynamicStepVector


def isiterable( obj ):
   try:
      iter( obj )
   except TypeError:
      return False
   return True


class DSVector:

   def __init__( self, typecode = 'i', start = 0, stop = sys.maxint, parent = None ):
      self.start = start
      self.stop = stop
      self.typecode = typecode
      if parent is None:
         if typecode == 'i':
            self._dsv = DynamicStepVector.intDSV()
         else:
            raise NotImplementedError, "So far, only type code 'i' works."
      else:
         self._dsv = parent._dsv
         
   def __getitem__( self, key ):

      if isinstance( key, int ):
         if key < self.start or key > self.stop:
            raise IndexError, "Index out of bounds"
         return self._dsv.get( key )

      elif isinstance( key, slice ):
         if key.start < self.start or key.stop > self.stop:
            raise IndexError, "Index out of bounds"
         if key.step is not None and key.step != 1:
            raise IndexError, "Slices with steps are not supported"
         return DSVector( self.typecode, key.start or self.start, 
            key.stop or self.stop, self )
         
      else:
         raise TypeError, "Illegal index"

   def __setitem__( self, key, value ):

      if isinstance( key, int ):
         if key < self.start or key > self.stop:
            raise IndexError, "Index out of bounds"
         next_val = self[key+1]
         self._dsv.set( key, value )
         self._dsv.set( key+1, next_val )

      elif isinstance( key, slice ):
         if key.start < self.start or key.stop > self.stop:
            raise IndexError, "Index out of bounds"
         if key.step is not None and key.step != 1:
            raise IndexError, "Slices with steps are not supported"
         key_start = key.start or self.start
         key_stop  = key.stop  or self.stop
         if key_stop - key_start > sys.maxint / 2 - 2:
            raise IndexError, "Attempt to assign to unlimited slice"

         if isinstance( value, self.__class__ ):
            if value.stop - value.start != key_stop - key_start:
               raise IndexError, "Attempt to assign a DSVector of wrong length"
            if self._dsv is value._dsv:
               # Copying unnecessary   
               return   
            for i in xrange( key_start, key_stop ):
               self[i] = value[i] # TODO: This is inefficient
         
         else:
            for i in xrange( key_start, key_stop ):
               self[i] = value # TODO: This is inefficient
         
      else:
         raise TypeError, "Illegal index type"
         
   def __iadd__( self, value ):
     self._dsv.add( self.start, self.stop, value )
     return self
     
   def __add__( self, value ):
      raise NotImplemented

   def __copy__( self ):
      raise NotImplemented
         
   def __iter__( self ):
      for i in xrange( self.start, self.stop ):
         yield self[i]
         
   def __len__( self ):
      return self.stop - self.start
      
   def __repr__( self ):
      return "<%s object, typecode '%s', range %s:%s>" % ( self.__class__.__name__,
         self.typecode, str(self.start) if self.start != -sys.maxint-1 else "-inf",
         str(self.stop) if self.stop != -sys.maxint else "inf" )
   
     
