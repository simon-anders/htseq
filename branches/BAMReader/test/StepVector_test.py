import unittest
import cPickle
import copy
from itertools import *

import StepVector


class MyTestCase( unittest.TestCase ):
   
   def checkPickleAndCopy( self, obj ):
      self.assertEqual( obj, copy.copy( obj ), "Copying failed" )
      self.assertEqual( obj, cPickle.loads( cPickle.dumps( obj ) ), "Pickling failed" )

   def checkPickle( self, obj ):
      self.assertEqual( obj, cPickle.loads( cPickle.dumps( obj ) ), "Pickling failed" )


class StepVector_Test( MyTestCase ):
   
   def test_float( self ):
      a = StepVector.StepVector( 100 )
      a[10:20] = 3.5
      a[40:75] = 2.75
      self.assertEqual(
         list( a.get_steps() ),
         [(0, 10, 0.0), (10, 20, 3.5), (20, 40, 0.0), (40, 75, 2.75), (75, 100, 0.0)] )
      self.assertEqual(
         list( a.get_steps( 15, 65 ) ),
         [(15, 20, 3.5), (20, 40, 0.0), (40, 65, 2.75) ] )
      self.assertEqual( a[19], 3.5 )
      self.assertEqual( a[20], 0 )
      self.assertEqual( 
         list( a.get_steps( 15, 65 ) ),
         list( a[ 15:65 ].get_steps() ) )
      self.assertEqual( 
         list( a )[ 55:79 ],
         list( a[ 55:79 ] ) )
      self.assertEqual( 
         list( a )[ 55: ],
         list( a[ 55: ] ) )

   def test_int( self ):
      a = StepVector.StepVector( 100, 'i', 5 )
      a[10:20] = 3
      a[40:75] = 2 
      self.assertEqual( 
         list( a.get_steps( 15, 65 ) ),
         list( a[ 15:65 ].get_steps() ) )
      self.assertEqual( 
         list( a )[ 50:74 ],
         list( a[ 55:79 ] ) )
      self.assertEqual( a.start_index(), 5 )
      self.assertEqual( len(a), 100 )
      self.assertEqual( a.num_steps(), 5 )
      self.checkPickleAndCopy( a )

   def test_cmp( self ):
      a = StepVector.StepVector( 100, 'i', 5 )
      a[10:20] = 3
      a[40:75] = 2 
      b = a[:]
      self.assertEqual( a, b )
      b[ 13 ] = 7
      self.assertNotEqual( a, b )
      self.assert_( a != b and not a == b)
      c = StepVector.StepVector( 100, 'i', 5 )
      c[10:15] = 3
      c[15:20] = 3
      c[40:75] = 2 
      self.assertEqual( a, c )
      self.assertEqual( a, a[:] )      

   def test_add( self ):
      a = StepVector.StepVector( 100, 'i', 5 )
      a[10:20] = 3
      a[40:75] = 2 
      b = a[:]
      b.add_value( 17, 20, 45 )
      b.add_value( 13, 25, 40 )
      self.assertNotEqual( a, b )
      b.add_value( -17, 20, 45 )
      b.add_value( -13, 25, 40 )
      self.assertEqual( a, b )
      self.assertRaises( IndexError, lambda: a.add_value( 10, 20, 17 ) )
      c = StepVector.StepVector( 100, 'O' )
      self.assertRaises( TypeError, lambda: c.add_value( 5, 10, 20 ) )

   def test_obj( self ):
      a = StepVector.StepVector( 100, 'O' )
      a[10:20] = "foo"
      a[40:75] = [ "bar", ( None, ), 17 ]
      a[45][2] = 19
      self.assertEqual( a[20], a[20] )
      self.assertEqual( 
         list( a.get_steps( 15, 65 ) ),
         list( a[ 15:65 ].get_steps() ) )
      self.assertEqual( 
         list( a )[ 50:74 ],
         list( a[ 50:74 ] ) )


unittest.main()      