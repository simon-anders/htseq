import numpy, random
import HTSeq

length = 10
start = 10

sv = HTSeq.StepVector.StepVector( length, 'd', start )
v = numpy.zeros( (length,), 'd' )

for i in xrange( 1000000 ):
   l = random.randrange( start, start+length )
   r = random.randrange( start, start+length )
   if not l<r:
      continue
   val = random.randrange( 3 )
   sv[l:r] = val
   v[l-start:r-start] = val
   for j in xrange(length):
      if( v[j] != sv[j+start] ):
         raise SystemError, "Failed!"
   if i % 10000 == 0:
      print v
      print list( sv.get_steps() )
