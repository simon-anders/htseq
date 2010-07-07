import sys, os, glob, os.path
import distutils.util
build_dir = "build/lib.%s-%s" % ( distutils.util.get_platform(), sys.version[0:3] )
sys.path.insert( 0, os.path.join( os.getcwd(), build_dir ) )


import HTSeq.DSVector

d = HTSeq.DSVector.DSVector()

for i in xrange(10, 20):
   d[i] = 100 + i
   
d[5:8] = 3   
   
print d
print list( d[:30] )

d2 = d[ 10:20 ]

print d2
print d2[0]
print d2[4]
print list( d2 )

d3 = d2[ 5:10 ]
print list( d3 )

print
print list( d[:30] )

print
for val in d[5:15].values_iter():
   print val

print
for start, stop, val in d[5:15].steps_iter():
   print start, stop, val
