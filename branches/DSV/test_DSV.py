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

it_f = d._dsv.get_step_iter(0, 20)
it_r = d._dsv.get_step_iter(0, 20, True) #reverse iterator
print "Iterating the steps form 0 to 20:"
while it_f.valid() : print it_f.next()
print "Now in reverse:"
while it_r.valid() : print it_r.prev()

print
print list( d[:30] )

print
for val in d[5:15].values_iter():
   print val

print
for start, stop, val in d[5:15].steps_iter():
   print start, stop, val
