import sys, os, glob, os.path
import distutils.util
build_dir = "build/lib.%s-%s" % ( distutils.util.get_platform(), sys.version[0:3] )
sys.path.insert( 0, os.path.join( os.getcwd(), build_dir ) )


import HTSeq.DSVector

d = HTSeq.DSVector.DSVector()

d[ 10:12 ] = 5
d[ 13 ] = 6 
d[ 14 ] = 3

for i in xrange( 5, 20 ): 
   print "%d:%d" % ( i, d[i] ),
print   

d2 = d[5:12] 
d2 += 100

print "Try A"
d[5:12] += 200

print "Try B"
d[5:12].__iadd__( 200 )

print "Try C"
d.__getitem__(slice(5,12)).__iadd__(200)

print list( d[:30] )
