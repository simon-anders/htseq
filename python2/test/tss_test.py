import sys, os
import distutils.util
import numpy

build_dir = "build/lib.%s-%s" % ( distutils.util.get_platform(), sys.version[0:3] )

sys.path.insert( 0, os.path.join( os.getcwd(), build_dir ) )
import HTSeq
os.chdir( "example_data" )

print "TSS test, scheme 1"
execfile( os.path.join( "..", "doc", "tss1.py" ) )
profile1 = profile.copy()
print "finished"

print "TSS test, scheme 2"
execfile( os.path.join( "..", "doc", "tss2.py" ) )
profile2 = profile.copy()
print "finished",

if ( profile2 == profile1 ).all():
   print "and matches result from scheme 1."
else:   
   print "and differs to result from scheme 1!  <<<!!!>>>"

print "TSS test, scheme 3"
execfile( os.path.join( "..", "doc", "tss3.py" ) )
profile3 = profile.copy()
print "finished",

if ( profile3 == profile1 ).all():
   print "and matches result from scheme 1."
else:   
   print "and differs to result from scheme 1!  <<<!!!>>>"

