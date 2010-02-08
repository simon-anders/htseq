import sys, os
import distutils.util
import doctest

build_dir = "build/lib.%s-%s" % ( distutils.util.get_platform(), sys.version[0:3] )

sys.path.insert( 0, os.path.join( os.getcwd(), build_dir ) )

print "Doctest of tutorial.rst:"
os.chdir( "example_data" )
(failure_count, test_count) = doctest.testfile( "../doc/tutorial.rst" )

if failure_count == 0:
   print "All %d tests passed." % test_count
else:   
   print "%d of %d tests failed." % (failure_count, test_count)
   exit( 1 )

