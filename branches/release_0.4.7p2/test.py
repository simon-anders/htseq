import sys, os, glob, os.path
import distutils.util
import doctest

build_dir = "build/lib.%s-%s" % ( distutils.util.get_platform(), sys.version[0:3] )

sys.path.insert( 0, os.path.join( os.getcwd(), build_dir ) )
os.chdir( "example_data" )

def test_rst_file( filename ):
   print "Doctest of %s:" % os.path.basename( filename )
   (failure_count, test_count) = doctest.testfile( filename )

   if failure_count == 0:
      print "All %d tests passed." % test_count
      return True
   else:   
      print "%d of %d tests failed." % (failure_count, test_count)
      return False

ok = True
if len(sys.argv) == 1:
   for fn in glob.glob( "../doc/*.rst" ):
      ok &= test_rst_file( fn )
      print
   if not ok:
      print "Not all tests passed."
      exit( 1 )
elif len(sys.argv) == 2:
   test_rst_file( "../doc/" + sys.argv[1] )
else:
   print "Wrong usage"
   print "Call without arguments to run all doctest, or with the (base) name"
   print "of one rst file from the doc directory to run doctest on it."
