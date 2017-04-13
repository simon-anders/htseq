import sys
import os
import glob
import distutils.util
import doctest

build_dir = "build/lib.%s-%s" % (distutils.util.get_platform(), sys.version[0:3])

sys.path.insert(0, os.path.join(os.getcwd(), build_dir))
import HTSeq


def test_rst_file(filename):
    print("Doctest of %s:" % os.path.basename(filename))
    os.chdir("example_data")
    (failure_count, test_count) = doctest.testfile(
            os.path.join("..", "doc", filename),
            module_relative=False)
    os.chdir("..")

    if failure_count == 0:
        print("All %d tests passed." % test_count)
        return True
    else:
        print("%d of %d tests failed." % (failure_count, test_count))
        return False


ok = True
if len(sys.argv) == 1:
    pathname = os.path.abspath(os.path.dirname(sys.argv[0]))
    rst_glob = os.path.join(pathname, '..', 'doc', '*.rst')
    print('RST files found in glob ', rst_glob+':', glob.glob(rst_glob))
    for fn in glob.glob(rst_glob):
        ok &= test_rst_file(os.path.basename(fn))
        print()
    if not ok:
        print("Not all tests passed.")
        exit(1)
elif len(sys.argv) == 2:
    test_rst_file(sys.argv[1])
else:
    print("Wrong usage")
    print("Call without arguments to run all doctest, or with the (base) name")
    print("of one rst file from the doc directory to run doctest on it.")
